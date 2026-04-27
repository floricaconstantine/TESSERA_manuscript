#' Calculates p-values for between-cell-type contrasts on kidney data using the empirical null estimation procedure from the TESSERA package
#' @author Florica Constantine


## Libraries

# General Utilities
library(dplyr)
# Custom
library(TESSERA)


## Paths
path <-  "/scratch/users/spatialseq/natgen_kidney/"
processed_data_path <- file.path(path, "processed")

## Set up
contrast_setup <- "all"
fit_model_TESSERA <- "Leroux"

## Read in model results
wald_df_all <- readRDS(file.path(processed_data_path, "real_results_waldstats_Leroux_between_celltypes_onlyinteract.rds"))

## Parameters

PVAL_ADJUST_METHOD <- "BH"
PVAL_THRESH <- 0.05
# For threshold selection
quantile_spacing <- 0.01
quantile_list <- seq(0 + quantile_spacing, 1.0 - quantile_spacing, quantile_spacing)

## Setup TESSERA data

# Subset fit model
wald_df_all <- wald_df_all %>% dplyr::filter(fit_model == fit_model_TESSERA)


## Find optimal threshold

# Note: The below code reproduces what TESSERA::select_Wald_threshold does and was written before that function existed

# Function to get errors and fits at a given threshold value
wrapper_function <- function(wald_thresh, wald_data) {
  # Subset to statistics under threshold
  wald_nulls <- wald_data[wald_data < wald_thresh]

  # Fit distribution (conditionally) to statistics under threshold
  chi2_fit <- TESSERA::fit_scaled_noncentral_chi2(wald_nulls, wald_thresh)
  # Just in case
  if (any(is.na(chi2_fit)))
    return(NULL)

  # Compute p-values
  # We use the fitted scaling (chi2_fit[1]) and ncp (chi2_fit[2])
  wald_null_pvals <- stats::pchisq(
    wald_nulls / chi2_fit[1],
    df = 1,
    ncp = chi2_fit[2],
    lower.tail = FALSE
  )

  # Prepare quantiles for comparison
  obs <- sort(wald_null_pvals)
  theo <- stats::ppoints(length(obs))

  # Logged versions
  eps <- 1e-10 # Small constant to prevent log(0)
  obs_log  <- -log10(obs + eps)
  theo_log <- -log10(theo + eps)

  # Helper function to calculate error metrics
  calc_errors <- function(o, t) {
    err <- abs(o - t)
    list(
      MSE   = mean(err^2),
      MAE   = mean(err),
      MedAE = median(err),
      MaxAE = max(err)
    )
  }
  # Calculate for both scales
  raw_errors <- calc_errors(obs, theo)
  log_errors <- calc_errors(obs_log, theo_log)

  # Combine into a single dataframe
  res <- data.frame(
    threshold = wald_thresh,
    scale     = as.numeric(chi2_fit[1]),
    shift     = as.numeric(chi2_fit[2]),
    n_obs     = length(wald_nulls)
  )
  # Append metrics with prefixes
  res[, paste0("Raw_", names(raw_errors))] <- as.list(raw_errors)
  res[, paste0("Log_", names(log_errors))] <- as.list(log_errors)

  return(res)
}


## Sweep over thresholds for each data model

if ("per" == contrast_setup) {
  fit_results_df_list <- list()
  for (ci in unique(wald_df_all$contrast_indices)) {
    print(ci)
    local_idx_tessera <- which(wald_df_all$contrast_indices == ci)

    # Find threshold and fit
    wald_F <- wald_df_all$wald_stat_t[local_idx_tessera]^2
    wald_F <- wald_F[!is.na(wald_F)]
    threshold_grid  <- stats::quantile(wald_F, quantile_list, na.rm = TRUE)

    # Record start time
    t0 <- Sys.time()
    # Standard Loop across Thresholds
    fit_results_list <- list()
    for (t_idx in 1:length(threshold_grid)) {
      cat(ci, t_idx, difftime(Sys.time(), t0, units = "sec"), "\n")
      fit_results_list[[t_idx]] <- wrapper_function(threshold_grid[t_idx], wald_F)
    }
    # Combine the list into one dataframe
    fit_results_df <- dplyr::bind_rows(fit_results_list)
    t1 <- Sys.time()

    # Clean up and Combine
    fit_results_df$fit_model <- fit_model_TESSERA
    fit_results_df$quantile <- quantile_list
    fit_results_df$contrast_indices <- ci

    fit_results_df_list[[ci]] <- fit_results_df
  }
  fit_results_df <- dplyr::bind_rows(fit_results_df_list)
} else if ("all" == contrast_setup) {
  # Find threshold and fit
  wald_F <- wald_df_all$wald_stat_t^2
  wald_F <- wald_F[!is.na(wald_F)]
  threshold_grid   <- stats::quantile(wald_F, quantile_list, na.rm = TRUE)

  # Record start time
  t0 <- Sys.time()
  # Standard Loop across Thresholds
  fit_results_list <- list()
  for (t_idx in 1:length(threshold_grid)) {
    cat(t_idx, difftime(Sys.time(), t0, units = "sec"), "\n")
    fit_results_list[[t_idx]] <- wrapper_function(threshold_grid[t_idx], wald_F)
  }
  # Combine the list into one dataframe
  fit_results_df <- dplyr::bind_rows(fit_results_list)
  t1 <- Sys.time()

  # Clean up and Combine
  fit_results_df$fit_model <- fit_model_TESSERA
  fit_results_df$quantile <- quantile_list
}


## Compute p-values

wald_df_all$wald_pval <- NA
wald_df_all$ncchi2_scale <- NA
wald_df_all$ncchi2_shift <- NA
wald_df_all$wald_thresh <- NA
wald_df_all$tessera_pi0 <- NA
if ("per" == contrast_setup) {
  for (ci in unique(wald_df_all$contrast_indices)) {
    local_idx_tessera <- which(wald_df_all$contrast_indices == ci)

    subset_fit <- fit_results_df[fit_results_df$contrast_indices == ci, ]
    min_idx <- which.min(subset_fit$Raw_MSE)
    # Extract fit
    wald_thresh <- subset_fit$threshold[min_idx]
    chi2_fit <- c(subset_fit$scale[min_idx], subset_fit$shift[min_idx])

    # Store results
    wald_df_all$ncchi2_scale[local_idx_tessera] <- chi2_fit[1]
    wald_df_all$ncchi2_shift[local_idx_tessera] <- chi2_fit[2]
    wald_df_all$wald_thresh[local_idx_tessera] <- wald_thresh
    wald_df_all$tessera_pi0[local_idx_tessera] <- mean(wald_df_all$wald_stat_t[local_idx_tessera]^2 < wald_thresh, na.rm = TRUE)

    # Compute p-values (unconditionally)
    wald_df_all$wald_pval[local_idx_tessera] <- stats::pchisq(
      wald_df_all$wald_stat_t[local_idx_tessera]^2 / chi2_fit[1],
      1,
      ncp = chi2_fit[2],
      lower.tail = FALSE
    )
  }
} else if ("all" == contrast_setup) {
  local_idx_tessera <- which(!is.na(wald_df_all$wald_stat_t))

  min_idx <- which.min(fit_results_df$Raw_MSE)
  # Extract fit
  wald_thresh <- fit_results_df$threshold[min_idx]
  chi2_fit <- c(fit_results_df$scale[min_idx], fit_results_df$shift[min_idx])

  # Store results
  wald_df_all$ncchi2_scale <- chi2_fit[1]
  wald_df_all$ncchi2_shift <- chi2_fit[2]
  wald_df_all$wald_thresh <- wald_thresh
  wald_df_all$tessera_pi0 <- mean(wald_df_all$wald_stat_t[local_idx_tessera]^2 < wald_thresh, na.rm = TRUE)

  # Compute p-values (unconditionally)
  wald_df_all$wald_pval <- stats::pchisq(wald_df_all$wald_stat_t^2 / chi2_fit[1],
                                         1,
                                         ncp = chi2_fit[2],
                                         lower.tail = FALSE)
}


## Process results

# Adjust p-values
wald_df_all$wald_pval_adj <- NA
for (cd in unique(wald_df_all$contrast_description)) {
  local_idx <- wald_df_all$contrast_description == cd
  wald_df_all$wald_pval_adj[local_idx] <- p.adjust(wald_df_all$wald_pval[local_idx], method =
                                                     PVAL_ADJUST_METHOD)
}

# Count number of significant genes per contrast
id_cols <- c("fit_model",
             "kernel_type",
             "contrast_indices",
             "contrast_description")
wald_gene_nsig_df <- data.frame(
  wald_df_all
  %>% dplyr::group_by(dplyr::across(dplyr::all_of(id_cols)))
  %>% dplyr::summarise(n_sig_genes = sum(wald_pval_adj < PVAL_THRESH, na.rm =
                                           TRUE))
)


## QQ for p-values

create_qq_df <- function(df, contrast_setup = "all") {
  if (contrast_setup == "per") {
    # Process each contrast index separately
    qq_df <- df %>%
      dplyr::group_by(contrast_indices) %>%
      dplyr::arrange(wald_pval, .by_group = TRUE) %>%
      dplyr::mutate(
        # Standard Uniform Quantiles: (i - 0.5) / n
        theoretical = (dplyr::row_number() - 0.5) / dplyr::n(),
        observed = wald_pval
      ) %>%
      dplyr::ungroup()

  } else {
    # Process all data as a single distribution
    qq_df <- df %>%
      dplyr::arrange(wald_pval) %>%
      dplyr::mutate(
        theoretical = (dplyr::row_number() - 0.5) / dplyr::n(),
        observed = wald_pval
      )
  }

  return(qq_df)
}
pval_qq_df <- create_qq_df(wald_df_all[wald_df_all$wald_stat_t^2 < wald_df_all$wald_thresh, ], contrast_setup)


## Save

# Create list of objects to save
out_list <- list(
  wald_df_all = wald_df_all,
  wald_gene_nsig_df = wald_gene_nsig_df,
  pval_qq_df = pval_qq_df,
  fit_results_df = fit_results_df
)

saveRDS(
  out_list,
  file.path(
    processed_data_path,
    "real_hyptest_tessera_between_celltypes.rds"
  )
)
