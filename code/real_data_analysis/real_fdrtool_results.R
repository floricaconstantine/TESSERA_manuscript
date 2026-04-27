#' Calculates p-values for within-cell-type, between-condition contrasts on kidney data using the FNDR empirical null estimation procedure, as implemented in the fdrtool package
#' @author Florica Constantine


## Libraries

# General Utilities
library(dplyr)
# Custom
library(TESSERA)
# P-values
library(fdrtool)


## Paths

path <-  "/scratch/users/spatialseq/natgen_kidney/"
processed_data_path <- file.path(path, "processed")

# Read in model results
wald_df_all <- readRDS(
  file.path(
    processed_data_path,
    "real_results_waldstats_Leroux_onlyinteract.rds"
  )
)
wald_df_all <- dplyr::select(wald_df_all, -c("contrast_string"))


## Parameters

PVAL_ADJUST_METHOD <- "BH"
PVAL_THRESH <- 0.05

contrast_setup_list <- c("all")
fit_model_list <- c("Leroux")
sim_df <- expand.grid(list(contrast_setup = contrast_setup_list, fit_model = fit_model_list))


## Clusterize

run_idx <- NA
offset <- 0
args <- commandArgs(trailingOnly = TRUE)
if (0 == length(args)) {
  stop("No arguments supplied: ERROR")
} else {
  for (idx in 1:length(args)) {
    eval(parse(text = args[idx]))
  }
  
  print(run_idx)
  print(offset)
  run_idx <- run_idx + offset
  print(run_idx)
  
  # run_idx must be defined in arguments
  if ((1 > run_idx) || (nrow(sim_df) < run_idx)) {
    stop("run_idx must be between 1 and # of possible jobs")
  }
}

contrast_setup <- sim_df$contrast_setup[run_idx]
fit_model_TESSERA <- sim_df$fit_model[run_idx]
print(sim_df[run_idx, ])


## Setup TESSERA data

# Subset to results from fitted model of interest
wald_df_all <- wald_df_all %>% dplyr::filter(fit_model == fit_model_TESSERA)


## Empirical Null Estimation (fdrtool)

wald_df_all$wald_pval <- NA
wald_df_all$locfdr_sd <- NA
wald_df_all$locfdr_pi0 <- NA
if ("per" == contrast_setup) {
  for (ci in unique(wald_df_all$contrast_indices)) {
    local_idx_tessera <- which(wald_df_all$contrast_indices == ci &
                                 !is.na(wald_df_all$wald_stat_t))
    
    # Run FNDR method per contrast
    locfdr_res <- fdrtool::fdrtool(
      wald_df_all$wald_stat_t[local_idx_tessera],
      statistic = "normal",
      cutoff.method = "fndr",
      plot = FALSE,
      verbose = FALSE
    )
    
    wald_df_all$locfdr_sd[local_idx_tessera] <- as.numeric(locfdr_res$param[1, "sd"])
    wald_df_all$locfdr_pi0[local_idx_tessera] <- as.numeric(locfdr_res$param[1, "eta0"])
    wald_df_all$wald_pval[local_idx_tessera] <- locfdr_res$pval
  }
} else if ("all" == contrast_setup) {
  local_idx_tessera <- which(!is.na(wald_df_all$wald_stat_t))
  
  # Run FNDR method across all contrasts
  locfdr_res <- fdrtool::fdrtool(
    wald_df_all$wald_stat_t[local_idx_tessera],
    statistic = "normal",
    cutoff.method = "fndr",
    plot = FALSE,
    verbose = FALSE
  )
  
  wald_df_all$locfdr_sd[local_idx_tessera] <- as.numeric(locfdr_res$param[1, "sd"])
  wald_df_all$locfdr_pi0[local_idx_tessera] <- as.numeric(locfdr_res$param[1, "eta0"])
  wald_df_all$wald_pval[local_idx_tessera] <- locfdr_res$pval
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
pval_qq_df <- create_qq_df(wald_df_all, contrast_setup)


## Save

# Create list of objects to save
out_list <- list(
  wald_df_all = wald_df_all,
  run_settings = sim_df[run_idx, ],
  wald_gene_nsig_df = wald_gene_nsig_df,
  pval_qq_df = pval_qq_df
)


saveRDS(out_list,
        file.path(output_data_path, "real_hyptest_fdrtool_fndr.rds"))
