#' Calculates p-values for within-cell-type between-condition contrasts on kidney data using the theoretical chi_1^2 distribution
#' @author Florica Constantine


## Libraries

# General Utilities
library(dplyr)
# Custom
library(TESSERA)


## Paths

path <- "/scratch/users/spatialseq/natgen_kidney/"
processed_data_path <- file.path(path, "processed")

# Read in model results
wald_df_all <- readRDS(file.path(
  processed_data_path,
  "real_results_waldstats_Leroux_onlyinteract.rds"
))
wald_df_all <- dplyr::select(wald_df_all, -c("contrast_string"))


## Parameters

PVAL_ADJUST_METHOD <- "BH"
PVAL_THRESH <- 0.05

fit_model_list <- c("Leroux")
sim_df <- expand.grid(
  list(
    fit_model = fit_model_list
  )
)


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

fit_model_TESSERA <- sim_df$fit_model[run_idx]
print(sim_df[run_idx, ])


## Setup TESSERA data

# Subset fit model
wald_df_all <- wald_df_all %>% dplyr::filter(fit_model == fit_model_TESSERA)


# P-values
wald_df_all$wald_pval <- pchisq(wald_df_all$wald_stat_t^2, df = 1, ncp = 0, lower.tail = FALSE)


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
pval_qq_df <- create_qq_df(wald_df_all, "all")


## Save

# Create list of objects to save
out_list <- list(
  wald_df_all = wald_df_all,
  run_settings = sim_df[run_idx, ],
  wald_gene_nsig_df = wald_gene_nsig_df,
  pval_qq_df = pval_qq_df
)


saveRDS(
  out_list,
  file.path(
    processed_data_path,
    "real_hyptest_theoretical_chisq.rds"
  )
)
