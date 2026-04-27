#' @description
#' This script computes p-values for the synthetic simulations using the real data fits.
#' Elsewhere, we have used the TESSERA procedure to find a cutoff and fit a scaled
#' chi_1^2(c) distribution. We use those parameters/that distribution to get
#' p-values.
#' @author Florica Constantine


## Libraries

# General Utilities
library(dplyr)
library(tibble)


## Paths

path <-  "/scratch/users/spatialseq/natgen_kidney/"
processed_data_path <- file.path(path, "processed")

# Load in synthetic power data
# Results from fitting gene IGHG2 at different effect sizes across 100 trials
tessera_wald_df <- readRDS(file.path(processed_data_path, "power_wald_df.rds"))
tessera_misspec_wald_df <- readRDS(
  file.path(
    processed_data_path,
    "power_misspec_results_waldstats_Leroux_onlyinteract.rds"
  )
)
tessera_wald_df <- dplyr::bind_rows(tessera_wald_df, tessera_misspec_wald_df)
# Subset to only TESSERA Poisson-Leroux fits
tessera_wald_df <- tessera_wald_df %>% dplyr::filter(fit_model %in% c("Leroux"))
rm(tessera_misspec_wald_df)

# Load in real data
nocutoff_res <- readRDS(
  file.path(
    processed_data_path,
    "real_hyptest_tessera.rds"
  )
)


## Parameters

PVAL_THRESH <- 0.05


## Process data

# Create celltype column
baseline_group <- sapply(as.character(tessera_wald_df$baseline_contrast), function (x) {
  gsub("Group", "", strsplit(x, ":")[[1]][1])
})
changed_group <- sapply(as.character(tessera_wald_df$changed_covariate), function (x) {
  gsub("Group", "", strsplit(x, ":")[[1]][1])
})
celltype_name <- sapply(as.character(tessera_wald_df$changed_covariate), function (x) {
  gsub("celltype", "", strsplit(x, ":")[[1]][2])
})
contrast_description <- paste0(celltype_name, ":", baseline_group, "-", changed_group)
tessera_wald_df$celltype <- celltype_name
tessera_wald_df$contrast_description <- contrast_description


## Extract parameters

# Data frame with chi_1^2 fits for each contrast, as well as threshold parameters
parameter_df <- distinct(nocutoff_res$wald_df_all[, c(
  "contrast_indices",
  "contrast_description",
  "ncchi2_scale",
  "ncchi2_shift",
  "wald_thresh",
  "tessera_pi0"
)])

# Merge the parameters into the synthetic dataframe
tessera_wald_df <- merge(tessera_wald_df,
                         parameter_df,
                         by = c("contrast_indices", "contrast_description"))


## Compute p-values

tessera_wald_df$wald_pval <- pchisq(
  tessera_wald_df$wald_stat_t^2 / tessera_wald_df$ncchi2_scale,
  df = 1,
  ncp = tessera_wald_df$ncchi2_shift,
  lower.tail = FALSE
)


## Compute power and Type-I error

groupby_columns <- c(
  "gene",
  "datagen_model",
  "fit_model",
  "celltype",
  "contrast_description",
  "contrast_indices",
  "true_contrast_value"
)
power_df <- (
  tessera_wald_df
  %>% dplyr::group_by(dplyr::across(dplyr::all_of(groupby_columns)))
  %>% dplyr::summarise(
    p_rej = mean(wald_pval < PVAL_THRESH, na.rm = TRUE),
    .groups = "drop_last"
  )
)

p_value_thresholds <- seq(0.0, 1.0, 0.01)
roc_df <- list()
for (thresh in p_value_thresholds) {
  # Compute FPR: true contrast == 0
  fpr_inner <- (
    tessera_wald_df %>% dplyr::filter(true_contrast_value == 0)
    %>% dplyr::group_by(dplyr::across(dplyr::all_of(groupby_columns)))
    %>% dplyr::summarise(
      fpr = mean(wald_pval < thresh, na.rm = TRUE),
      .groups = "drop_last"
    )
  )
  # Compute TPR: true contrast != 0
  tpr_inner <- (
    tessera_wald_df %>% dplyr::filter(true_contrast_value != 0)
    %>% dplyr::group_by(dplyr::across(dplyr::all_of(groupby_columns)))
    %>% dplyr::summarise(
      tpr = mean(wald_pval < thresh, na.rm = TRUE),
      .groups = "drop_last"
    )
  )
  # Combine into one dataframe
  df_inner <- merge(fpr_inner, tpr_inner, by = groupby_columns[!grepl("true_contrast_value", groupby_columns)])
  # Add in threshold information
  df_inner <- tibble::add_column(df_inner,
                                 threshold = rep(thresh, nrow(df_inner)),
                                 .before = "fpr")
  roc_df[[1 + length(roc_df)]] <- df_inner
}
roc_df <- dplyr::bind_rows(roc_df)


## Save

# Create list of objects to save
out_list <- list(
  tessera_wald_df = tessera_wald_df,
  power_df = power_df,
  roc_df = roc_df
)

saveRDS(
  out_list,
  file.path(
    processed_data_path,
    "synth_from_real_hyptest_tessera.rds"
  )
)
