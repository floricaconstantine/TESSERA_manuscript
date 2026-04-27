#' @description
#' This script computes p-values for the synthetic simulations using a chi_1^2 distribution.
#' @author Florica Constantine


## Libraries

# General Utilities
library(dplyr)
# Custom
library(TESSERA)


## Paths

path <-  "/scratch/users/spatialseq/natgen_kidney/"
processed_data_path <- file.path(path, "processed")

# Results from fitting gene IGHG2 at different effect sizes across 100 trials
tessera_wald_df <- readRDS(paste0(processed_data_path, "power_wald_df.rds"))
tessera_misspec_wald_df <- readRDS(
  paste0(
    processed_data_path,
    "power_misspec_results_waldstats_Leroux_onlyinteract.rds"
  )
)
tessera_wald_df <- dplyr::bind_rows(tessera_wald_df, tessera_misspec_wald_df)
# Subset to only TESSERA fits
tessera_wald_df <- tessera_wald_df %>% dplyr::filter(fit_model %in% c("Leroux"))
rm(tessera_misspec_wald_df)


## Parameters

PVAL_THRESH <- 0.05
quantile_spacing <- 0.01
quantile_list <- seq(0 + quantile_spacing, 1.0 - quantile_spacing, quantile_spacing)

fit_model_TESSERA <- "Leroux"


## Setup TESSERA data

# Subset fit model
tessera_wald_df <- tessera_wald_df %>% dplyr::filter(fit_model == fit_model_TESSERA)
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


## Compute p-values

tessera_wald_df$wald_pval <- pchisq(
  tessera_wald_df$wald_stat_t^2,
  df = 1,
  ncp = 0,
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
out_list <- list(roc_df = roc_df,
                 power_df = power_df,
                 tessera_wald_df = tessera_wald_df)

saveRDS(
  out_list,
  paste0(
    processed_data_path,
    "synth_from_real_hyptest_theoretical.rds"
  )
)
