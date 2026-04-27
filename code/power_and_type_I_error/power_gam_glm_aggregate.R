#' Aggregate files for GAM, NB GLM, and Poisson GLM power results on both Poisson-Leroux and Poisson-spNNGP generated data
#' @author Florica Constantine

# Libraries
library(dplyr)

# Paths
path <- "/scratch/users/spatialseq/natgen_kidney/"
data_path <- file.path(path, "glm_offset_power_output_onlyinteract")
processed_data_path <- file.path(path, "processed")

# List of files to read in
pattern_type <- ".rds"
file_list <- base::list.files(data_path, pattern = pattern_type)

# Aggregate files
if (file.exists(file.path(
  processed_data_path,
  "glm_offset_power_perfsummary_onlyinteract.rds"
))
&
file.exists(file.path(
  processed_data_path,
  "glm_offset_power_waldstats_onlyinteract.rds"
))) {
  summary_out_all <- readRDS(file.path(
    processed_data_path,
    "glm_offset_power_perfsummary_onlyinteract.rds"
  ))
  wald_df_all <- readRDS(file.path(
    processed_data_path,
    "glm_offset_power_waldstats_onlyinteract.rds"
  ))
} else {
  summary_out_list <- list()
  wald_out_list <- list()
  for (idx in 1:length(file_list)) {
    if (0 == (idx %% 1000)) {
      cat(idx, "\n")
    }
    
    # Read in output
    out <- readRDS(paste0(data_path, file_list[idx]))
    
    summary_out_list[[idx]] <- out$summary_out_all
    wald_out_list[[idx]] <- out$wald_df
  }
  summary_out_all <- dplyr::bind_rows(summary_out_list)
  saveRDS(
    summary_out_all,
    paste0(
      processed_data_path,
      "glm_offset_power_perfsummary_onlyinteract.rds"
    ),
    compress = "xz"
  )
  
  wald_df_all <- dplyr::bind_rows(wald_out_list)
  wald_df_all$wald_pval <- stats::pchisq(wald_df_all$wald_stat_t^2, 1, lower.tail = FALSE)
  saveRDS(
    wald_df_all,
    file.path(
      processed_data_path,
      "glm_offset_power_waldstats_onlyinteract.rds"
    ),
    compress = "xz"
  )
}
