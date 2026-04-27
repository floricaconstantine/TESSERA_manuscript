#' Aggregate files for TESSERA power results on both Poisson-Leroux and Poisson-spNNGP generated data
#' @author Florica Constantine

# Libraries
library(dplyr)

# Paths
path <- "/scratch/users/spatialseq/natgen_kidney/"
modelfit_data_path  <- file.path(path, "power_output_onlyinteract")
processed_data_path <- file.path(path, "processed")

# List of files to read in
pattern_type <- "fitmodel_Leroux"
file_list <- base::list.files(modelfit_data_path, pattern = pattern_type)

# Aggregate files
if (file.exists(file.path(processed_data_path, "power_perf_df.rds"))
    &
    file.exists(file.path(processed_data_path, "power_wald_df.rds"))) {
  summary_out_all <- readRDS(file.path(processed_data_path, "power_perf_df.rds"))
  wald_df_all <- readRDS(file.path(processed_data_path, "power_wald_df.rds"))
} else{
  summary_out_list <- list()
  wald_out_list <- list()
  for (idx in 1:length(file_list)) {
    if (0 == (idx %% 250)) {
      cat(idx, "\n")
    }
    
    out <- readRDS(file.path(modelfit_data_path, file_list[idx]))
    
    beta_names <- names(out$beta_hat)
    data_model <- strsplit(strsplit(file_list[idx], "datamodel_")[[1]][2], "_fitmodel")[[1]][1]
    trial_idx <- as.numeric(strsplit(strsplit(file_list[idx], "trial_")[[1]][2], "_gene")[[1]][1])
    covariate_name <- beta_names[stringdist::amatch(strsplit(strsplit(file_list[idx], "covariate_")[[1]][2], ".rds")[[1]][1], beta_names, maxDist = Inf)]
    new_beta_val <- as.numeric(strsplit(strsplit(file_list[idx], "beta_")[[1]][2], "_covariate")[[1]][1])
    
    summary_df <- out$performanceSummary
    summary_df <- cbind(
      data.frame(
        datagen_model = data_model,
        trial = trial_idx,
        changed_covariate = covariate_name,
        new_beta_val = new_beta_val
      ),
      summary_df
    )
    summary_out_list[[idx]] <- summary_df
    
    wald_df <- out$wald_df
    wald_out_list[[idx]] <- wald_df
  }
  
  # Model performance results
  summary_out_all <- dplyr::bind_rows(summary_out_list)
  saveRDS(summary_out_all,
          file.path(processed_data_path, "power_perf_df.rds"),
          compress = "xz")
  
  # Wald test statistics for the within-cell type between-conditions comparison
  wald_df_all <- dplyr::bind_rows(wald_out_list)
  saveRDS(wald_df_all,
          file.path(processed_data_path, "power_wald_df.rds"),
          compress = "xz")
}
