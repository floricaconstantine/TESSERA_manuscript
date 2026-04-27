#' Aggregate TESSERA and GAM/GLM results on non-spatial synthetic data based on real-world kidney data
#' @author Florica Constantine

# Libraries
library(TESSERA)
library(dplyr)
library(stringdist)

# Paths
path <-  "/scratch/users/spatialseq/natgen_kidney/"
modelfit_data_path <- file.path(path, "synth_nonspatial_output_onlyinteract")
processed_data_path <- file.path(path, "processed")

# Load in prep_data object
prep_obj <- readRDS(paste0(processed_data_path, "prepData_lattice_onlyinteract.rds"))


# Aggregate files for TESSERA results
pattern_type <- ".rds"
file_list <- base::list.files(modelfit_data_path, pattern = pattern_type)
file_list <- file_list[grepl("fitmodel_(Leroux|SAR|CAR|spNNGP)", file_list)]
if(file.exists(paste0(processed_data_path, "nonspatial_results_perfsummary_onlyinteract.rds"))
   & file.exists(paste0(processed_data_path, "nonspatial_results_waldstats_onlyinteract.rds"))){
  summary_out_all <- readRDS(paste0(processed_data_path, "nonspatial_results_perfsummary_onlyinteract.rds"))
  wald_df_all <- readRDS(paste0(processed_data_path, "nonspatial_results_waldstats_onlyinteract.rds"))
}else{
  summary_out_list <- list()
  wald_out_list <- list()
  for(idx in 1:length(file_list)){
    if (0 == (idx %% 250)) {
      cat(idx, "\n")
    }

    out <- readRDS(paste0(modelfit_data_path, file_list[idx]))

    beta_names <- names(out$beta_hat)
    data_model <- strsplit(strsplit(file_list[idx], "datamodel_")[[1]][2], "_")[[1]][1]
    trial_idx <- as.numeric(strsplit(strsplit(file_list[idx], "trial_")[[1]][2], "_gene")[[1]][1])
    covariate_name <- beta_names[stringdist::amatch(strsplit(strsplit(file_list[idx], "covariate_")[[1]][2], ".rds")[[1]][1],
                                                    beta_names, maxDist = Inf)]
    new_beta_val <- as.numeric(strsplit(strsplit(file_list[idx], "beta_")[[1]][2], "_covariate")[[1]][1])


    summary_df <- out$performanceSummary
    summary_df <- cbind(data.frame(datagen_model=data_model,
                                   trial=trial_idx,
                                   changed_covariate=covariate_name,
                                   new_beta_val=new_beta_val
    ), summary_df)
    summary_out_list[[idx]] <- summary_df

    wald_df <- out$wald_df
    wald_df$new_beta_val <- new_beta_val
    wald_out_list[[idx]] <- wald_df
  }
  summary_out_all <- dplyr::bind_rows(summary_out_list)
  saveRDS(summary_out_all, paste0(processed_data_path, "nonspatial_results_perfsummary_onlyinteract.rds"))
  wald_df_all <- dplyr::bind_rows(wald_out_list)
  saveRDS(wald_df_all, paste0(processed_data_path, "nonspatial_results_waldstats_onlyinteract.rds"))
}

# Aggregate files for GAM/GLM results
pattern_type <- ".rds"
file_list <- base::list.files(modelfit_data_path, pattern = pattern_type)
file_list <- file_list[grepl("fitmodel_GAM_GLM", file_list)]
if (file.exists(paste0(processed_data_path, "nonspatial_gam_glm_perfsummary_onlyinteract.rds"))) {
  summary_out_all <- readRDS(paste0(processed_data_path, "nonspatial_gam_glm_perfsummary_onlyinteract.rds"))
} else {
  summary_out_list <- list()
  for(idx in 1:length(file_list)){
    if (0 == (idx %% 250)) {
      cat(idx, "\n")
    }

    out <- readRDS(paste0(modelfit_data_path, file_list[idx]))
    summary_out_list[[idx]] <- out
  }
  summary_out_all <- dplyr::bind_rows(summary_out_list)
  saveRDS(summary_out_all, paste0(processed_data_path, "nonspatial_gam_glm_perfsummary_onlyinteract.rds"))
}


