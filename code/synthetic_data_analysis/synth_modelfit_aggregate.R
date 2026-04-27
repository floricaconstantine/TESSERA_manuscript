#' Aggregate TESSERA and GAM/GLM results from synthetic model fitting
#' @author Florica Constantine

# Libraries
library(TESSERA)
library(dplyr)

# Paths
path <-  "/scratch/users/spatialseq/natgen_kidney/"
processed_data_path <- file.path(path, "processed")
real_data_path <- file.path(path, "real_output_onlyinteract")
modelfit_data_path <- file.path(path, "synth_output_onlyinteract")
glm_gam_modelfit_data_path <- file.path(path, "glm_offset_synth_output_onlyinteract")

# Load prep_data object
prep_obj <- readRDS(file.path(processed_data_path, "prepData_lattice_onlyinteract.rds"))

# Genes
gene_keep <- c("IGKC", "DDX24", "SMTN", "CRLS1", "ERP44")

### TESSERA ###

# List of files to read in
pattern_type <- ".rds"
full_file_list <- base::list.files(modelfit_data_path, pattern = pattern_type)
file_list <- c()
for (f_idx in full_file_list) {
  gene_name <- sub(".*gene_(.*)\\.rds", "\\1", f_idx)
  if (gene_name %in% gene_keep) {
    file_list <- c(file_list, f_idx)
  }
}
length(file_list)

# Aggregate files
if (file.exists(file.path(
  processed_data_path,
  "synth_results_perfsummary_onlyinteract.rds"
))
&
file.exists(file.path(
  processed_data_path,
  "synth_results_waldstats_onlyinteract.rds"
))) {
  summary_out_all <- readRDS(file.path(
    processed_data_path,
    "synth_results_perfsummary_onlyinteract.rds"
  ))
  wald_df_all <- readRDS(file.path(
    processed_data_path,
    "synth_results_waldstats_onlyinteract.rds"
  ))
} else {
  real_data_list <- list()

  summary_out_list <- list()
  wald_out_list <- list()
  for (idx in 1:length(file_list)) {
    if (0 == (idx %% 250)) {
      cat(idx, "\n")
    }

    out <- readRDS(file.path(modelfit_data_path, file_list[idx]))
    # Add in information that wasn't stored
    trial_idx <- as.numeric(strsplit(strsplit(file_list[idx], "_trial_")[[1]][2], "_gene_")[[1]][1])

    if (grepl("_datamodel_Leroux", file_list[idx])) {
      data_model <- "Leroux"
    } else if (grepl("_datamodel_CAR", file_list[idx])) {
      data_model <- "CAR"
    } else if (grepl("_datamodel_SAR", file_list[idx])) {
      data_model <- "SAR"
    } else if (grepl("_datamodel_spNNGP", file_list[idx])) {
      data_model <- "spNNGP"
    } else {
      data_model <- strsplit(strsplit(strsplit(file_list[idx], "_datamodel_")[[1]][2], "_fitmodel_")[[1]][1], "_kernel_")[[1]][1]
      cat("Unknown datamodel",
          file_list[idx],
          "; Assigning: ",
          data_model,
          "\n")
    }

    # data_model <- strsplit(strsplit(strsplit(file_list[idx], "_datamodel_")[[1]][2], "_fitmodel_")[[1]][1], "_kernel_")[[1]][1]

    out$performanceSummary <- cbind(data.frame(data_model = data_model, trial = trial_idx),
                                    out$performanceSummary)


    # Add in parameter fit information
    # Extract real data fit information
    gene <- out$run_settings$gene_name
    kernel_type <- "Mat"

    # Need to include both gene and data model to look up parameters
    list_lookup_name <- paste0(data_model, "_", gene)
    if (list_lookup_name %in% names(real_data_list)) {

    } else {
      cat("Loading data for gene ", gene, idx, " model ", data_model, list_lookup_name, "\n")
      if (data_model == "spNNGP") {
        real_alg_out <- readRDS(
          file.path(
            real_data_path,
            "poisECM_real_natgen_kidney_",
            "model_",
            data_model,
            "_kernel_",
            kernel_type,
            "_gene_",
            gene,
            ".rds"
          )
        )
        # Extract covariance parameters
        cov_params_hat <- real_alg_out$cov_param_hat
        # Extract beta
        beta_hat <- real_alg_out$beta_hat

        gamma_hat <- NA
        tau2_hat <- NA
      } else{
        real_alg_out <- readRDS(
          file.path(
            real_data_path,
            "poisECM_real_natgen_kidney_",
            "model_",
            data_model,
            "_gene_",
            gene,
            ".rds"
          )
        )

        # Extract gamma and tau^2
        gamma_hat <- real_alg_out$gamma_hat
        tau2_hat <- real_alg_out$tau2_hat
        # Extract beta
        beta_hat <- real_alg_out$beta_hat

        cov_params_hat <- NA
      }

      real_data_list[[list_lookup_name]] <- list(
        beta_hat = beta_hat,
        gamma_hat = gamma_hat,
        tau2_hat = tau2_hat,
        cov_params_hat = cov_params_hat
      )
    }

    if (data_model == "spNNGP") {
      out$performanceSummary$nugget_true <- real_data_list[[list_lookup_name]]$cov_param_hat[, 1]
      out$performanceSummary$sill_true <- real_data_list[[list_lookup_name]]$cov_param_hat[, 2]
      out$performanceSummary$range_true <- real_data_list[[list_lookup_name]]$cov_param_hat[, 3]
      out$performanceSummary$smoothness_true <- real_data_list[[list_lookup_name]]$cov_param_hat[, 4]
    } else {
      out$performanceSummary$gamma_true <- real_data_list[[list_lookup_name]]$gamma_hat
      out$performanceSummary$tau2_true <- real_data_list[[list_lookup_name]]$tau2_hat
    }
    out$performanceSummary$SE_beta <- mean((real_data_list[[list_lookup_name]]$beta_hat - out$beta_hat)^2)
    out$performanceSummary$Norm2_beta <- mean((real_data_list[[list_lookup_name]]$beta_hat)^2)

    summary_out_list[[idx]] <- out$performanceSummary
    wald_out_list[[idx]] <- out$wald_df
  }
  summary_out_all <- dplyr::bind_rows(summary_out_list)
  saveRDS(
    summary_out_all,
    file.path(
      processed_data_path,
      "synth_results_perfsummary_onlyinteract.rds"
    ),
    compress = "xz"
  )
  wald_df_all <- dplyr::bind_rows(wald_out_list)
  saveRDS(
    wald_df_all,
    file.path(
      processed_data_path,
      "synth_results_waldstats_onlyinteract.rds"
    ),
    compress = "xz"
  )
}


### GLMs and GAMs ###

# List of files to read in
pattern_type <- ".rds"
full_file_list <- base::list.files(glm_gam_modelfit_data_path, pattern = pattern_type)
file_list <- c()
for (f_idx in full_file_list) {
  gene_name <- sub(".*gene_(.*)\\.rds", "\\1", f_idx)
  if (gene_name %in% gene_keep) {
    file_list <- c(file_list, f_idx)
  }
}
length(file_list)

# Aggregate files
if (file.exists(
  file.path(
    processed_data_path,
    "glm_gam_synth_results_perfsummary_onlyinteract.rds"
  )
)
&
file.exists(file.path(
  processed_data_path,
  "glm_gam_synth_results_waldstats_onlyinteract.rds"
))) {
  summary_out_all <- readRDS(
    file.path(
      processed_data_path,
      "glm_gam_synth_results_perfsummary_onlyinteract.rds"
    )
  )
  wald_df_all <- readRDS(
    file.path(
      processed_data_path,
      "glm_gam_synth_results_waldstats_onlyinteract.rds"
    )
  )
} else{
  summary_out_list <- list()
  wald_out_list <- list()
  for (idx in 1:length(file_list)) {
    if (0 == (idx %% 250)) {
      cat(idx, "\n")
    }

    out <- readRDS(file.path(glm_gam_modelfit_data_path, file_list[idx]))
    summary_out_list[[idx]] <- out$performanceSummary
    wald_out_list[[idx]] <- out$wald_df

  }
  summary_out_all <- dplyr::bind_rows(summary_out_list)
  saveRDS(
    summary_out_all,
    file.path(
      processed_data_path,
      "glm_gam_synth_results_perfsummary_onlyinteract.rds"
    ),
    compress = "xz"
  )
  wald_df_all <- dplyr::bind_rows(wald_out_list)
  saveRDS(
    wald_df_all,
    file.path(
      processed_data_path,
      "glm_gam_synth_results_waldstats_onlyinteract.rds"
    ),
    compress = "xz"
  )
}
