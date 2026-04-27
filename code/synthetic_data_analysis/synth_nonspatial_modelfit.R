#' Run TESSERA on non-spatial synthetic data based on real-world kidney data
#' @author Florica Constantine


## Scripts
source("TESSERA_comparison_methods.R")

## Libraries
library(TESSERA)
library(dplyr)
library(reshape2)
library(Rfast)
library(Matrix)
library(spatstat.geom)
library(ggplot2)

## Path
path <-  "/scratch/users/spatialseq/natgen_kidney/"
processed_data_path <- file.path(path, "processed")
real_out_path <- file.path(path, "real_output_onlyinteract")
out_path <- file.path(path, "synth_nonspatial_output_onlyinteract")

## Simulation setup
datagen_model_list <- c("Leroux", "spNNGP")
fit_model_list <- c("Leroux", "spNNGP", "CAR", "SAR")
gene_list <- c("IGKC", "DDX24", "SMTN", "CRLS1", "ERP44")
n_trials <- 100
sim_df <- expand.grid(1:n_trials, fit_model_list, datagen_model_list, gene_list)
colnames(sim_df) <- c("trial", "fit_model", "datagen_model", "gene")

# Clusterize
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
  if ((1 > run_idx) | (nrow(sim_df) < run_idx)) {
    stop("run_idx must be between 1 and # of possible jobs")
  }
}

kernel_type <- "Mat"
data_model <- sim_df[run_idx, ]$datagen_model
fit_model <- sim_df[run_idx, ]$fit_model
gene_name <- sim_df[run_idx, ]$gene
trial_idx <- sim_df[run_idx, ]$trial
print(sim_df[run_idx, ])

## Global parameters
GLOBAL_SEED <- 2025
DATA_SEED <- GLOBAL_SEED + trial_idx
set.seed(DATA_SEED)

## Load in real data results
if (data_model == "spNNGP") {
  real_alg_out <- readRDS(
    file.path(
      real_out_path,
      "poisECM_real_natgen_kidney_",
      "model_",
      data_model,
      "_kernel_",
      kernel_type,
      "_gene_",
      gene_name,
      ".rds"
    )
  )
  # Extract covariance parameters
  cov_params_hat <- real_alg_out$cov_param_hat
  # Extract beta
  beta_hat <- real_alg_out$beta_hat
} else{
  real_alg_out <- readRDS(
    file.path(
      real_out_path,
      "poisECM_real_natgen_kidney_",
      "model_",
      data_model,
      "_gene_",
      gene_name,
      ".rds"
    )
  )

  # Extract gamma and tau^2
  gamma_hat <- real_alg_out$gamma_hat
  tau2_hat <- real_alg_out$tau2_hat
  # Extract beta
  beta_hat <- real_alg_out$beta_hat
}

# Make data non-spatial
if (data_model == "spNNGP") {
  cov_params_hat[, 2] <- 0
} else {
  gamma_hat <- 0 * gamma_hat
}


## Load in prep_data object for data
prep_real_out <- readRDS(file.path(processed_data_path, "prepData_lattice_onlyinteract.rds"))

## Create new prep_data object with counts sampled from Poisson lattice model
if (data_model == "spNNGP") {
  prep_synth_out <- TESSERA::prep_synth_data(
    TESSERAData_obj = prep_real_out,
    gene_list = gene_name,
    data_gen_model = data_model,
    cov_params = cov_params_hat,
    cov_type = kernel_type,
    nngp_k = 20,
    beta_true = beta_hat
  )
} else{
  prep_synth_out <- TESSERA::prep_synth_data(
    TESSERAData_obj = prep_real_out,
    gene_list = gene_name,
    data_gen_model = data_model,
    tau2_true = tau2_hat,
    gamma_true = pmin(gamma_hat, 0.999),
    beta_true = beta_hat
  )
}

## Run TESSERA algorithm
if (fit_model == "spNNGP") {
  alg_out <- TESSERA::TESSERA_spNNGP(
    TESSERAData_obj = prep_synth_out$new_TESSERAData_obj,
    gene_name = gene_name,
    cov_type = kernel_type,
    nngp_k = 20,
    em_iters = 200,
    opt_iters = 5,
    em_min_iters = 30,
    em_tol = 1e-3,
    em_stopping = "rel_loglike",
    beta_init = "glm",
    cov_fit_method = "BRISC",
    cov_init = "BRISC",
    verbose = TRUE
  )
} else if (fit_model %in% c("CAR", "SAR", "Leroux")) {
  alg_out <- TESSERA::TESSERA_lattice(
    TESSERAData_obj = prep_synth_out$new_TESSERAData_obj,
    gene_name = gene_name,
    model_type = fit_model,
    em_iters = 200,
    opt_iters = 5,
    em_min_iters = 30,
    em_tol = 1e-3,
    em_stopping = "rel_loglike",
    beta_init = "glm",
    gamma_init = "moran",
    tau2_init = "var",
    verbose = TRUE
  )
}

# Store error in beta
alg_out$performanceSummary$beta_SE <- sum((alg_out$beta_hat - beta_hat)^2)
alg_out$performanceSummary$beta_norm2 <- sum(beta_hat^2)

# Store Moran of X beta_true
alg_out$performanceSummary$Moran_Xbeta_true <- 0.0
for (idx in 1:length(prep_synth_out$new_TESSERAData_obj$counts_list)) {
  alg_out$performanceSummary$Moran_Xbeta_true[idx] <- calc_moran(
    prep_synth_out$new_TESSERAData_obj$X_list[[idx]] %*% beta_hat,
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 1],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 2]
  )
}


## Fit GLM

glm_out <- glm_wrapper(
  prep_synth_out$new_TESSERAData_obj$counts_list,
  prep_synth_out$new_TESSERAData_obj$X_list
)
glm_fits <- predict(glm_out, type = "response")
glm_resid <- resid(glm_out, type = "response")
# Store Moran of X beta_true
alg_out$performanceSummary$Moran_Xbeta_glm <- 0.0
alg_out$performanceSummary$Moran_glm_fit <- 0.0
alg_out$performanceSummary$Moran_glm_resid <- 0.0
alg_out$performanceSummary$MSE_counts_sample_glm <- 0.0
glm_beta <- coef(glm_out)
for (idx in 1:length(prep_synth_out$new_TESSERAData_obj$counts_list)) {
  alg_out$performanceSummary$Moran_Xbeta_glm[idx] <- calc_moran(
    prep_synth_out$new_TESSERAData_obj$X_list[[idx]] %*% glm_beta,
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 1],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 2]
  )

  if (length(prep_synth_out$new_TESSERAData_obj$counts_list) == idx) {
    end_idx <- length(glm_fits)
  } else {
    end_idx <- alg_out$start_idx_list[idx + 1]
  }

  alg_out$performanceSummary$Moran_glm_fit[idx] <- calc_moran(
    glm_fits[alg_out$start_idx_list[idx]:end_idx],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 1],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 2]
  )

  alg_out$performanceSummary$Moran_glm_resid[idx] <- calc_moran(
    glm_resid[alg_out$start_idx_list[idx]:end_idx],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 1],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 2]
  )

  alg_out$performanceSummary$MSE_counts_sample_glm[idx] <- mean(glm_resid[alg_out$start_idx_list[idx]:end_idx]^2)

}
alg_out$performanceSummary$MSE_counts_total_glm <- mean(glm_resid^2)

# Store error in beta
alg_out$performanceSummary$beta_SE_glm <- sum((coef(glm_out) - beta_hat)^2)
alg_out$beta_hat_glm <- coef(glm_out)

## Fit NB GLM

glm_out <- glm_nb_wrapper(
  prep_synth_out$new_TESSERAData_obj$counts_list,
  prep_synth_out$new_TESSERAData_obj$X_list
)
glm_fits <- predict(glm_out, type = "response")
glm_resid <- resid(glm_out, type = "response")
# Store Moran of X beta_true
alg_out$performanceSummary$Moran_Xbeta_nbglm <- 0.0
alg_out$performanceSummary$Moran_nbglm_fit <- 0.0
alg_out$performanceSummary$Moran_nbglm_resid <- 0.0
alg_out$performanceSummary$MSE_counts_sample_nbglm <- 0.0
glm_beta <- coef(glm_out)
for (idx in 1:length(prep_synth_out$new_TESSERAData_obj$counts_list)) {
  alg_out$performanceSummary$Moran_Xbeta_nbglm[idx] <- calc_moran(
    prep_synth_out$new_TESSERAData_obj$X_list[[idx]] %*% glm_beta,
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 1],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 2]
  )

  if (length(prep_synth_out$new_TESSERAData_obj$counts_list) == idx) {
    end_idx <- length(glm_fits)
  } else {
    end_idx <- alg_out$start_idx_list[idx + 1]
  }

  alg_out$performanceSummary$Moran_nbglm_fit[idx] <- calc_moran(
    glm_fits[alg_out$start_idx_list[idx]:end_idx],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 1],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 2]
  )

  alg_out$performanceSummary$Moran_nbglm_resid[idx] <- calc_moran(
    glm_resid[alg_out$start_idx_list[idx]:end_idx],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 1],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 2]
  )

  alg_out$performanceSummary$MSE_counts_sample_nbglm[idx] <- mean(glm_resid[alg_out$start_idx_list[idx]:end_idx]^2)

}
alg_out$performanceSummary$MSE_counts_total_nbglm <- mean(glm_resid^2)

# Store error in beta
alg_out$performanceSummary$beta_SE_nbglm <- sum((coef(glm_out) - beta_hat)^2)
alg_out$beta_hat_nbglm <- coef(glm_out)

## Fields to keep in alg_out
keep_list <- c(
  "beta_hat",
  "gamma_hat",
  "tau2_hat",
  "cov_param_hat",
  "performanceSummary",
  "beta_hat_nbglm",
  "beta_hat_glm"
)

## Save to file

if (fit_model == "spNNGP") {
  outfile_base <- paste0(
    "poisECM_synth_natgen_kidney_non_spatial_",
    "datamodel_",
    data_model,
    "_kernel_",
    kernel_type,
    "_fitmodel_",
    fit_model,
    "_kernel_",
    kernel_type,
    "_trial_",
    trial_idx,
    "_gene_",
    gene_name
  )
} else{
  outfile_base <- paste0(
    "poisECM_synth_natgen_kidney_nonspatial_",
    "datamodel_",
    data_model,
    "_fitmodel_",
    fit_model,
    "_trial_",
    trial_idx,
    "_gene_",
    gene_name
  )
}
outfile_base <- file.path(out_path, outfile_base)
saveRDS(alg_out[keep_list],
        file = paste0(outfile_base, ".rds"),
        compress = TRUE)

## End
print("Done.")
