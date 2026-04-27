#' Run GAM, NB GLM, and Poisson GLM methods on non-spatial synthetic data based on real-world kidney data
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
datagen_model_list <- c("spNNGP")
fit_model_list <- c("GAM_GLM")
gene_list <- c("IGKC", "DDX24", "SMTN", "CRLS1", "ERP44")
n_trials <- 100
sim_df <- expand.grid(1:n_trials, fit_model_list, datagen_model_list, gene_list)
colnames(sim_df) <- c("trial", "fit_model", "datagen_model", "gene")

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


## Load in prepData object for data

prep_real_out <- readRDS(file.path(processed_data_path, "prepData_lattice_onlyinteract.rds"))

## Create new prepData object with counts sampled from Poisson lattice model

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


## Counts Analysis

Mean_counts2_sample <- sapply(prep_synth_out$new_TESSERAData_obj$counts_list, function (x) {
  mean(x^2)
})
Mean_counts2_total <- mean(Reduce(c, prep_synth_out$new_TESSERAData_obj$counts_list)^2)

# Store Moran of X beta_true and counts
n_s <- length(prep_synth_out$new_TESSERAData_obj$X_list)
Moran_Xbeta_true <- rep(0.0, n_s)
Moran_counts <- rep(0.0, n_s)
for (idx in 1:length(prep_synth_out$new_TESSERAData_obj$counts_list)) {
  Moran_Xbeta_true[idx] <- calc_moran(
    prep_synth_out$new_TESSERAData_obj$X_list[[idx]] %*% beta_hat,
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 1],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 2]
  )
  Moran_counts[idx] <- calc_moran(
    prep_synth_out$new_TESSERAData_obj$counts_list[[idx]][1, ],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 1],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 2]
  )
}

beta_norm2 <- sum(beta_hat^2)
sample_names <- names(prep_synth_out$new_TESSERAData_obj$X_list)

## Run algorithms and store results

beta_outs <- list()
perf_df <- list()
for (alg in c("GAM", "NB_GLM", "Poisson_GLM")) {
  if ("GAM" == alg) {
    alg_out <- mgcv_gam_wrapper(
      z_list = prep_synth_out$new_TESSERAData_obj$counts_list,
      X_list = prep_synth_out$new_TESSERAData_obj$X_list,
      coords_list = prep_synth_out$new_TESSERAData_obj$coords_list,
      library_size_list = prep_synth_out$new_TESSERAData_obj$library_size_list
    )
  } else if ("NB_GLM" == alg) {
    alg_out <- glm_nb_wrapper(
      z_list = prep_synth_out$new_TESSERAData_obj$counts_list,
      X_list = prep_synth_out$new_TESSERAData_obj$X_list,
      library_size_list = prep_synth_out$new_TESSERAData_obj$library_size_list
    )
  } else if ("Poisson_GLM" == alg) {
    alg_out <- glm_wrapper(
      z_list = prep_synth_out$new_TESSERAData_obj$counts_list,
      X_list = prep_synth_out$new_TESSERAData_obj$X_list,
      library_size_list = prep_synth_out$new_TESSERAData_obj$library_size_list
    )
  }

  fits <- predict(alg_out, type = "response")
  resids <- resid(alg_out, type = "response")

  # Store Moran of X beta_true
  n_s <- length(prep_synth_out$new_TESSERAData_obj$X_list)
  Moran_Xbeta <- rep(0.0, n_s)
  Moran_fit <- rep(0.0, n_s)
  Moran_resid <- rep(0.0, n_s)
  MSE_counts_sample <- rep(0.0, n_s)
  beta_out <- coef(alg_out)[1:ncol(prep_synth_out$new_TESSERAData_obj$X_list[[1]])]
  for (idx in 1:n_s) {
    Moran_Xbeta[idx] <- calc_moran(
      prep_synth_out$new_TESSERAData_obj$X_list[[idx]] %*% beta_out,
      prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 1],
      prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 2]
    )

    if (1 == idx) {
      start_idx <- 1
    } else {
      start_idx <- cumsum(sapply(prep_synth_out$new_TESSERAData_obj$X_list, nrow))[idx - 1] + 1
    }
    if (n_s == idx) {
      end_idx <- length(fits)
    } else {
      end_idx <- cumsum(sapply(prep_synth_out$new_TESSERAData_obj$X_list, nrow))[idx]
    }

    Moran_fit[idx] <- calc_moran(
      fits[start_idx:end_idx],
      prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 1],
      prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 2]
    )

    Moran_resid[idx] <- calc_moran(
      resids[start_idx:end_idx],
      prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 1],
      prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 2]
    )

    MSE_counts_sample[idx] <- mean(resids[start_idx:end_idx]^2)
  }
  MSE_counts_total <- mean(resids^2)

  # Store error in beta
  beta_SE <- sum((beta_out - beta_hat)^2)
  beta_outs[[alg]] <- beta_out

  # Create dataframe
  perf_df[[alg]] <- data.frame(
    sample=sample_names,
    fit_model = alg,
    trial = trial_idx,
    datagen_model = data_model,
    gene = gene_name,
    Moran_counts = Moran_counts,
    Moran_Xbeta_true = Moran_Xbeta_true,
    Moran_Xbeta = Moran_Xbeta,
    Moran_fit = Moran_fit,
    Moran_resid = Moran_resid,
    MSE_counts_sample = MSE_counts_sample,
    MSE_counts_total = MSE_counts_total,
    Mean_counts2_sample = Mean_counts2_sample,
    Mean_counts2_total = Mean_counts2_total,
    beta_SE = beta_SE,
    beta_norm2 = beta_norm2
  )
}
perf_df <- dplyr::bind_rows(perf_df)


## Save to file

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
outfile_base <- file.path(out_path, outfile_base)
saveRDS(perf_df,
        file = paste0(outfile_base, ".rds"),
        compress = TRUE)

## End
print("Done.")
