#' Run TESSERA on synthetic data based on real-world kidney data
#' @author Florica Constantine


## Libraries
library(TESSERA)
library(dplyr)
library(reshape2)
library(Rfast)
library(Matrix)
library(spatstat.geom)
library(ggplot2)

## Paths
path <-  "/scratch/users/spatialseq/natgen_kidney/"
processed_data_path <- file.path(path, "processed")
real_out_path <- file.path(path, "real_output_onlyinteract")
out_path <- file.path(path, "synth_output_onlyinteract")

## Simulation setup
datagen_model_list <- c("Leroux", "spNNGP")
fit_model_list <- c("Leroux", "spNNGP", "CAR", "SAR")
gene_list <- c("IGKC", "DDX24", "SMTN", "CRLS1", "ERP44") 
n_trials <- 1
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
} else {
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

## Wald statistics
# Get dataframe of Wald statistics
wald_df <- TESSERA::calc_Wald_statistics(alg_out, prep_synth_out$new_TESSERAData_obj$contrast_mat)
# Add in run parameters that aren't in wald_df (gene, fit_model present)
wald_df <- cbind(data.frame(datagen_model=data_model,
                            trial=trial_idx), wald_df)
# Add to alg_out
alg_out$wald_df <- wald_df


## Save to file
if (fit_model == "spNNGP") {
  outfile_base <- paste0(
    "poisECM_synth_natgen_kidney_",
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
    "poisECM_synth_natgen_kidney_",
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
saveRDS(alg_out,
        file = paste0(outfile_base, ".rds"),
        compress = TRUE)

## End
print("Done.")
