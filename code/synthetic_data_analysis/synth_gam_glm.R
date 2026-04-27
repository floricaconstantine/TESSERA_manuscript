#' Run GAM, NB-GLM, and Poisson-GLM methods on synthetic data based on real-world kidney data
#' @author Florica Constantine


## Scripts
source("TESSERA_comparison_methods.R")


## Libraries

library(dplyr)
library(Matrix)
library(TESSERA)


## Paths

path <-  "/scratch/users/spatialseq/natgen_kidney/"
processed_data_path <- file.path(path, "processed")
real_out_path <- file.path(path, "real_output_onlyinteract")
out_path <- file.path(path, "glm_offset_synth_output_onlyinteract")

## Simulation setup

datagen_model_list <- c("Leroux", "spNNGP")
fit_model_list <- c("Leroux", "spNNGP", "Poisson_GLM", "NB_GLM", "GAM")
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

data_model <- as.character(sim_df[run_idx, ]$datagen_model)
fit_model <- as.character(sim_df[run_idx, ]$fit_model)
gene_name <- as.character(sim_df[run_idx, ]$gene)
trial_idx <- sim_df[run_idx, ]$trial
print(sim_df[run_idx, ])

if ("spNNGP" == as.character(data_model)) {
  data_kernel_type <- "Mat"
} else {
  data_kernel_type <- "NA"
}
if ("spNNGP" == as.character(fit_model)) {
  fit_kernel_type <- "Mat"
} else {
  fit_kernel_type <- "NA"
}


## Output file

outfile_base <- "poisECM_synth_natgen_kidney"
if ("spNNGP" == data_model) {
  dm_string <- paste0("data_model_", data_model, "_datakernel_", data_kernel_type)
} else {
  dm_string <- paste0("data_model_", data_model, "_datakernel_", "NA")
}
if ("spNNGP" == fit_model) {
  fm_string <- paste0("fit_model_", fit_model, "_fitkernel_", fit_kernel_type)
} else {
  fm_string <- paste0("fit_model_", fit_model, "_fitkernel_", "NA")
}
outfile <- paste0(
  out_path,
  outfile_base,
  "_",
  dm_string,
  "_",
  fm_string,
  "_trial_",
  trial_idx,
  "_gene_",
  gene_name,
  ".rds"
)
print(outfile)


## Global parameters

GLOBAL_SEED <- 2025
DATA_SEED <- GLOBAL_SEED + trial_idx
set.seed(DATA_SEED)


## TESSERA parameters

nngp_k <- 20
em_iters <- 200
opt_iters <- 5
em_min_iters <- 30
em_tol <- 1e-3
em_stopping <- "rel_loglike"
beta_init <- "glm"
cov_fit_method <- "BRISC"
cov_init <- "BRISC"
TESSERA_verbose <- TRUE
gamma_init <- "moran"
tau2_init <- "var"
max_gamma <- 0.999


## Load in real data results

if (data_model == "spNNGP") {
  real_alg_out <- readRDS(
    paste0(
      real_out_path,
      "poisECM_real_natgen_kidney_",
      "model_",
      data_model,
      "_kernel_",
      data_kernel_type,
      "_gene_",
      gene_name,
      ".rds"
    )
  )
  # Extract covariance parameters
  cov_params_hat <- real_alg_out$cov_param_hat
  # Extract beta
  beta_hat <- real_alg_out$beta_hat

  # Unused parameters
  gamma_hat <- NA
  tau2_hat <- NA
} else{
  real_alg_out <- readRDS(
    paste0(
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

  # Unused parameters
  cov_params_hat <- NA
}


## Load in prep_data object for data

if (data_model == "spNNGP") {
  prep_real_out <- readRDS(paste0(processed_data_path, "prepData_spNNGP_onlyinteract.rds"))
} else{
  prep_real_out <- readRDS(paste0(processed_data_path, "prepData_lattice_onlyinteract.rds"))
}


## Create new prep_data object with counts sampled from Poisson lattice model

if (data_model == "spNNGP") {
  prep_synth_out <- TESSERA::prep_synth_data(
    TESSERAData_obj = prep_real_out,
    gene_list = gene_name,
    data_gen_model = data_model,
    cov_params = cov_params_hat,
    cov_type = data_kernel_type,
    nngp_k = nngp_k,
    beta_true = beta_hat
  )
} else{
  prep_synth_out <- TESSERA::prep_synth_data(
    TESSERAData_obj = prep_real_out,
    gene_list = gene_name,
    data_gen_model = data_model,
    tau2_true = tau2_hat,
    gamma_true = pmin(gamma_hat, max_gamma),
    beta_true = beta_hat
  )
}


## Run algorithm

if (fit_model == "spNNGP") {
  alg_out <- TESSERA::TESSERA_spNNGP(
    TESSERAData_obj = prep_synth_out$new_TESSERAData_obj,
    gene_name = gene_name,
    cov_type = fit_kernel_type,
    nngp_k = nngp_k,
    em_iters = em_iters,
    opt_iters = opt_iters,
    em_min_iters = em_min_iters,
    em_tol = em_tol,
    em_stopping = em_stopping,
    beta_init = beta_init,
    cov_fit_method = cov_fit_method,
    cov_init = cov_init,
    verbose = TESSERA_verbose
  )
} else if (fit_model %in% c("Leroux", "CAR", "SAR")) {
  alg_out <- TESSERA::TESSERA_lattice(
    TESSERAData_obj = prep_synth_out$new_TESSERAData_obj,
    gene_name = gene_name,
    model_type = fit_model,
    em_iters = em_iters,
    opt_iters = opt_iters,
    em_min_iters = em_min_iters,
    em_tol = em_tol,
    em_stopping = em_stopping,
    beta_init = beta_init,
    gamma_init = gamma_init,
    tau2_init = tau2_init,
    verbose = TESSERA_verbose
  )
} else if (fit_model == "GAM") {
  alg_out <- mgcv_gam_wrapper(
    z_list = prep_synth_out$new_TESSERAData_obj$counts_list,
    X_list = prep_synth_out$new_TESSERAData_obj$X_list,
    coords_list = prep_synth_out$new_TESSERAData_obj$coords_list,
    library_size_list = prep_synth_out$new_TESSERAData_obj$library_size_list
  )
  alg_out$df_total <- sum(alg_out$edf)
} else if (fit_model == "Poisson_GLM") {
  alg_out <- glm_wrapper(
    prep_synth_out$new_TESSERAData_obj$counts_list,
    prep_synth_out$new_TESSERAData_obj$X_list,
    prep_synth_out$new_TESSERAData_obj$library_size_list
  )
} else if (fit_model == "NB_GLM") {
  alg_out <- glm_nb_wrapper(
    prep_synth_out$new_TESSERAData_obj$counts_list,
    prep_synth_out$new_TESSERAData_obj$X_list,
    prep_synth_out$new_TESSERAData_obj$library_size_list
  )
}


## Wald statistics

contrast_mat <- prep_synth_out$new_TESSERAData_obj$contrast_mat
if (fit_model %in% c("Leroux", "CAR", "SAR", "spNNGP")) {
  # Get dataframe of Wald statistics
  wald_df <- calc_Wald_statistics(alg_out, contrast_mat)
  # Add in run parameters that aren't in wald_df (gene, fit_model present)
  wald_df <- cbind(
    data.frame(
      datagen_model = data_model,
      data_kernel_type = data_kernel_type,
      trial = trial_idx
    ),
    wald_df
  )
  # Drop pointless columns
  wald_df <- wald_df[, which(colnames(wald_df) != "contrast_string")]
  # Add to alg_out
  alg_out$wald_df <- wald_df
} else if (fit_model %in% c("Poisson_GLM", "NB_GLM", "GAM")) {
  wald_df <- list()
  for (c_idx in 1:nrow(contrast_mat)) {
    # Contrast
    R_beta <- as.numeric(t(contrast_mat[c_idx, ]) %*% alg_out$coefficients[1:ncol(contrast_mat)])
    # SE^2
    RVR <- as.numeric(t(contrast_mat[c_idx, ]) %*% vcov(alg_out)[1:ncol(contrast_mat), 1:ncol(contrast_mat)] %*% (contrast_mat[c_idx, ]))

    wald_df[[c_idx]] <- data.frame(
      datagen_model = data_model,
      trial = trial_idx,
      gene = gene_name,
      fit_model = fit_model,
      data_kernel_type = data_kernel_type,
      fit_kernel_type = fit_kernel_type,
      contrast_indices = paste0(which(contrast_mat[c_idx, ] != 0), collapse =
                                  "_"),
      contrast_val = R_beta,
      contrast_se = sqrt(RVR),
      wald_stat_t = R_beta / sqrt(RVR)
    )
  }

  # Add to alg_out
  alg_out$wald_df <- dplyr::bind_rows(wald_df)

  alg_out$beta_hat <- alg_out$coefficients[1:length(beta_hat)]
  alg_out$beta_vcov <- vcov(alg_out)[1:length(beta_hat), 1:length(beta_hat)]
}


## Get statistics for GLM/GAM

if (fit_model %in% c("Poisson_GLM", "NB_GLM", "GAM")) {
  MSE <- rep(NA,
             length(prep_synth_out$new_TESSERAData_obj$counts_list))
  counts2 <- rep(NA, length(MSE))
  Moran_pred <- rep(NA, length(MSE))
  Moran_counts <- rep(NA, length(MSE))
  Moran_resid <- rep(NA, length(MSE))
  sample_indexing <- c(0, cumsum(
    sapply(prep_synth_out$new_TESSERAData_obj$counts_list, length)
  ))

  # Extract predictions
  alg_pred <- predict(alg_out, type = "response")
  lib_size_list <- Reduce(c, prep_synth_out$new_TESSERAData_obj$library_size_list)
  # alg_pred <- alg_pred * lib_size_list # Uncomment if not using offset
  alg_resid <- Reduce(c, prep_synth_out$new_TESSERAData_obj$counts_list) - alg_pred
  for (s_idx in 1:length(MSE)) {
    fit_local <- alg_pred[(1 + sample_indexing[s_idx]):sample_indexing[s_idx + 1]]
    resid_local <- alg_resid[(1 + sample_indexing[s_idx]):sample_indexing[s_idx + 1]]
    MSE[s_idx] <- mean(resid_local^2, na.rm = TRUE)
    counts2[s_idx] <- mean(prep_synth_out$new_TESSERAData_obj$counts_list[[s_idx]]^2)

    Moran_resid[s_idx] <- calc_moran(
      resid_local,
      prep_synth_out$new_TESSERAData_obj$coords_list[[s_idx]][, 1],
      prep_synth_out$new_TESSERAData_obj$coords_list[[s_idx]][, 2]
    )[1]
    Moran_pred[s_idx] <- calc_moran(
      fit_local,
      prep_synth_out$new_TESSERAData_obj$coords_list[[s_idx]][, 1],
      prep_synth_out$new_TESSERAData_obj$coords_list[[s_idx]][, 2]
    )[1]
    Moran_counts[s_idx] <- calc_moran(
      prep_synth_out$new_TESSERAData_obj$counts_list[[s_idx]],
      prep_synth_out$new_TESSERAData_obj$coords_list[[s_idx]][, 1],
      prep_synth_out$new_TESSERAData_obj$coords_list[[s_idx]][, 2]
    )[1]
  }
  MSE_total <- mean(alg_out$residuals^2)
  counts2_total <- mean(Reduce(c, prep_synth_out$new_TESSERAData_obj$counts_list)^2)
  beta_SE <- mean((beta_hat - coef(alg_out)[1:length(beta_hat)])^2)
  beta_hat_norm2 <- mean(beta_hat^2)

  # Add to object
  summary_out_all <- data.frame(
    gene = sapply(prep_synth_out$synthetic_count_summary$gene, function(x) {
      x[1]
    }),
    fit_model = fit_model,
    data_model = data_model,
    data_kernel_type = data_kernel_type,
    trial_idx = trial_idx,
    sample = sapply(prep_synth_out$synthetic_count_summary$sample, function(x) {
      x[1]
    }),
    n_cells = sapply(prep_synth_out$new_TESSERAData_obj$coords_list, nrow),
    MSE_counts_sample = MSE,
    MSE_counts_total = MSE_total,
    Mean_counts2_sample = counts2,
    Mean_counts2_total = counts2_total,
    Moran_counts = Moran_counts,
    Moran_predictions = Moran_pred,
    Moran_residuals = Moran_resid,
    SE_beta = beta_SE,
    Norm2_beta = beta_hat_norm2,
    Group = sapply(prep_synth_out$new_TESSERAData_obj$covariates_list, function(x) {
      x$Group[1]
    })
  )
  alg_out$performanceSummary <- summary_out_all

  # Store likelihood
  alg_out$data_ll <- as.numeric(logLik(alg_out))
} else {
  beta_SE <- mean((beta_hat - alg_out$beta_hat)^2)
  beta_hat_norm2 <- mean(beta_hat^2)

  alg_out$performanceSummary$SE_beta <- beta_SE
  alg_out$performanceSummary$Norm2_beta <- beta_hat_norm2
  alg_out$performanceSummary$Group <- sapply(prep_synth_out$new_TESSERAData_obj$covariates_list, function(x) {
    x$Group[1]
  })
  alg_out$performanceSummary$data_model <- data_model
  alg_out$performanceSummary$data_kernel_type <- data_kernel_type
  alg_out$performanceSummary$trial_idx <- trial_idx
}

# Add in gamma/tau^2
if (data_model %in% c("Leroux", "CAR", "SAR")) {
  alg_out$performanceSummary$gamma_true <- pmin(gamma_hat, max_gamma)
  alg_out$performanceSummary$tau2_true <- tau2_hat
} else {
  alg_out$performanceSummary$nugget_true <- cov_params_hat[, 1]
  alg_out$performanceSummary$sill_true <- cov_params_hat[, 2]
  alg_out$performanceSummary$range_true <- cov_params_hat[, 3]
  alg_out$performanceSummary$smoothness_true <- cov_params_hat[, 4]
}


## Save to file

TESSERA_keep_list <- c(
  "beta_hat",
  "gamma_hat",
  "tau2_hat",
  "data_log_like_tracker",
  "expected_log_like_tracker",
  "tau2_neghessian",
  "gamma_neghessian",
  "beta_neghessian",
  "run_settings",
  "performanceSummary",
  "wald_df",
  "time",
  "start_idx_list"
)
TESSERA_spNNGP_keep_list <- c(
  "beta_hat",
  "cov_param_hat",
  "data_log_like_tracker",
  "expected_log_like_tracker",
  "beta_neghessian",
  "run_settings",
  "performanceSummary",
  "wald_df",
  "time",
  "start_idx_list"
)
Poisson_GLM_keep_list <- c(
  "coefficients",
  "time",
  "wald_df",
  "beta_hat",
  "beta_vcov",
  "performanceSummary",
  "data_ll"
)
NB_GLM_keep_list <- c(
  "coefficients",
  "time",
  "wald_df",
  "beta_hat",
  "beta_vcov",
  "performanceSummary",
  "data_ll",
  "theta"
)
GAM_keep_list <- c(
  "coefficients",
  "time",
  "wald_df",
  "beta_hat",
  "beta_vcov",
  "performanceSummary",
  "data_ll",
  "df.null",
  "df.residual",
  "min.edf",
  "df_total"
)

output_object <- list()
if (fit_model == "spNNGP") {
  for (nm in TESSERA_spNNGP_keep_list) {
    output_object[[nm]] <- alg_out[[nm]]
  }
} else if (fit_model %in% c("Leroux", "CAR", "SAR")) {
  for (nm in TESSERA_keep_list) {
    output_object[[nm]] <- alg_out[[nm]]
  }
} else if (fit_model == "GAM") {
  for (nm in GAM_keep_list) {
    output_object[[nm]] <- alg_out[[nm]]
  }
} else if (fit_model == "Poisson_GLM") {
  for (nm in Poisson_GLM_keep_list) {
    output_object[[nm]] <- alg_out[[nm]]
  }
} else if (fit_model == "NB_GLM") {
  for (nm in NB_GLM_keep_list) {
    output_object[[nm]] <- alg_out[[nm]]
  }
}


saveRDS(output_object,
        file = outfile,
        compress = TRUE)


## End

print("Done.")
