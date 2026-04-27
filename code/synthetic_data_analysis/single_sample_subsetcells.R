#' Create and run TESSERA and competing methods on a single sample, using all but 4 cells, from synthetic data based on real-world kidney data
#' @author Florica Constantine


## Scripts
source("TESSERA_comparison_methods.R")

## Libraries
library(Matrix)
library(TESSERA)

# Path
path <-  "/scratch/users/spatialseq/natgen_kidney/"
processed_data_path <- file.path(path, "processed")
real_out_path <- file.path(path, "real_output_onlyinteract")
out_path <- file.path(path, "singlesample_output_onlyinteract")

# Sample to keep
sample_name <- "HK3035_ST"

## Simulation setup
datagen_model_list <- c("Leroux", "spNNGP")
fit_model_list <- c(
  "Leroux",
  "spNNGP",
  "CAR",
  "SAR",
  "GLM_Poisson",
  "GLM_NB",
  "MGCV",
  "LM",
  "BRISC",
  "MCMC_Leroux",
  "MCMC_CAR"
)
gene_list <- c(
  "IGKC",
  "DDX24",
  "SMTN",
  "CRLS1",
  "ERP44"
)

n_trials <- 100
sim_df <- expand.grid(1:n_trials, fit_model_list, gene_list, datagen_model_list)
colnames(sim_df) <- c("trial", "fit_model", "gene", "datagen_model")

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

  # Subset to sample
  cov_params_hat <- cov_params_hat[which(rownames(cov_params_hat) == sample_name), ]

  gamma_hat <- NA
  tau2_hat <- NA
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

  # Subset to sample
  gamma_hat <- gamma_hat[which(names(gamma_hat) == sample_name)]
  tau2_hat <- tau2_hat[which(names(tau2_hat) == sample_name)]

  cov_params_hat <- c(NA, NA, NA, NA)
}

## Load in prep_data object for data
prep_out <- readRDS(file.path(processed_data_path, "prepData_lattice_onlyinteract.rds"))

# Subset to sample
for (ky in names(prep_out)) {
  if (is.list(prep_out[[ky]])) {
    prep_out[[ky]] <- list(prep_out[[ky]][[sample_name]])
    names(prep_out[[ky]])[1] <- sample_name
  }
}

# Cell types to drop
drop_cell_types <- c("RBC", "Endo_Lym")
drop_cols_ct <- which(apply(sapply(drop_cell_types, function (x) {
  grepl(x, names(beta_hat))
}), 1, any))

# Edit design matrix: Drop columns that are not estimable
drop_cols_0 <- which(0 == colSums(prep_out$X_list[[sample_name]]))
drop_cols_int <- which(nrow(prep_out$X_list[[sample_name]])
                       == colSums(prep_out$X_list[[sample_name]]))
drop_cols <- sort(unique(c(drop_cols_ct, drop_cols_0, drop_cols_int)))
prep_out$X_list[[sample_name]] <- prep_out$X_list[[sample_name]][, -drop_cols]
prep_out$contrast_mat <- prep_out$contrast_mat[, -drop_cols]
# Check that design is full rank
stopifnot(Matrix::rankMatrix(prep_out$X_list[[sample_name]])[1]
          == ncol(prep_out$X_list[[sample_name]]))
# Drop entries of beta that correspond to dropped columns
beta_hat <- beta_hat[-drop_cols]
stopifnot(length(beta_hat) == ncol(prep_out$X_list[[sample_name]]))

# Drop cells
drop_rows <- which(prep_out$covariates_list[[sample_name]]$celltype %in% drop_cell_types)
prep_out$coords_list[[sample_name]] <- prep_out$coords_list[[sample_name]][-drop_rows, ]
prep_out$covariates_list[[sample_name]] <- prep_out$covariates_list[[sample_name]][-drop_rows, ]
prep_out$library_size_list[[sample_name]] <- prep_out$library_size_list[[sample_name]][-drop_rows]
# Subset design matrix
prep_out$X_list[[sample_name]] <- prep_out$X_list[[sample_name]][-drop_rows, ]
stopifnot(0 < min(rowSums(prep_out$X_list[[sample_name]])))
# Check that design is full rank
stopifnot(Matrix::rankMatrix(prep_out$X_list[[sample_name]])[1]
          == ncol(prep_out$X_list[[sample_name]]))
# Subset counts
prep_out$counts_list[[sample_name]] <- prep_out$counts_list[[sample_name]][, -drop_rows, drop =
                                                                             FALSE]
# Subset adjacency matrix and check
prep_out$W_list[[sample_name]] <- prep_out$W_list[[sample_name]][-drop_rows, -drop_rows]
stopifnot(0 == sum(abs(prep_out$W_list[[sample_name]] - t(prep_out$W_list[[sample_name]]))))
stopifnot(1 == max(prep_out$W_list[[sample_name]]))
stopifnot(0 == min(prep_out$W_list[[sample_name]]))
# Recompute degree matrix
prep_out$D_list[[sample_name]] <- Matrix::Diagonal(nrow(prep_out$W_list[[sample_name]]), Matrix::rowSums(prep_out$W_list[[sample_name]]), )
stopifnot(0 == sum(abs(
  diag(prep_out$D_list[[sample_name]]) - rowSums(prep_out$W_list[[sample_name]])
)))

# Make sure no cells are isolated now
stopifnot(0 < min(diag(prep_out$D_list[[sample_name]])))

# Recompute eigenvalues
prep_out$eig_CS_list[[sample_name]] <- Re(eigen(
  Matrix::solve(prep_out$D_list[[sample_name]], prep_out$W_list[[sample_name]]),
  FALSE,
  only.values = TRUE
)$values)
prep_out$eig_L_list[[sample_name]] <- Re(eigen(prep_out$D_list[[sample_name]] - prep_out$W_list[[sample_name]], TRUE, only.values = TRUE)$values)

## Create new prep_data object with counts sampled from Poisson lattice model
if (data_model == "spNNGP") {
  prep_synth_out <- TESSERA::prep_synth_data(
    TESSERAData_obj = prep_out,
    gene_list = gene_name,
    data_gen_model = data_model,
    cov_params = cov_params_hat,
    cov_type = kernel_type,
    nngp_k = 20,
    beta_true = beta_hat
  )
} else{
  prep_synth_out <- TESSERA::prep_synth_data(
    TESSERAData_obj = prep_out,
    gene_list = gene_name,
    data_gen_model = data_model,
    tau2_true = tau2_hat,
    gamma_true = pmin(gamma_hat, 0.999),
    beta_true = beta_hat
  )
}

# Methods
if (fit_model == "spNNGP") {
  ## Run TESSERA algorithm using spNNGP model
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
} else if (fit_model == "CAR" ||
           fit_model == "SAR" || fit_model == "Leroux") {
  ## Run TESSERA algorithm using a lattice model
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
} else if (fit_model == "BRISC") {
  # BRISC
  alg_out <- BRISC_wrapper(
    z_list = prep_synth_out$new_TESSERAData_obj$counts_list,
    X_list = prep_synth_out$new_TESSERAData_obj$X_list,
    coords_list = prep_synth_out$new_TESSERAData_obj$coords_list,
    k = 20,
    cov_type = "Mat",
    transform_z = TRUE,
    z_offset = 0.5,
    verbose = FALSE
  )

  # Parameter estimates
  alg_out$beta_hat <- alg_out$Beta
  alg_out$cov_param_hat <- alg_out$Theta
  # Fitted values
  brisc_preds <- BRISC::BRISC_prediction(alg_out, as.matrix(
    Reduce(rbind, prep_synth_out$new_TESSERAData_obj$coords_list)
  ), as.matrix(Reduce(
    rbind, prep_synth_out$new_TESSERAData_obj$X_list
  )))
  alg_out$predictions <- exp(as.vector(brisc_preds$prediction)) - 0.5
  # SE of parameter estimates
  brisc_boot <- BRISC::BRISC_bootstrap(alg_out, n_boot = 200)
  alg_out$beta_cov <- cov(brisc_boot$boot.Beta)
  alg_out$cov_param_cov <- cov(brisc_boot$boot.Theta)

  # Delete extra stuff that we don't need or want to save space
  keep_list <- c("predictions",
                 "beta_hat",
                 "cov_param_hat",
                 "beta_cov",
                 "cov_param_cov",
                 "time")
  alg_out <- alg_out[keep_list]
} else if (fit_model == "MGCV") {
  # MGCV
  alg_out <- mgcv_gam_wrapper(
    z_list = prep_synth_out$new_TESSERAData_obj$counts_list,
    X_list = prep_synth_out$new_TESSERAData_obj$X_list,
    coords_list = prep_synth_out$new_TESSERAData_obj$coords_list,
    model_family = "poisson",
    spline_k = -1,
    spline_basis = "tp",
    library_size_list = prep_synth_out$new_TESSERAData_obj$library_size_list
  )

  # Parameter estimates
  alg_out$beta_hat <- as.vector(coef(alg_out))[1:length(beta_hat)]
  # Fitted values
  alg_out$predictions <- as.vector(
    predict(alg_out, type = "response")
  )
  # SE of parameter estimates
  alg_out$beta_cov <- vcov(alg_out)[1:length(beta_hat), 1:length(beta_hat)]

  # Delete extra stuff that we don't need or want to save space
  keep_list <- c("predictions", "beta_hat", "beta_cov", "time")
  alg_out <- alg_out[keep_list]
} else if (fit_model == "LM") {
  # LM
  alg_out <- lm_wrapper(
    z_list = prep_synth_out$new_TESSERAData_obj$counts_list,
    X_list = prep_synth_out$new_TESSERAData_obj$X_list,
    transform_z = TRUE,
    z_offset = 0.5,
    library_size_list = prep_synth_out$new_TESSERAData_obj$library_size_list
  )

  # Parameter estimates
  alg_out$beta_hat <- coef(alg_out)
  # Fitted values
  alg_out$predictions <- exp(predict(alg_out)) - 0.5
  # SE of parameter estimates
  alg_out$beta_cov <- vcov(alg_out)

  # Delete extra stuff that we don't need or want to save space
  keep_list <- c("predictions", "beta_hat", "beta_cov", "time")
  alg_out <- alg_out[keep_list]
} else if (fit_model == "INLA_BYM2") {
  alg_out <- INLA_BYM2_wrapper(
    z_list = prep_synth_out$new_TESSERAData_obj$counts_list,
    X_list = prep_synth_out$new_TESSERAData_obj$X_list,
    W_list = prep_synth_out$new_TESSERAData_obj$W_list,
    model_family = "poisson",
    #z_offset = 0.5,
    num_threads = 1,
    compute_dic = TRUE,
    compute_waic = TRUE,
    compute_cpo = FALSE,
    verbose = TRUE,
    library_size_list = prep_synth_out$new_TESSERAData_obj$library_size_list
  )

  # Store stuff not already stired
  alg_out$predictions <- alg_out$predicted
  alg_out$gamma_hat_se <- alg_out$summary.hyperpar$sd[2]
  # Var(1/x) \approx (1/mean(x)^4) var(x)
  # https://stats.stackexchange.com/questions/41896/varx-is-known-how-to-calculate-var1-x
  alg_out$tau2_hat_se <-sqrt(alg_out$summary.hyperpar$sd[1]^2 / alg_out$summary.hyperpar$mean[1]^4)

  # Memory
  alg_out$all.hyper <- NA
  alg_out$model.matrix <- NA
  alg_out$.args <- NA
  alg_out$marginals.random <- NA
  alg_out$summary.linear.predictor <- NA
  alg_out$summary.fitted.values <- NA
  alg_out$summary.random <- NA
  alg_out$predicted <- NA
} else if (fit_model == "MCMC_CAR") {
  # MCMC CAR / SAR
  alg_out <- BRMS_CAR_SAR_wrapper(
    z_list = prep_synth_out$new_TESSERAData_obj$counts_list,
    X_list = prep_synth_out$new_TESSERAData_obj$X_list,
    W_list = prep_synth_out$new_TESSERAData_obj$W_list,
    model_type = gsub("MCMC_", "", fit_model),
    model_family = "poisson",
    chains = 4,
    iter = 2000,
    warmup = 200,
    thin = 2,
    cores = 1,
    library_size_list = prep_synth_out$new_TESSERAData_obj$library_size_list
  )

  # Store predictions to get residuals
  alg_out$predictions <- as.vector(predict(alg_out)[, 1])

  # Store estimates
  alg_out$beta_hat <- as.vector(brms::fixef(alg_out)[, 1])
  alg_out$gamma_hat <- as.data.frame(brms::posterior_summary(alg_out))["car", 1]
  alg_out$tau2_hat <- as.data.frame(brms::posterior_summary(alg_out))["sdcar", 1]^2

  # Store covariances
  alg_out$beta_cov <- vcov(alg_out)
  alg_out$gamma_hat_se <- as.data.frame(brms::posterior_summary(alg_out))["car", 2]
  tmp_se <- as.data.frame(brms::posterior_summary(alg_out))["sdcar", 2]
  alg_out$tau2_hat_se <- sqrt(4 * alg_out$tau2_hat^2 * tmp_se^2 + 2 * tmp_se^4)

  # Delete extra stuff that we don't need or want to save space
  keep_list <- c(
    "predictions",
    "beta_hat",
    "gamma_hat",
    "tau2_hat",
    "beta_cov",
    "gamma_hat_se",
    "tau2_hat_se",
    "time"
  )
  alg_out <- alg_out[keep_list]
} else if (fit_model == "MCMC_Leroux") {
  # MCMC Leroux
  alg_out <- CARBayes_Leroux_wrapper(
    z_list = prep_synth_out$new_TESSERAData_obj$counts_list,
    X_list = prep_synth_out$new_TESSERAData_obj$X_list,
    W_list = prep_synth_out$new_TESSERAData_obj$W_list,
    model_family = "poisson",
    chains = 4,
    iter = 2000,
    warmup = 200,
    thin = 2,
    cores = 1,
    library_size_list = prep_synth_out$new_TESSERAData_obj$library_size_list
  )

  # Store parameter estimates
  alg_out$beta_hat <- alg_out$summary.results[1:length(beta_hat), 1]
  alg_out$tau2_hat <- alg_out$summary.results[1 + length(beta_hat), 1]
  alg_out$gamma_hat <- alg_out$summary.results[2 + length(beta_hat), 1]

  # Store predictions
  alg_out$predictions <- alg_out$fitted.values

  # SE of parameter estimates
  alg_out$beta_cov <- cov(Reduce(rbind, alg_out$samples$beta))
  alg_out$gamma_hat_se <- sd(Reduce(c, alg_out$samples$rho))
  alg_out$tau2_hat_se <- sd(Reduce(c, alg_out$samples$tau2))

  # Delete extra stuff that we don't need or want to save space
  keep_list <- c(
    "predictions",
    "beta_hat",
    "gamma_hat",
    "tau2_hat",
    "beta_cov",
    "gamma_hat_se",
    "tau2_hat_se",
    "time"
  )
  alg_out <- alg_out[keep_list]
} else if (fit_model == "GLM_Poisson") {
  # GLM Poisson
  alg_out <- glm_wrapper(
    z_list = prep_synth_out$new_TESSERAData_obj$counts_list,
    X_list = prep_synth_out$new_TESSERAData_obj$X_list,
    library_size_list = prep_synth_out$new_TESSERAData_obj$library_size_list
  )

  # Parameter estimates
  alg_out$beta_hat <- coef(alg_out)
  # Fitted values
  alg_out$predictions <- as.vector(
    predict(alg_out, type = "response")
  )
  # SE of parameter estimates
  alg_out$beta_cov <- vcov(alg_out)

  # Delete extra stuff that we don't need or want to save space
  keep_list <- c("predictions", "beta_hat", "beta_cov", "time")
  alg_out <- alg_out[keep_list]
} else if (fit_model == "GLM_NB") {
  # GLM NB
  alg_out <- glm_nb_wrapper(
    z_list = prep_synth_out$new_TESSERAData_obj$counts_list,
    X_list = prep_synth_out$new_TESSERAData_obj$X_list,
    library_size_list = prep_synth_out$new_TESSERAData_obj$library_size_list
  )

  # Parameter estimates
  alg_out$beta_hat <- coef(alg_out)
  alg_out$tau2_hat <- alg_out$theta
  # Fitted values
  alg_out$predictions <- as.vector(
    predict(alg_out, type = "response")
  )
  # SE of parameter estimates
  alg_out$beta_cov <- vcov(alg_out)
  alg_out$tau2_hat_se <- alg_out$SE.theta

  # Delete extra stuff that we don't need or want to save space
  keep_list <- c("predictions",
                 "beta_hat",
                 "beta_cov",
                 "tau2_hat",
                 "tau2_hat_se",
                 "time")
  alg_out <- alg_out[keep_list]
}

## Save to file
if (fit_model == "spNNGP" || fit_model == "BRISC") {
  outfile_base <- paste0(
    "poisECM_synth_natgen_kidney_",
    "datamodel_",
    data_model,
    "_fitmodel_",
    fit_model,
    "_kernel_",
    kernel_type,
    "_trial_",
    trial_idx,
    "_gene_",
    gene_name,
    "_sample_",
    sample_name
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
    gene_name,
    "_sample_",
    sample_name
  )
}

MSE_total <- mean((as.vector(
  Reduce(c, prep_synth_out$new_TESSERAData_obj$counts_list)
) - alg_out$predictions)^2)
Mean_counts2 <- mean((as.vector(
  Reduce(c, prep_synth_out$new_TESSERAData_obj$counts_list)
))^2)
end_idx <- as.vector(cumsum(
  lapply(prep_synth_out$new_TESSERAData_obj$counts_list, length)
))
MSE_sample <- rep(0, length(end_idx))
summary_df <- list()
for (idx in 1:length(prep_synth_out$new_TESSERAData_obj$counts_list)) {
  start_idx <- end_idx[idx] - length(prep_synth_out$new_TESSERAData_obj$counts_list[[idx]]) + 1

  true_local <- as.vector(prep_synth_out$new_TESSERAData_obj$counts_list[[idx]])
  pred_local <- alg_out$predictions[start_idx:end_idx]
  resid_local <- true_local - pred_local
  MSE_sample[idx] <- mean((resid_local)^2)

  if ("gamma_hat" %in% names(alg_out)) {
    fitted_gamma_hat <- alg_out$gamma_hat[idx]
  } else {
    fitted_gamma_hat <- NA
  }
  if ("tau2_hat" %in% names(alg_out)) {
    fitted_tau2_hat <- alg_out$tau2_hat[idx]
  } else {
    fitted_tau2_hat <- NA
  }
  if ("cov_param_hat" %in% names(alg_out)) {
    nugget_hat <- alg_out$cov_param_hat[1]
    sill_hat <- alg_out$cov_param_hat[2]
    range_hat <- alg_out$cov_param_hat[3]
    smoothness_hat <- alg_out$cov_param_hat[4]
  } else {
    nugget_hat <- NA
    sill_hat <- NA
    range_hat <- NA
    smoothness_hat <- NA
  }
  Moran_counts <- calc_moran(
    true_local,
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 1],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 2]
  )[1]
  Moran_predictions <- calc_moran(
    pred_local,
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 1],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 2]
  )[1]
  Moran_residuals <- calc_moran(
    resid_local,
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 1],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 2]
  )[1]
  Moran_Xbeta <- calc_moran(
    as.matrix(prep_synth_out$new_TESSERAData_obj$X_list[[idx]]) %*% alg_out$beta_hat[1:length(beta_hat)],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 1],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 2]
  )[1]
  Moran_librarysize <- calc_moran(
    prep_synth_out$new_TESSERAData_obj$library_size_list[[idx]],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 1],
    prep_synth_out$new_TESSERAData_obj$coords_list[[idx]][, 2]
  )[1]
  if ("performanceSummary" %in% names(alg_out)) {
    Moran_phi <- alg_out$performanceSummary$Moran_phi[idx]
    Moran_eta <- alg_out$performanceSummary$Moran_eta[idx]
    Moran_theta <- alg_out$performanceSummary$Moran_theta[idx]
  } else {
    Moran_theta <- NA
    Moran_eta <- NA
    Moran_phi <- NA
  }
  tmp_df <- data.frame(
    gene = gene_name,
    data_model = data_model,
    fit_model = fit_model,
    sample = sample_name,
    n_cells = length(true_local),
    gamma_hat = fitted_gamma_hat,
    tau2_hat = fitted_tau2_hat,
    kernel_type = kernel_type,
    nugget_hat = nugget_hat,
    sill_hat = sill_hat,
    range_hat = range_hat,
    smoothness_hat = smoothness_hat,
    MSE_counts_sample = MSE_sample[idx],
    MSE_counts_total = MSE_total,
    Mean_counts2_sample = mean(true_local^2),
    Mean_counts2_total = Mean_counts2,
    Moran_counts = Moran_counts,
    Moran_predictions = Moran_predictions,
    Moran_residuals = Moran_residuals,
    Moran_phi = Moran_phi,
    Moran_eta = Moran_eta,
    Moran_Xbeta = Moran_Xbeta,
    Moran_theta = Moran_theta,
    Moran_librarysize = Moran_librarysize,
    beta_SE = sum((alg_out$beta_hat - beta_hat)^2),
    beta_norm2 = sum((beta_hat)^2),
    true_gamma = gamma_hat,
    true_tau2 = tau2_hat,
    true_nugget = cov_params_hat[1],
    true_sill = cov_params_hat[2],
    true_range = cov_params_hat[3],
    true_smoothness = cov_params_hat[4]
  )
  summary_df[[idx]] <- tmp_df
}
summary_df <- dplyr::bind_rows(summary_df)
alg_out$summary_df <- summary_df


# Extra information to store for later
alg_out$sample_name <- sample_name
alg_out$sim_settings <- sim_df[run_idx, ]
alg_out$GLOBAL_SEED <- GLOBAL_SEED
if (fit_model == "spNNGP" || fit_model == "BRISC") {
  alg_out$fit_model_kernel_type <- kernel_type
} else{
  alg_out$fit_model_kernel_type <- NA
}
if (data_model == "spNNGP" || data_model == "BRISC") {
  alg_out$data_model_kernel_type <- kernel_type
} else{
  alg_out$data_model_kernel_type <- NA
}

outfile_base <- file.path(out_path, outfile_base)
saveRDS(alg_out,
        file = paste0(outfile_base, ".rds"),
        compress = TRUE)

## End
print("Done.")
