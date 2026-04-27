#' Create and run TESSERA and competing methods for timing results from synthetic data based on real-world kidney data
#' @author Florica Constantine


## Scripts
source("TESSERA_comparison_methods.R")

# Libraries: Custom
library(TESSERA)

# Libraries: CRAN
library(dplyr)
library(Matrix)
library(fitdistrplus)

# Paths
path <- "/scratch/users/spatialseq/natgen_kidney/"
processed_data_path <- file.path(path, "processed")
real_out_path <- file.path(path, "real_output_onlyinteract")
out_path <- file.path(path, "timing")

## Real data

chosen_gene <- "IGHG2"

# Read in beta values for gene at 75th percentile for highest variance
real_data_results <- readRDS(
  paste0(real_out_path, "poisECM_real_natgen_kidney_model_Leroux_gene_", chosen_gene, ".rds")
)
real_data_input <- readRDS(
  processed_data_path, "prepData_lattice_onlyinteract.rds"
)

## Common Parameters

GLOBAL_SEED <- 2025 # Common seed
# Generate data from a Leroux model: the actual values don't really matter
data_model <- "Leroux"
gene_name <- paste0("synthetic_", chosen_gene)

# TESSERA Parameters
em_iters <- 30
opt_iters <- 5
em_min_iters <- em_iters
em_tol <- 1e-3
em_stopping <- NULL
beta_init <- "glm"
gamma_init <- "moran"
tau2_init <- "var"
# spNNGP Parameters, shared with BRISC
cov_init <- "BRISC"
cov_fit_method <- "BRISC"
cov_type <- "Mat"
nngp_k <- 20
verbose <- FALSE

# MCMC Parameters
mcmc_chains <- 4
mcmc_iters <- 2000
mcmc_thin <- 2
mcmc_warmup <- 200

# Log-Transform Parameters
z_offset <- 0.5

# MGCV Parameters
spline_k <- -1
spline_basis <- "tp"


## Simulation Parameters

# What models to run
fit_model_list <- c(
  "CAR",
  "SAR",
  "Leroux",
  "spNNGP",
  "MGCV",
  "Poisson_GLM",
  "NB_GLM",
  "LM",
  "BRISC",
  "INLA_BYM2",
  "MCMC_CAR",
  "MCMC_Leroux"
)
# What numbers of samples to run
n_samples_list <- c(1, 6, 12, 24)
# What numbers of cells per sample to run
n_cells_per_sample_list <- c(500, 1000, 2000, 4000)
# Dimensionality of covariate matrix
beta_dim_list <- c(25, 50, 100)
# Number of trials
trial_list <- c(1:5)

# Create dataframe with runs
sim_df <- expand.grid(n_samples_list,
                      n_cells_per_sample_list,
                      beta_dim_list,
                      trial_list,
		      fit_model_list)
colnames(sim_df) <- c("n_samples", "n_cells", "beta_dim", "trial_idx", "fit_model")

# Subset
sim_df <- sim_df[which(sim_df$fit_model %in% c("INLA_BYM2")), ]


## Select Simulation Parameters

# Define so linter doesn't warn later
# run_idx goes from 1 to nrow(sim_df)
# but most clusters/slurm setups can't handle large array jobs (e.g., a max of 5k)
run_idx <- NA
offset <- 0
args <- commandArgs(trailingOnly = TRUE)
print("Original run_idx and offset")
if (0 == length(args)) {
  stop("No arguments supplied: ERROR")
} else {
  for (idx in 1:length(args)) {
    eval(parse(text = args[idx]))
  }
  print(run_idx)
  print(offset)
}

# Shift run_idx
run_idx <- run_idx + offset
if ((1 > run_idx) | (nrow(sim_df) < run_idx)) {
  stop("run_idx must be between 1 and #rows in sim_df")
}

fit_model <- sim_df[run_idx, ]$fit_model
n_samples <- sim_df[run_idx, ]$n_samples
n_cells_per_sample <- sim_df[run_idx, ]$n_cells
beta_dim <- sim_df[run_idx, ]$beta_dim
trial_idx <- sim_df[run_idx, ]$trial_idx
print(sim_df[run_idx, ])


## Derived parameters

# Set seed for data generation
DATA_SEED <- GLOBAL_SEED + trial_idx
set.seed(DATA_SEED)

# Sample from real data fits
real_gamma_tau2_df <- data.frame(gamma = real_data_results$gamma_hat, tau2 =
                                   real_data_results$tau2_hat)
real_gamma_tau2_df$gamma <- pmin(real_gamma_tau2_df$gamma, 0.95)
g_t_rows <- sample(1:nrow(real_gamma_tau2_df),
                   size = n_samples,
                   replace = TRUE)
gamma_true <- real_gamma_tau2_df$gamma[g_t_rows]
tau2_true <- real_gamma_tau2_df$tau2[g_t_rows]


## Generate Data

# https://math.stackexchange.com/questions/1995011/within-a-unit-square-given-n-random-uniform-points-what-is-the-average-distanc
# See comment and https://www3.nd.edu/~mhaenggi/pubs/tvt09.pdf
# Points in the square are approximately a Poisson Process
# The number of points in a ball of radius r is Poisson, with parameter \pi r^2 n,
# for n points
# If we want on average k points in a ball of radius r, for a given n,
# We set k = \pi r^2 n, or, \sqrt{k / (\pi n)} = r
n_nb_lattice <- 5
nb_dist <- sqrt(n_nb_lattice / (pi * n_cells_per_sample)) # Distance to be considered a neighbor

# Generate covariates
# Sample effects, two diseases, and rest are cell types
# Even split for diseases and cell types
covariates_list <- list()
disease_A_current_idx <- 0
disease_B_current_idx <- 0
covariate_levels <- c(
  LETTERS,
  paste0(LETTERS, LETTERS),
  paste0(LETTERS, LETTERS, LETTERS),
  paste0(LETTERS, LETTERS, LETTERS, LETTERS)
)

set.seed(DATA_SEED)

for (idx in 1:n_samples) {
  # Generate a more realistic set of covariates
  # E.g., 5 diseases and 5 cell types and 5 groups ...
  if (0 == idx %% 2) {
    disease_level <- "A"
    disease_A_current_idx <- disease_A_current_idx + 1
    sample_idx_local <- disease_A_current_idx
  } else {
    disease_level <- "B"
    disease_B_current_idx <- disease_B_current_idx + 1
    sample_idx_local <- disease_B_current_idx
  }
  # For a single sample, just simulate cell types
  if (1 < n_samples) {
    covariates <- data.frame(
      disease = rep(disease_level, n_cells_per_sample),
      sample = factor(rep(
        sample_idx_local, n_cells_per_sample
      ))
    )
    # 2 for disease and 2 * (n_samples / 2 - 1) for samples
    n_cols_used_up <- 2 + n_samples - 2
    n_cols_left <- beta_dim - n_cols_used_up
    if (0 < n_cols_left) {
      source_levels <- covariate_levels[1:(n_cols_left + 1)]
      covariates$cell_type <- as.character(sample(source_levels, n_cells_per_sample, replace = TRUE))
    }
  } else {
    source_levels <- covariate_levels[1:beta_dim]
    covariates <- data.frame(cell_type = as.character(sample(
      source_levels, n_cells_per_sample, replace = TRUE
    )))
  }

  covariates_list[[idx]] <- covariates
}
covariates_list <- dplyr::bind_rows(covariates_list)
if (2 < n_samples) {
  X_mat <- as.matrix(model.matrix(~ 0 + disease + sample:disease + cell_type, data = covariates_list))
} else if (2 == n_samples) {
  X_mat <- as.matrix(model.matrix(~ 0 + disease + cell_type, data = covariates_list))
} else {
  X_mat <- as.matrix(model.matrix(~ 0 + cell_type, data = covariates_list))
}
stopifnot(ncol(X_mat) == beta_dim)
stopifnot(ncol(X_mat) == Matrix::rankMatrix(X_mat)[1])

# Sample beta: beta drawn to have same mean and variance as real data

set.seed(DATA_SEED)

beta_true <- rep(0, beta_dim)
names(beta_true) <- colnames(X_mat)
disease_A_lb <- "Control"
disease_B_lb <- "DKD"
# Sample cell type effects by averaging over all cell type effects
ct_cols_real <- which(grepl("celltype", names(real_data_results$beta_hat)))
ct_cols_sim <- which(grepl("cell_type", names(beta_true)))
beta_true[ct_cols_sim] <- rnorm(
  length(ct_cols_sim),
  mean = mean(real_data_results$beta_hat[ct_cols_real]),
  sd = sd(real_data_results$beta_hat[ct_cols_real])
)
# Sample disease A effect by averaging over control effects
ct_cols_real <- which(grepl(disease_A_lb, names(real_data_results$beta_hat)))
ct_cols_sim <- which("diseaseA" == names(beta_true))
beta_true[ct_cols_sim] <- rnorm(
  length(ct_cols_sim),
  mean = mean(real_data_results$beta_hat[ct_cols_real]),
  sd = sd(real_data_results$beta_hat[ct_cols_real])
)
# Sample disease B effect by averaging over DKD effects
ct_cols_real <- which(grepl(disease_B_lb, names(real_data_results$beta_hat)))
ct_cols_sim <- which("diseaseB" == names(beta_true))
beta_true[ct_cols_sim] <- rnorm(
  length(ct_cols_sim),
  mean = mean(real_data_results$beta_hat[ct_cols_real]),
  sd = sd(real_data_results$beta_hat[ct_cols_real])
)
# Sample disease A sample effects
ct_cols_real <- which(grepl(disease_A_lb, names(real_data_results$beta_hat)) &
                        grepl("_ST", names(real_data_results$beta_hat)))
ct_cols_sim <- which(grepl("diseaseA:sample", names(beta_true)))
beta_true[ct_cols_sim] <- rnorm(
  length(ct_cols_sim),
  mean = mean(real_data_results$beta_hat[ct_cols_real]),
  sd = sd(real_data_results$beta_hat[ct_cols_real])
)
# Sample disease B sample effects
ct_cols_real <- which(grepl(disease_B_lb, names(real_data_results$beta_hat)) &
                        grepl("_ST", names(real_data_results$beta_hat)))
ct_cols_sim <- which(grepl("diseaseB:sample", names(beta_true)))
beta_true[ct_cols_sim] <- rnorm(
  length(ct_cols_sim),
  mean = mean(real_data_results$beta_hat[ct_cols_real]),
  sd = sd(real_data_results$beta_hat[ct_cols_real])
)

# Fit distribution to library sizes
lib_size_dist <- fitdistrplus::fitdist(Reduce(c, real_data_input$library_size_list), "nbinom")

# Sample data
data_list <- list()
library_size_list <- list()
set.seed(DATA_SEED)
for (idx in 1:n_samples) {
  if (1 < n_samples) {
    idx_start <- 1 + (idx - 1) * n_cells_per_sample
    idx_end <- idx * n_cells_per_sample
    X_mat_local <- X_mat[idx_start:idx_end, ]
  } else {
    X_mat_local <- X_mat
  }

  # Sample library size
  lib_size_local <- rnbinom(n_cells_per_sample,
                            size = lib_size_dist$estimate["size"],
                            mu = lib_size_dist$estimate["mu"])
  library_size_list[[idx]] <- lib_size_local

  data_list[[idx]] <- TESSERA::generate_data_one_area(
    n_cells_per_sample,
    nb_dist,
    data_model,
    beta_true,
    gamma_true[idx],
    tau2_true[idx],
    X = X_mat_local,
    library_size = lib_size_local
  )
}
# Memory
rm(covariates,
   covariates_list,
   source_levels,
   X_mat,
   X_mat_local,
   lib_size_local)
rm(real_data_input, real_data_results)
gc()

# Convert to a format amenable to the prepData function
count_matrix <- matrix(Reduce(c, lapply(data_list, function (x) {
  x$z
})), nrow = 1)
rownames(count_matrix) <- c(gene_name)
colnames(count_matrix) <- 1:ncol(count_matrix)
meta_data <- data.frame(sample = Reduce(c, lapply(1:n_samples, function (x) {
  rep(x, length(data_list[[x]]$z))
})))
rownames(meta_data) <- colnames(count_matrix)
X_mat <- Reduce(rbind, lapply(data_list, function (x) {
  x$X
}))
rownames(X_mat) <- colnames(count_matrix)
coord_mat <- Reduce(rbind, lapply(data_list, function (x) {
  cbind(x$x_coords, x$y_coords)
}))
rownames(coord_mat) <- colnames(count_matrix)
adj_mat <- Matrix::bdiag(lapply(data_list, function (x) {
  x$W
}))
rownames(adj_mat) <- colnames(count_matrix)
colnames(adj_mat) <- colnames(count_matrix)
data_obj <- TESSERA::prep_data(
  count_matrix,
  meta_data,
  "sample",
  design_formula = NULL,
  design_mat = X_mat,
  coord_data = coord_mat,
  D_THRESH = NULL,
  k_search = NULL,
  adj_mat = adj_mat,
  model_type = "ALL"
)
# Add in library size list
names(library_size_list) <- names(data_obj$counts_list)
for (idx in 1:length(data_obj$counts_list)) {
  names(library_size_list[[idx]]) <- colnames(data_obj$counts_list[[idx]])
}
data_obj$library_size_list <- library_size_list

# Memory
rm(data_list,
   count_matrix,
   meta_data,
   X_mat,
   coord_mat,
   adj_mat,
   library_size_list)
gc()


## Run algorithm
set.seed(DATA_SEED)

if (("CAR" == fit_model) ||
    ("SAR" == fit_model) || ("Leroux" == fit_model)) {
  invisible(capture.output(
    alg_out <- TESSERA::TESSERA_lattice(
      TESSERAData_obj = data_obj,
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
      verbose = verbose
    )
  ))
} else if ("spNNGP" == fit_model) {
  invisible(capture.output(
    alg_out <- TESSERA::TESSERA_spNNGP(
      TESSERAData_obj = data_obj,
      gene_name = gene_name,
      cov_type = cov_type,
      nngp_k = nngp_k,
      em_iters = em_iters,
      opt_iters = opt_iters,
      em_min_iters = em_min_iters,
      em_tol = em_tol,
      em_stopping = em_stopping,
      beta_init = beta_init,
      cov_fit_method = cov_fit_method,
      cov_init = cov_init,
      verbose = verbose
    )
  ))
} else if ("INLA_BYM2" == fit_model) {
  alg_out <- INLA_BYM2_wrapper(
    z_list = data_obj$counts_list,
    X_list = data_obj$X_list,
    W_list = data_obj$W_list,
    model_family = "poisson",
    #z_offset = 0.5,
    num_threads = 1,
    compute_dic = TRUE,
    compute_waic = TRUE,
    compute_cpo = FALSE,
    verbose = FALSE,
    library_size_list = data_obj$library_size_list
  )
} else if ("MGCV" == fit_model) {
  alg_out <- mgcv_gam_wrapper(
    z_list = data_obj$counts_list,
    X_list = data_obj$X_list,
    coords_list = data_obj$coords_list,
    model_family = "poisson",
    spline_k = spline_k,
    spline_basis = spline_basis,
    library_size_list = data_obj$library_size_list
  )

  # Parameter estimates
  alg_out$beta_hat <- as.vector(coef(alg_out))[1:length(beta_true)]
} else if ("Poisson_GLM" == fit_model) {
  alg_out <- glm_wrapper(
    z_list = data_obj$counts_list,
    X_list = data_obj$X_list,
    library_size_list = data_obj$library_size_list
  )

  # Parameter estimates
  alg_out$beta_hat <- coef(alg_out)
} else if ("NB_GLM" == fit_model) {
  alg_out <- glm_nb_wrapper(
    z_list = data_obj$counts_list,
    X_list = data_obj$X_list,
    library_size_list = data_obj$library_size_list
  )

  # Parameter estimates
  alg_out$beta_hat <- coef(alg_out)
} else if ("LM" == fit_model) {
  alg_out <- lm_wrapper(
    z_list = data_obj$counts_list,
    X_list = data_obj$X_list,
    transform_z = TRUE,
    z_offset = z_offset
  )

  # Parameter estimates
  alg_out$beta_hat <- coef(alg_out)
} else if ("BRISC" == fit_model) {
  alg_out <- BRISC_wrapper(
    z_list = data_obj$counts_list,
    X_list = data_obj$X_list,
    coords_list = data_obj$coords_list,
    k = nngp_k,
    cov_type = cov_type,
    transform_z = TRUE,
    z_offset = z_offset,
    verbose = verbose
  )

  # Parameter estimates
  alg_out$beta_hat <- alg_out$Beta
} else if (("MCMC_CAR" == fit_model) || ("MCMC_SAR" == fit_model)) {
  alg_out <- BRMS_CAR_SAR_wrapper(
    z_list = data_obj$counts_list,
    X_list = data_obj$X_list,
    W_list = data_obj$W_list,
    model_type = gsub("MCMC_", "", fit_model),
    model_family = "poisson",
    chains = mcmc_chains,
    iter = mcmc_iters,
    warmup = mcmc_warmup,
    thin = mcmc_thin,
    cores = 1,
    library_size_list = data_obj$library_size_list
  )

  # Parameter estimates
  alg_out$beta_hat <- as.vector(brms::fixef(alg_out)[, 1])
} else if ("MCMC_Leroux" == fit_model) {
  alg_out <- CARBayes_Leroux_wrapper(
    z_list = data_obj$counts_list,
    X_list = data_obj$X_list,
    W_list = data_obj$W_list,
    model_family = "poisson",
    chains = mcmc_chains,
    iter = mcmc_iters,
    warmup = mcmc_warmup,
    thin = mcmc_thin,
    cores = 1,
    library_size_list = data_obj$library_size_list
  )

  # Store parameter estimates
  alg_out$beta_hat <- alg_out$summary.results[1:length(beta_true), 1]
}


## Save and end
alg_run_time <- as.numeric(alg_out$time, units = "secs")

alg_out_list <- list(run_info = sim_df[run_idx, ],
     alg_run_time = alg_run_time,
     beta_hat = alg_out$beta_hat,
     beta_true = beta_true,
     gamma_true = gamma_true,
     tau2_true = tau2_true,
     data_model = data_model,
     gene_name = gene_name
     )

## Save to file
if (fit_model == "spNNGP" | fit_model == "BRISC") {
  outfile_base <- paste0(
    "poisECM_timing_",
    "data_model_",
    data_model,
    "_fitmodel_",
    fit_model,
    "_kernel_",
    "Mat",
    "_nsamp_",
    n_samples,
    "_ncells_",
    n_cells_per_sample,
    "_betadim_",
    beta_dim,
    "_trial_",
    trial_idx
  )
} else{
  outfile_base <- paste0(
    "poisECM_timing_",
    "data_model_",
    data_model,
    "_fitmodel_",
    fit_model,
    "_nsamp_",
    n_samples,
    "_ncells_",
    n_cells_per_sample,
    "_betadim_",
    beta_dim,
    "_trial_",
    trial_idx
  )
}
outfile_base <- paste0(out_path, outfile_base)
saveRDS(alg_out_list,
        file = paste0(outfile_base, ".rds"),
        compress = TRUE)

## End
print("Done.")
