#' Conduct a power analysis using pseudobulk methods on synthetic data based on real-world kidney data
#' @author Florica Constantine


## Libraries

# General data manipulation
library(dplyr)
library(Matrix)

# Generating synthetic data
library(TESSERA)

# Fitting methods
library(DESeq2)
library(limma)
library(edgeR)


## Scripts: For pseudobulking data

source("aggregate_data.R")


## Paths

path <- "/scratch/users/spatialseq/natgen_kidney/"
processed_data_path <- file.path(path, "processed")
real_out_path <- file.path(path, "real_output_onlyinteract")
out_path <- file.path(path, "power_contrasts")


## Simulation setup

# Data generation
datagen_model_list <- c("spNNGP")
# Gene
gene_list <- c("IGHG2")
# Number of trials
n_trials <- 100
# Create new beta values
pos_betas <- c(
  seq(from = 0.01, to = 0.1, by = 0.01),
  seq(from = 0.15, to = 0.3, by = 0.05),
  seq(from = 0.4, to = 0.5, by = 0.1),
  seq(from = 0.6, to = 1, by = 0.2),
  seq(from = 1.25, to = 1.5, by = 0.25),
  2
)
neg_betas <- -1 * rev(pos_betas)
beta_list <- c(neg_betas, 0, pos_betas)

# Which covariates to modify
covariate_list <- c(
  "GroupDKD:celltypePodo",
  "GroupDKD:celltypeIC",
  "GroupDKD:celltypeGS_Stromal",
  "GroupDKD:celltypeC_TAL"
)

sim_df <- expand.grid(1:n_trials,
                      datagen_model_list,
                      gene_list,
                      beta_list,
                      covariate_list)
colnames(sim_df) <- c("trial", "datagen_model", "gene", "beta", "covariate")


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
gene_name <- sim_df[run_idx, ]$gene
trial_idx <- sim_df[run_idx, ]$trial
beta_val <- sim_df[run_idx, ]$beta
covariate_name <- sim_df[run_idx, ]$covariate
print(sim_df[run_idx, ])

# Change DKD to Control in string
covariate_pair_function <- function (x) {
  return (gsub("DKD", "Control", x))
}

covariate_contrast <- covariate_pair_function(covariate_name)
print(covariate_contrast)


## Global parameters

GLOBAL_SEED <- 2025
DATA_SEED <- GLOBAL_SEED + trial_idx
set.seed(DATA_SEED)


## TESSERA parameters

nngp_k <- 20
max_gamma <- 0.999


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
  # Unused parameters
  gamma_hat <- NA
  tau2_hat <- NA
} else if (data_model %in% c("Leroux", "CAR", "SAR")) {
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

  # Unused parameters
  cov_params_hat <- NA
}


## Load in prepData object for data

if (data_model == "spNNGP") {
  prep_real_out <- readRDS(file.path(processed_data_path, "prepData_spNNGP_onlyinteract.rds"))
} else{
  prep_real_out <- readRDS(file.path(processed_data_path, "prepData_lattice_onlyinteract.rds"))
}


## Set up beta and contrasts

# Create new beta vector
beta_hat_power <- beta_hat
beta_baseline_idx <- which(names(beta_hat_power) == covariate_contrast)
beta_change_idx <- which(names(beta_hat_power) == covariate_name)
beta_baseline_val <- mean(beta_hat_power[which(grepl("celltype", names(beta_hat_power)))])
beta_hat_power[beta_baseline_idx] <- beta_baseline_val
beta_hat_power[beta_change_idx] <- beta_baseline_val + beta_val
print(beta_hat[c(beta_baseline_idx, beta_change_idx)])
print(beta_hat_power[c(beta_baseline_idx, beta_change_idx)])

# Contrast matrix
contrast_mat <- rep(0, length(beta_hat_power))
contrast_mat[beta_baseline_idx] <- 1
contrast_mat[beta_change_idx] <- -1
contrast_mat <- matrix(contrast_mat, nrow = 1, ncol = length(contrast_mat))
colnames(contrast_mat) <- names(beta_hat_power)
rownames(contrast_mat)[1] <- paste0(covariate_contrast, "-", covariate_name)


## Create new prepData object with counts sampled using new beta vector

if (data_model == "spNNGP") {
  prep_synth_out <- TESSERA::prep_synth_data(
    TESSERAData_obj = prep_real_out,
    gene_list = gene_name,
    data_gen_model = data_model,
    cov_params = cov_params_hat,
    cov_type = kernel_type,
    nngp_k = nngp_k,
    beta_true = beta_hat_power
  )
} else if (data_model %in% c("Leroux", "CAR", "SAR")) {
  prep_synth_out <- TESSERA::prep_synth_data(
    TESSERAData_obj = prep_real_out,
    gene_list = gene_name,
    data_gen_model = data_model,
    tau2_true = tau2_hat,
    gamma_true = pmin(gamma_hat, max_gamma),
    beta_true = beta_hat_power
  )
}

# Sometimes spNNGP generation leads to massive outliers in the sampled counts
# E.g., there'll be one sample with one or two counts that are 10^6 or greater
# So, we resample until this no longer happens
if (data_model == "spNNGP") {
  for (s_idx in 1:length(prep_synth_out$new_TESSERAData_obj$counts_list)) {
    loop_counter <- 0
    while (max(prep_synth_out$new_TESSERAData_obj$counts_list[[s_idx]]) > 1e4) {
      cat("Resampling", s_idx, "\n")
      tmp <- TESSERA::sample_Poisson_spNNGP(
        cov_type = kernel_type,
        X = prep_real_out$X_list[[s_idx]],
        library_size = prep_real_out$library_size_list[[s_idx]],
        coords = prep_real_out$coords_list[[s_idx]],
        cov_params = cov_params_hat[s_idx, ],
        nngp_k = nngp_k,
        beta_true = beta_hat_power
      )$z
      prep_synth_out$new_TESSERAData_obj$counts_list[[s_idx]][1, ] <- tmp
      loop_counter <- loop_counter + 1
      if (10 < loop_counter) {
        cat("Failed to sample, breaking\n")
        break
      }
    }
  }
}

## Make pseudobulk data

# Edit the counts matrix so all genes/edited counts are in the same matrix
# Replace original gene counts with resampled counts

gene_idx <- which(rownames(prep_real_out$counts_list[[1]]) == gene_name)
prep_new <- prep_real_out
for (idx in 1:length(prep_new$counts_list)) {
  prep_new$counts_list[[idx]][gene_idx, ] <- prep_synth_out$new_TESSERAData_obj$counts_list[[idx]][1, ]
}

# Clean up (memory)
rm(prep_real_out)
rm(prep_synth_out)
gc()

# Aggregate to single data matrix/frame
count_matrix <- Reduce(cbind, prep_new$counts_list)[1:3000, ]
meta_data <- Reduce(rbind, prep_new$covariates_list)

# Clean up (memory)
rm(prep_new)
rm(real_alg_out)
gc()

# Call function to create pseudobulk and bulk data
sample_col <- "orig.ident"
celltype_col <- "celltype"
cell_id_col <- "cell_id_metadata"
group_col <- "Group"
pb_data <- make_pseudobulk_bulk(
  count_matrix,
  meta_data,
  individual_colname = sample_col,
  cell_type_colname = celltype_col,
  cell_id_colname = cell_id_col
)
pseudobulk_meta <- pb_data$pseudobulk_meta
pseudobulk_meta <- pseudobulk_meta[, c("pb_id", sample_col, celltype_col, group_col)]
pseudobulk_meta <- dplyr::distinct(pseudobulk_meta)
rownames(pseudobulk_meta) = pseudobulk_meta$pb_id
pseudobulk_meta <- pseudobulk_meta[, -which(colnames(pseudobulk_meta) == "pb_id")]

control_ids <- unique(pseudobulk_meta[, sample_col][pseudobulk_meta[, group_col] == "Control"])
DKD_ids <- unique(pseudobulk_meta[, sample_col][pseudobulk_meta[, group_col] == "DKD"])
HKD_ids <- unique(pseudobulk_meta[, sample_col][pseudobulk_meta[, group_col] == "HKD"])
new_id <- rep(0, nrow(pseudobulk_meta))
for (idx in 1:length(control_ids)) {
  new_id[pseudobulk_meta[, sample_col] == control_ids[idx]] <- idx
}
for (idx in 1:length(DKD_ids)) {
  new_id[pseudobulk_meta[, sample_col] == DKD_ids[idx]] <- idx
}
for (idx in 1:length(HKD_ids)) {
  new_id[pseudobulk_meta[, sample_col] == HKD_ids[idx]] <- idx
}
pseudobulk_meta$nested_id <- factor(new_id)
new_nested_sample_col <- "nested_id"

# Design formula
design_formula <- as.formula(paste0("~ 0 + ", paste(
  group_col,
  c(celltype_col, new_nested_sample_col),
  sep = ":",
  collapse = " + "
)))
print(design_formula)

# Design matrix
design_mat <- model.matrix(design_formula, pseudobulk_meta)
# Drop zero columns
design_mat <- design_mat[, -which(colSums(design_mat) == 0)]
stopifnot(Matrix::rankMatrix(design_mat)[1] == ncol(design_mat))

# Rename sample/patient effects for clarity later
for (idx in 2:length(DKD_ids)) {
  colnames(design_mat) <- gsub(
    paste0("GroupDKD:nested_id", idx),
    paste0("GroupDKD:", DKD_ids[idx]),
    colnames(design_mat)
  )
}
for (idx in 2:length(HKD_ids)) {
  colnames(design_mat) <- gsub(
    paste0("GroupHKD:nested_id", idx),
    paste0("GroupHKD:", HKD_ids[idx]),
    colnames(design_mat)
  )
}
for (idx in 2:length(control_ids)) {
  colnames(design_mat) <- gsub(
    paste0("GroupControl:nested_id", idx),
    paste0("GroupControl:", control_ids[idx]),
    colnames(design_mat)
  )
}

stopifnot(all.equal(colnames(contrast_mat), colnames(design_mat)))

# Clean up (memory)
rm(count_matrix)
rm(meta_data)
gc()


## Run DESeq2

t0 <- Sys.time()

# Put data into a DESeq2 object
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = pb_data$pseudobulk_data,
  colData = pseudobulk_meta,
  design = design_mat
)
# Run DESeq2
dds <- DESeq2::DESeq(dds)

# Get test statistics
res_deseq2 <- DESeq2::results(dds, contrast = contrast_mat[1, ], parallel = FALSE)
res_deseq2 <- as.data.frame(res_deseq2)
# Subset to synthetic gene of interest
res_deseq2 <- res_deseq2[gene_idx, ]

t1 <- Sys.time()

# Create dataframe with results
beta_hat_deseq2 <- coef(dds)[gene_idx, ]
res_deseq2 <- data.frame(
  datagen_model = data_model,
  trial = trial_idx,
  changed_covariate = covariate_name,
  baseline_contrast = covariate_contrast,
  true_contrast_value = beta_val,
  baseline_coef_value = beta_hat_power[beta_baseline_idx],
  changed_coef_value = beta_hat_power[beta_change_idx],
  estimated_baseline_coef = beta_hat_deseq2[beta_baseline_idx],
  estimated_changed_coef = beta_hat_deseq2[beta_change_idx],
  gene = gene_name,
  fit_model = "DESeq2",
  kernel_type = NA,
  contrast_indices = paste0(beta_baseline_idx, "_", beta_change_idx),
  contrast_val = res_deseq2$log2FoldChange,
  contrast_se = res_deseq2$lfcSE,
  wald_stat_t = res_deseq2$stat,
  wald_pval = res_deseq2$pvalue,
  deseq_pval_adj = res_deseq2$padj,
  deseq2_baseMean = res_deseq2$baseMean,
  time = as.numeric(difftime(t1, t0, units = "secs"))
)


## Run Limma

t0 <- Sys.time()

limma_data <- edgeR::DGEList(pb_data$pseudobulk_data)
limma_data <- edgeR::calcNormFactors(limma_data)

# VOOM it up
v_dat <- limma::voom(limma_data, design_mat)

# Fit model
limma_fit <- limma::lmFit(v_dat, design_mat)

# Compute contrasts
limma_fit2 <- limma::contrasts.fit(limma_fit, contrast_mat[1, ])
limma_fit2 <- limma::eBayes(limma_fit2)

# Get test statistics
res_limma <- limma::topTable(limma_fit2, n = Inf)
# Subset to synthetic gene of interest
res_limma <- res_limma[which(rownames(res_limma) == gene_name), ]

t1 <- Sys.time()

# Create dataframe with results
beta_hat_limma <- coef(limma_fit)[gene_idx, ]
res_limma <- data.frame(
  datagen_model = data_model,
  trial = trial_idx,
  changed_covariate = covariate_name,
  baseline_contrast = covariate_contrast,
  true_contrast_value = beta_val,
  baseline_coef_value = beta_hat_power[beta_baseline_idx],
  changed_coef_value = beta_hat_power[beta_change_idx],
  estimated_baseline_coef = beta_hat_limma[beta_baseline_idx],
  estimated_changed_coef = beta_hat_limma[beta_change_idx],
  gene = gene_name,
  fit_model = "limma_voom",
  kernel_type = NA,
  contrast_indices = paste0(beta_baseline_idx, "_", beta_change_idx),
  contrast_val = res_limma$logFC,
  contrast_se = res_limma$logFC / res_limma$t,
  wald_stat_t = res_limma$t,
  wald_pval = res_limma$P.Value,
  limma_pval_adj = res_limma$adj.P.Val,
  limma_avgXPR = res_limma$AveExpr,
  limma_B = res_limma$B,
  time = as.numeric(difftime(t1, t0, units = "secs"))
)


## Save to file

out_df <- dplyr::bind_rows(res_deseq2, res_limma)


if (data_model == "spNNGP") {
  outfile_base <- paste0(
    "poisECM_power_natgen_kidney_",
    "datamodel_",
    data_model,
    "_kernel_",
    kernel_type,
    "_fitmodel_",
    "deseq2_limma_voom",
    "_trial_",
    trial_idx,
    "_gene_",
    gene_name,
    "_beta_",
    beta_val,
    "_covariate_",
    gsub(":", "_", gsub("/", "_", covariate_name))
  )
} else if (data_model %in% c("Leroux", "CAR", "SAR")) {
  outfile_base <- paste0(
    "poisECM_power_natgen_kidney_",
    "datamodel_",
    data_model,
    "_fitmodel_",
    "deseq2_limma_voom",
    "_trial_",
    trial_idx,
    "_gene_",
    gene_name,
    "_beta_",
    beta_val,
    "_covariate_",
    gsub(":", "_", gsub("/", "_", covariate_name))
  )
}
outfile_base <- paste0(out_path, outfile_base)

saveRDS(out_df,
        file = file.path(outfile_base, ".rds"),
        compress = TRUE)


## End

print("Done.")
