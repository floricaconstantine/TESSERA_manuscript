#' Run TESSERA on real-world kidney data
#' @author Florica Constantine

## Libraries
library(TESSERA)
library(dplyr)
library(reshape2)
library(Rfast)
library(Matrix)
library(spatstat.geom)
library(ggplot2)

## Global parameters
GLOBAL_SEED <- 2025
set.seed(GLOBAL_SEED)

## Paths
path <-  "/scratch/users/spatialseq/natgen_kidney/"
processed_data_path <- file.path(path, "processed")

## Model types
model_list <- c("Leroux", "spNNGP", "CAR", "SAR")

## Read in kidney data
seu_kid <- readRDS(file.path(processed_data_path, "seu_kidney.rds"))
# Subset to top 3000 genes by variance
dim(seu_kid)
seu_kid <- seu_kid[1:3000, ]
dim(seu_kid)

# Cluster arguments
job_df <- expand.grid(rownames(seu_kid[["RNA"]]$counts), model_list)
colnames(job_df) <- c("gene", "model")

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
  if ((1 > run_idx) | (nrow(job_df) < run_idx)) {
    stop("run_idx must be between 1 and # of possible jobs")
  }
}

# Gene
gene_name <- job_df[run_idx, "gene"]
print(gene_name)
# Model type
model_type <- job_df[run_idx, "model"]
# If an spNNGP model, type of kernel
kernel_type <- "Mat"

## Output from prep_data
if (model_type == "spNNGP") {
  prep_out <- readRDS(file.path(processed_data_path, "prepData_spNNGP_onlyinteract.rds"))
} else {
  prep_out <- readRDS(file.path(processed_data_path, "prepData_lattice_onlyinteract.rds"))
}

### Run TESSERA algorithm ###

if (model_type == "spNNGP") {
  alg_out <- TESSERA::TESSERA_spNNGP(
    poisECMData_obj = prep_out,
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
    poisECMData_obj = prep_out,
    gene_name = gene_name,
    model_type = model_type,
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
wald_df <- TESSERA::calc_Wald_statistics(alg_out, prep_out$contrast_mat)
# Add to alg_out
alg_out$wald_df <- wald_df

## Save to file
if (model_type == "spNNGP") {
  outfile_base <- paste0(
    "poisECM_real_natgen_kidney_",
    "model_",
    model_type,
    "_kernel_",
    kernel_type,
    "_gene_",
    gene_name
  )
} else{
  outfile_base <- paste0("poisECM_real_natgen_kidney_",
                         "model_",
                         model_type,
                         "_gene_",
                         gene_name)
}
outfile_base <- file.path(processed_data_path, outfile_base)
saveRDS(alg_out,
        file = paste0(outfile_base, ".rds"),
        compress = TRUE)

## End
print("Done.")
