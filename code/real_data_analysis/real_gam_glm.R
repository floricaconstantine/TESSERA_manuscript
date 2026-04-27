#' Run GAM, NB-GLM, and Poisson-GLM methods on real-world kidney data
#' @author Florica Constantine


## Libraries

library(dplyr)
library(Matrix)
library(poisECM)


## Global parameters

GLOBAL_SEED <- 2025
set.seed(GLOBAL_SEED)


## Paths

path <-  "/scratch/users/spatialseq/natgen_kidney/"
# Input data path
processed_data_path <- file.path(path, "processed")
# Path to write to
out_path <- file.path(path, "glm_offset_real_output_onlyinteract")

## Scripts

source("TESSERA_comparison_methods.R")

## Output from prepData

prep_out <- readRDS(file.path(processed_data_path, "prepData_lattice_onlyinteract.rds"))


## Run parameters

N_TOP_GENES <- 3000 # Only look at top 3000 genes
model_list <- c("Poisson_GLM", "NB_GLM", "GAM")
gene_names <- rownames(prep_out$counts_list[[1]])[1:N_TOP_GENES]

# Cluster arguments
job_df <- expand.grid(gene_names, model_list)
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
# Model type
model_type <- job_df[run_idx, "model"]


## Run model

if ("Poisson_GLM" == model_type) {
  alg_out <- glm_wrapper(
    lapply(prep_out$counts_list, function (x) {
      x[gene_name, ]
    }),
    prep_out$X_list,
    library_size_list = prep_out$library_size_list
  )
} else if ("NB_GLM" == model_type) {
  alg_out <- glm_nb_wrapper(
    lapply(prep_out$counts_list, function (x) {
      x[gene_name, ]
    }),
    prep_out$X_list,
    library_size_list = prep_out$library_size_list
  )
} else if ("GAM" == model_type) {
  alg_out <- mgcv_gam_wrapper(
    lapply(prep_out$counts_list, function (x) {
      x[gene_name, ]
    }),
    prep_out$X_list,
    prep_out$coords_list,
    library_size_list = prep_out$library_size_list
  )
}


## Wald statistics

contrast_mat <- prep_out$contrast_mat
if (model_type %in% c("Poisson_GLM", "NB_GLM", "GAM")) {
  wald_df <- list()
  for (c_idx in 1:nrow(contrast_mat)) {
    # Contrast
    R_beta <- as.numeric(t(contrast_mat[c_idx, ]) %*% alg_out$coefficients[1:ncol(contrast_mat)])
    # SE^2
    RVR <- as.numeric(t(contrast_mat[c_idx, ]) %*% vcov(alg_out)[1:ncol(contrast_mat), 1:ncol(contrast_mat)] %*% (contrast_mat[c_idx, ]))

    wald_df[[c_idx]] <- data.frame(
      gene = as.character(gene_name),
      fit_model = as.character(model_type),
      contrast_name = rownames(contrast_mat)[c_idx],
      contrast_indices = paste0(which(contrast_mat[c_idx, ] != 0), collapse =
                                  "_"),
      contrast_val = R_beta,
      contrast_se = sqrt(RVR),
      wald_stat_t = R_beta / sqrt(RVR)
    )
  }

  # Add to alg_out
  alg_out$wald_df <- dplyr::bind_rows(wald_df)
  rownames(alg_out$wald_df) <- rownames(contrast_mat)

  # Compute p-values via Wilks' Theorem
  alg_out$wald_df$wald_pval <- pchisq(alg_out$wald_df$wald_stat_t^2, 1, lower.tail = FALSE)

  alg_out$beta_hat <- alg_out$coefficients[1:ncol(contrast_mat)]
  alg_out$beta_vcov <- vcov(alg_out)[1:ncol(contrast_mat), 1:ncol(contrast_mat)]
}


## Get statistics for GLM/GAM

if (model_type %in% c("Poisson_GLM", "NB_GLM", "GAM")) {
  MSE <- rep(NA, length(prep_out$counts_list))
  counts2 <- rep(NA, length(MSE))
  data_ll <- rep(NA, length(MSE))
  Moran_pred <- rep(NA, length(MSE))
  Moran_counts <- rep(NA, length(MSE))
  Moran_resid <- rep(NA, length(MSE))
  sample_indexing <- c(0, cumsum(sapply(prep_out$library_size_list, length)))

  alg_pred <- predict(alg_out, type = "response")
  lib_size_list <- Reduce(c, prep_out$library_size_list)
  # alg_pred <- alg_pred * lib_size_list  # Uncomment if not using offset
  alg_resid <- Reduce(c, lapply(prep_out$counts_list, function (x) {
    x[gene_name, ]
  })) - alg_pred
  for (s_idx in 1:length(MSE)) {
    fit_local <- alg_pred[(1 + sample_indexing[s_idx]):sample_indexing[s_idx + 1]]
    resid_local <- alg_resid[(1 + sample_indexing[s_idx]):sample_indexing[s_idx + 1]]
    MSE[s_idx] <- mean(resid_local^2, na.rm = TRUE)
    counts2[s_idx] <- mean(prep_out$counts_list[[s_idx]][gene_name, ]^2)

    if ("Poisson_GLM" == model_type || "GAM" == model_type) {
      data_ll[s_idx] <- sum(dpois(prep_out$counts_list[[s_idx]][gene_name, ], fit_local, log =
                                    TRUE))
    } else if ("NB_GLM" == model_type) {
      data_ll[s_idx] <- sum(dnbinom(
        prep_out$counts_list[[s_idx]][gene_name, ],
        size = alg_out$theta,
        mu = fit_local,
        log = TRUE
      ))
    }


    Moran_resid[s_idx] <- TESSERA::calc_moran(resid_local,
                                              prep_out$coords_list[[s_idx]][, 1],
                                              prep_out$coords_list[[s_idx]][, 2])[1]
    Moran_pred[s_idx] <- TESSERA::calc_moran(fit_local,
                                             prep_out$coords_list[[s_idx]][, 1],
                                             prep_out$coords_list[[s_idx]][, 2])[1]
    Moran_counts[s_idx] <- TESSERA::calc_moran(prep_out$counts_list[[s_idx]][gene_name, ],
                                               prep_out$coords_list[[s_idx]][, 1],
                                               prep_out$coords_list[[s_idx]][, 2])[1]
  }
  MSE_total <- mean(alg_out$residuals^2)
  counts2_total <- mean(Reduce(c, lapply(prep_out$counts_list, function (x) {
    x[gene_name, ]
  }))^2)

  if ("NB_GLM" == model_type) {
    theta_local <- alg_out$theta
  } else {
    theta_local <- NA
  }
  # Add to object
  summary_out_all <- data.frame(
    gene = as.character(gene_name),
    fit_model = as.character(model_type),
    sample = names(prep_out$coords_list),
    n_cells = sapply(prep_out$coords_list, nrow),
    MSE_counts_sample = MSE,
    MSE_counts_total = MSE_total,
    Mean_counts2_sample = counts2,
    Mean_counts2_total = counts2_total,
    Moran_counts = Moran_counts,
    Moran_predictions = Moran_pred,
    Moran_residuals = Moran_resid,
    data_ll = data_ll,
    time = as.numeric(alg_out$time, units = "secs"),
    theta = theta_local,
    Group = sapply(prep_out$covariates_list, function(x) {
      x$Group[1]
    })
  )
  alg_out$performanceSummary <- summary_out_all
}


## Save to file

output_list <- list(
  performanceSummary = alg_out$performanceSummary,
  wald_df = alg_out$wald_df,
  beta_hat = alg_out$beta_hat,
  coefficients = alg_out$coefficients,
  beta_vcov = alg_out$beta_vcov,
  beta_vcov_full = vcov(alg_out),
  gene_name = as.character(gene_name),
  model_type = model_type
)

outfile_base <- paste0("real_natgen_kidney",
                       "_model_",
                       model_type,
                       "_gene_",
                       gene_name)
outfile_base <- file.path(out_path, outfile_base)
saveRDS(output_list,
        file = paste0(outfile_base, ".rds"),
        compress = TRUE)


## End

print("Done.")
