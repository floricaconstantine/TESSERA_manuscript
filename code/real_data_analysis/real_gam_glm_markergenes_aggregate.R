#' Aggregate GAM, NB GLM, and Poisson GLM results from between-cell-type comparisons in a real-world kidney dataset
#' @author Florica Constantine

# Libraries
library(dplyr)

# Paths
path <-  "/scratch/users/spatialseq/natgen_kidney/"
processed_data_path <- file.path(path, "processed")
glm_gam_modelfit_data_path <- file.path(path, "glm_offset_real_output_onlyinteract")

# Subset to top 3000 genes by variance
seu_kid <- readRDS(file.path(processed_data_path, "seu_kidney.rds"))
dim(seu_kid)
seu_kid <- seu_kid[1:3000, ]
gene_keep <- rownames(seu_kid)
length(gene_keep)


# Contrast matrix
prep_obj <- readRDS(file.path(processed_data_path, "prepData_lattice_onlyinteract.rds"))
X_mat <- Reduce(rbind, prep_obj$X_list)
cell_types_contrasts <- unique(seu_kid@meta.data$celltype)
contrast_combs <- combn(cell_types_contrasts, 2)
contrast_mat <- matrix(data=0, nrow=ncol(contrast_combs), ncol=ncol(X_mat))
rownames(contrast_mat) <- seq(1, ncol(contrast_combs))
for (idx in 1:ncol(contrast_combs)) {
  col_idx <- which(grepl(contrast_combs[1, idx], colnames(X_mat)))
  contrast_mat[idx, col_idx] <- 1 / length(col_idx)
  col_idx <- which(grepl(contrast_combs[2, idx], colnames(X_mat)))
  contrast_mat[idx, col_idx] <- -1 / length(col_idx)
  rownames(contrast_mat)[idx] <- paste(contrast_combs[, idx], collapse = " v. ")
}
colnames(contrast_mat) <- colnames(X_mat)


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
    "glm_gam_real_results_waldstats_marker_onlyinteract.rds"
  )
)) {
  wald_df_marker <- readRDS(
    file.path(
      processed_data_path,
      "glm_gam_real_results_waldstats_marker_onlyinteract.rds"
    )
  )
} else{
  wald_out_list <- list()
  for (idx in 1:length(file_list)) {
    if (0 == (idx %% 250)) {
      cat(idx, "\n")
    }

    out <- readRDS(file.path(glm_gam_modelfit_data_path, file_list[idx]))

    for (c_idx in 1:nrow(contrast_mat)) {
      # Contrast
      R_beta <- as.numeric(t(contrast_mat[c_idx, ]) %*% out$beta_hat[1:ncol(contrast_mat)])
      # SE^2
      RVR <- as.numeric(t(contrast_mat[c_idx, ]) %*% out$beta_vcov[1:ncol(contrast_mat), 1:ncol(contrast_mat)] %*% (contrast_mat[c_idx, ]))

      wald_out_list[[1 + length(wald_out_list)]] <- data.frame(
        gene = as.character(out$wald_df$gene[1]),
        fit_model = as.character(out$wald_df$fit_model[1]),
        contrast_name = rownames(contrast_mat)[c_idx],
        contrast_indices = paste0(which(contrast_mat[c_idx, ] != 0), collapse =
                                    "_"),
        contrast_val = R_beta,
        contrast_se = sqrt(RVR),
        wald_stat_t = R_beta / sqrt(RVR)
      )
    }
  }
  wald_df_marker <- dplyr::bind_rows(wald_out_list)
  saveRDS(
    wald_df_marker,
    file.path(
      processed_data_path,
      "glm_gam_real_results_waldstats_markergene_onlyinteract.rds"
    ),
    compress = "xz"
  )
}
