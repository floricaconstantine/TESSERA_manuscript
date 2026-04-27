#' Aggregate TESSERA results for between-cell-type comparisons in a real-world kidney dataset
#' @author Florica Constantine

# Libraries
library(dplyr)
library(TESSERA)

# Paths
path <- "/scratch/users/spatialseq/natgen_kidney/"
modelfit_data_path <- file.path(path, "real_output_onlyinteract")
processed_data_path <- file.path(path, "processed")

# Which covariates to use
celltype_col <- "celltype"
sample_col <- "orig.ident"
group_col <- "Group"
covariates_keep <- c(group_col, celltype_col) # Use when creating design formula
coord_cols <- c("pxl_row_in_fullres", "pxl_col_in_fullres")

### Create design matrix ###

# Reference celltype, sample, and condition
ref_celltype <- "C_TAL"
ref_sample <- "HK2874_ST"
ref_group <- "Control"

# Subset to top 3000 genes by variance
seu_kid <- readRDS(file.path(processed_data_path, "seu_kidney.rds"))
dim(seu_kid)
seu_kid <- seu_kid[1:3000, ]
gene_keep <- rownames(seu_kid)
length(gene_keep)

# Re-level factors
seu_kid@meta.data[, sample_col] <- relevel(factor(seu_kid@meta.data[, sample_col]), ref = ref_sample)
seu_kid@meta.data[, celltype_col] <- relevel(factor(seu_kid@meta.data[, celltype_col]), ref = ref_celltype)
seu_kid@meta.data[, group_col] <- relevel(factor(seu_kid@meta.data[, group_col]), ref = ref_group)

## If we want a sample/patient effect; we need to nest this inside disease
# Create a new 'nested' variable as in the DESeq2 vignette
control_ids <- unique(seu_kid@meta.data[, sample_col][seu_kid@meta.data[, group_col] == "Control"])
DKD_ids <- unique(seu_kid@meta.data[, sample_col][seu_kid@meta.data[, group_col] == "DKD"])
HKD_ids <- unique(seu_kid@meta.data[, sample_col][seu_kid@meta.data[, group_col] == "HKD"])
new_id <- rep(0, nrow(seu_kid@meta.data))
for (idx in 1:length(control_ids)) {
  new_id[seu_kid@meta.data[, sample_col] == control_ids[idx]] <- idx
}
for (idx in 1:length(DKD_ids)) {
  new_id[seu_kid@meta.data[, sample_col] == DKD_ids[idx]] <- idx
}
for (idx in 1:length(HKD_ids)) {
  new_id[seu_kid@meta.data[, sample_col] == HKD_ids[idx]] <- idx
}
seu_kid@meta.data$nested_id <- factor(new_id)
new_nested_sample_col <- "nested_id"

data_mat <- seu_kid@meta.data
data_mat$z <- seu_kid@assays$RNA@layers$counts[1, ]

# Design formula
design_formula <- as.formula(paste0(
  "z ~ ",
  "0 + ",
  paste(
    covariates_keep[1],
    c(covariates_keep[-c(1)], new_nested_sample_col),
    sep = ":",
    collapse = " + "
  )
))
print(design_formula)
# Create design matrix
X_mat <- model.matrix(design_formula, data_mat)

# Rename sample/patient effects for clarity later
for (idx in 2:length(DKD_ids)) {
  colnames(X_mat) <- gsub(
    paste0("GroupDKD:nested_id", idx),
    paste0("GroupDKD:", DKD_ids[idx]),
    colnames(X_mat)
  )
}
for (idx in 2:length(HKD_ids)) {
  colnames(X_mat) <- gsub(
    paste0("GroupHKD:nested_id", idx),
    paste0("GroupHKD:", HKD_ids[idx]),
    colnames(X_mat)
  )
}
for (idx in 2:length(control_ids)) {
  colnames(X_mat) <- gsub(
    paste0("GroupControl:nested_id", idx),
    paste0("GroupControl:", control_ids[idx]),
    colnames(X_mat)
  )
}


cell_types_contrasts <- unique(seu_kid@meta.data[, celltype_col]) # ALL cell type pairs
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

# Drop columns without entries
drop_cols <- which(colSums(X_mat) == 0)
X_mat <- X_mat[, -drop_cols]
contrast_mat <- contrast_mat[, -drop_cols]


### Process TESSERA data ###

# List of files to read in
pattern_type <- "Leroux"
file_list <- base::list.files(modelfit_data_path, pattern = pattern_type)

# List of files to read in
pattern_type <- "Leroux"
full_file_list <- base::list.files(modelfit_data_path, pattern = pattern_type)
file_list <- c()
for(f_idx in full_file_list){
  gene_name <- sub(".*gene_(.*)\\.rds", "\\1", f_idx)
  if(gene_name %in% gene_keep){
    file_list <- c(file_list, f_idx)
  }
}
length(file_list)

# Get Wald stats
wald_celltype_list <- list()
for(idx in 1:length(file_list)){
  if (0 == (idx %% 250)) {
    cat(idx, "\n")
  }
  
  out <- readRDS(file.path(modelfit_data_path, file_list[idx]))
  wald_out <- TESSERA::calc_Wald_statistics(out, contrast_mat)
  wald_out$contrast_description <- rownames(wald_out)
  wald_celltype_list[[idx]] <- wald_out
}
wald_df_celltype <- dplyr::bind_rows(wald_celltype_list)


### Save wald_df_celltype
saveRDS(wald_df_celltype, file.path(processed_data_path, "real_results_waldstats_Leroux_between_celltypes_onlyinteract.rds"))

