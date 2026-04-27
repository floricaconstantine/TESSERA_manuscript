# Prepare data to run TESSERA algorithm

library(TESSERA)
library(Matrix)
library(ggplot2)


# Data path
processed_data_path <- "../data/processed/"
# Read in kidney data
seu_kid <- readRDS(paste0(processed_data_path, "seu_kidney.rds"))

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

# Create dataframe
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


### Find distance threshold for cell-cell adjacency matrix ###

k_search <- 20
dist_out <- TESSERA::visualizeNeighborDistances(
  meta_data = seu_kid@meta.data,
  sample_col = sample_col,
  coord_data = seu_kid@meta.data[, coord_cols],
  k_search = k_search
)

# Make a plot to visualize (look for gaps)
# After 6 neighbors we see a large gap
# This is in line with the grid for the Visium technology (cells typically have 6 neighbors)
p <- dist_out$distances_plot
p

# Distance threshold to decide which cells are adjacent neighbors
# Compute mean distance across all samples for the 6th neighbor and 7th neighbor
# Then compute the mean distance BETWEEN the 6th and 7th neighbor
D_THRESH <- mean(colMeans(dist_out$mean_nb_dist)[6:7]) # Distance threshold


### Set up contrasts ###

# Unique cell types and groups
cell_types <- levels(seu_kid@meta.data[, celltype_col])
group_types <- levels(seu_kid@meta.data[, group_col])

# Get all pairs of groups
group_comparisons <- combn(group_types, 2)

# Loop over cell types and group pairs
contrasts_list <- list()
contrasts_names <- c()
for(ct_idx in 1:length(cell_types)){
  for(gc_idx in 1:ncol(group_comparisons)){
    tmp <- rep(0, ncol(X_mat))
    # Group 1 - Group 2 for one cell type
    tmp[which(colnames(X_mat) == paste0(group_col, group_comparisons[1, gc_idx], ":", celltype_col, cell_types[ct_idx]))] <- 1
    tmp[which(colnames(X_mat) == paste0(group_col, group_comparisons[2, gc_idx], ":", celltype_col, cell_types[ct_idx]))] <- -1

    contrasts_names <- c(contrasts_names, paste0(cell_types[ct_idx], ":", group_comparisons[1, gc_idx], "-", group_comparisons[2, gc_idx]))
    contrasts_list[[1 + length(contrasts_list)]] <- tmp
  }
}
contrast_mat <- Reduce(rbind, contrasts_list)
colnames(contrast_mat) <- colnames(X_mat)
rownames(contrast_mat) <- contrasts_names

# Which rows and columns to drop
# Drop columns in X corresponding to cell type, group combinations with no cells
# Also drop samples that do not exist
# Drop rows in contrast matrix that involve dropped columns
drop_cols <- which(colSums(X_mat) == 0)
drop_rows <- which(rowSums(abs(contrast_mat[, drop_cols])) > 0)
contrast_mat <- contrast_mat[-drop_rows, ]
contrast_mat <- contrast_mat[, -drop_cols]
X_mat <- X_mat[, -drop_cols]

stopifnot(ncol(X_mat) == ncol(contrast_mat))
stopifnot(ncol(X_mat) == Matrix::rankMatrix(X_mat)[1])
print(dim(contrast_mat))
print(dim(X_mat))


### Prepare data for TESSERA algorithm ###

prep_out_lattice <- TESSERA::prep_data(
  x = seu_kid[["RNA"]]$counts,
  meta_data = seu_kid@meta.data,
  sample_col = sample_col,
  design_mat = X_mat,
  coord_data = seu_kid@meta.data[, coord_cols],
  D_THRESH = D_THRESH,
  k_search = k_search,
  model_type = "ALL"
)
prep_out_lattice$contrast_mat <- contrast_mat
saveRDS(prep_out_lattice,
        paste0(processed_data_path, "prepData_lattice_onlyinteract.rds"))

prep_out_spnngp <- TESSERA::prep_data(
  x = seu_kid[["RNA"]]$counts,
  meta_data = seu_kid@meta.data,
  sample_col = sample_col,
  design_mat = X_mat,
  coord_data = seu_kid@meta.data[, coord_cols],
  D_THRESH = NULL,
  k_search = k_search,
  model_type = "ALL"
)
prep_out_spnngp$contrast_mat <- contrast_mat
saveRDS(prep_out_spnngp,
        paste0(processed_data_path, "prepData_spNNGP_onlyinteract.rds"))
