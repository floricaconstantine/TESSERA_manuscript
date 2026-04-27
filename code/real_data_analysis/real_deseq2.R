#' Run DESeq2 on real-world kidney data and calculate test statistics for within-cell type between-condition comparisons
#' @author Florica Constantine

# Scripts
source("aggregate_data.R")

## Libraries
library(dplyr)
library(reshape2)
library(Rfast)
library(Matrix)
library(spatstat.geom)
library(ggplot2)
library(DESeq2)

## Global parameters
GLOBAL_SEED <- 2025
set.seed(GLOBAL_SEED)

## Paths
path <-  "/scratch/users/spatialseq/natgen_kidney/"
processed_data_path <- file.path(path, "processed")

## Read in kidney data
seu_kid <- readRDS(file.path(processed_data_path, "seu_kidney.rds")) 

# Subset to top 3000 genes by variance
dim(seu_kid)
seu_kid <- seu_kid[1:3000, ]
dim(seu_kid)

# Call function to create pseudobulk and bulk data
sample_col <- "orig.ident"
celltype_col <- "celltype"
cell_id_col <- "cell_id_metadata"
group_col <- "Group"
pb_data <- make_pseudobulk_bulk(
  seu_kid@assays$RNA@layers$counts,
  seu_kid@meta.data,
  individual_colname = sample_col,
  cell_type_colname = celltype_col,
  cell_id_colname = cell_id_col
)
pseudobulk_meta <- pb_data$pseudobulk_meta
pseudobulk_meta <- pseudobulk_meta[, c("pb_id", sample_col, celltype_col, group_col)]
pseudobulk_meta <- dplyr::distinct(pseudobulk_meta)
rownames(pseudobulk_meta) = pseudobulk_meta$pb_id
pseudobulk_meta <- pseudobulk_meta[,-which(colnames(pseudobulk_meta) == "pb_id")]

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
design_formula <- as.formula(paste0(
  "~ 0 + ",
  paste(
    group_col,
    c(celltype_col, new_nested_sample_col),
    sep = ":",
    collapse = " + "
  )
))
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

# Unique cell types and groups
cell_types <- unique(pseudobulk_meta[, celltype_col])
group_types <- c("Control", "DKD", "HKD")
# Get all pairs of groups
group_comparisons <- combn(group_types, 2)
# Loop over cell types and group pairs
contrasts_list <- list()
contrasts_names <- c()
for(ct_idx in 1:length(cell_types)){
  for(gc_idx in 1:ncol(group_comparisons)){
    tmp <- rep(0, ncol(design_mat))
    # Group 1 - Group 2 for one cell type
    tmp[which(colnames(design_mat) == paste0(group_col, group_comparisons[1, gc_idx], ":", celltype_col, cell_types[ct_idx]))] <- 1
    tmp[which(colnames(design_mat) == paste0(group_col, group_comparisons[2, gc_idx], ":", celltype_col, cell_types[ct_idx]))] <- -1

    contrasts_names <- c(contrasts_names, paste0(cell_types[ct_idx], ":", group_comparisons[1, gc_idx], "-", group_comparisons[2, gc_idx]))
    contrasts_list[[1 + length(contrasts_list)]] <- tmp
  }
}
contrast_mat <- Reduce(rbind, contrasts_list)
colnames(contrast_mat) <- colnames(design_mat)
rownames(contrast_mat) <- contrasts_names
contrast_mat <- contrast_mat[-which(rowSums(abs(contrast_mat)) != 2), ]

# Put data into a DESeq2 object
dds <- DESeq2::DESeqDataSetFromMatrix(countData = pb_data$pseudobulk_data,
                                     colData = pseudobulk_meta,
                                     design = design_mat)
# Run DESeq2
dds <- DESeq2::DESeq(dds)

# Extract DESeq2 results
# Loop over contrasts
res_list <- list()
for (c_idx in 1:nrow(contrast_mat)) {
  cat(c_idx, rownames(contrast_mat)[c_idx], Sys.time(), "\n")
  res <- DESeq2::results(dds, contrast = contrast_mat[c_idx, ], parallel = TRUE)
  res <- as.data.frame(res)

  # Add in gene and model information, as well as contrast information
  res_meta <- data.frame(gene=rownames(seu_kid),
                         fit_model="DESeq2",
                         contrast_description=rownames(contrast_mat)[c_idx]
                         )

  res <- cbind(res_meta, res)
  res_list[[c_idx]] <- res
}
res <- dplyr::bind_rows(res_list)
# Globally adjust p-values (not used)
res$pval_adj <- p.adjust(res$pvalue, method = "BH")

# Inspect number of significant genes per contrast
data.frame(res
           %>% dplyr::group_by(fit_model, contrast_description)
           %>% dplyr::summarise(n_sig=sum(pval_adj < 0.05, na.rm=TRUE))
)

## Save to file
saveRDS(res,
        file = file.path(processed_data_path, "DESeq2_real_natgen_kidney.rds"),
        compress = TRUE)

## End
print("Done.")
