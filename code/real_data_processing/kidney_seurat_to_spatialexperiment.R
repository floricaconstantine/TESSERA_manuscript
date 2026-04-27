# This script converts the Seurat object (storing the Kidney data)
# to a SpatialExperiment object

# Libraries for Bioconductor data objects
library(S4Vectors)
library(MatrixGenerics)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)

# Seurat
library(Seurat)

# Load object
seu <- readRDS(file.path(".", "..", "data", "processed", "seu_kidney.rds"))

# Which covariates to use
celltype_col <- "celltype"
sample_col <- "orig.ident"
group_col <- "Group"
coord_cols <- c("pxl_row_in_fullres", "pxl_col_in_fullres")

# Convert to SCE
sce <- Seurat::as.SingleCellExperiment(seu)

# Extract spatial coordinates
spatialCoords <- as.matrix(seu@meta.data[, coord_cols])

# Convert to SpatialExperiment
spe <- SpatialExperiment::SpatialExperiment(
  assays = SummarizedExperiment::assays(sce),
  rowData = SummarizedExperiment::rowData(sce),
  colData = SummarizedExperiment::colData(sce),
  metadata = S4Vectors::metadata(sce),
  reducedDims = SingleCellExperiment::reducedDims(sce),
  altExps = SingleCellExperiment::altExps(sce),
  sample_id = sample_col,
  spatialCoords = spatialCoords
)

# Make sure names match
rownames(SpatialExperiment::spatialCoords(spe)) <- rownames(SummarizedExperiment::colData(spe))
colnames(SingleCellExperiment::counts(spe)) <- rownames(SummarizedExperiment::colData(spe))

# Calculate variance across all spots for each gene
rv <- MatrixGenerics::rowVars(SingleCellExperiment::counts(spe))
names(rv) <- rownames(spe)
# Order the entire object by raw count variance (highest to lowest)
spe <- spe[order(rv, decreasing = TRUE), ]

# Save SpatialExperiment object
# Use "xz" compression for smallest file size
saveRDS(spe, file.path(".", "kidney_data_spe.rds"), compress = "xz")
