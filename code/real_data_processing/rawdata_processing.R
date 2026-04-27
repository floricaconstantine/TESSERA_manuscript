## Purpose of file: Process raw data (csvs of count data, meta data,
## and spatial coordinates) into a Seurat object

## Florica Constantine (florica AT berkeley DOT edu)

# Raw data sources
# Paper: https://www.nature.com/articles/s41588-024-01802-x
# GEO submission: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211785
# Spatial coordinate files: https://susztaklab.com/samui/files/

# Libraries
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(Matrix))
suppressMessages(library(sparseMatrixStats))

# Data paths
raw_data_path <- "../data/raw/"
processed_data_path <- "../data/processed/"

# Count matrix
count_mat <- readRDS(paste0(raw_data_path, "GSE211785_EXPORT_ST_counts.rds"))
# Meta data
meta_df <- read.csv(paste0(raw_data_path, "GSE211785_ST_metadata.txt"), sep =
                      "\t")

# Supplementary Table 1
supp_table1 <- read.csv(paste0(raw_data_path, "supp_table1.csv"), sep =
                          ",")
# Create column with sample names that match those in the meta data
supp_table1$orig.ident <- paste0(supp_table1$ID, "_ST")
# Expand meta data to include information from Supplementary table 1
meta_df <- dplyr::left_join(meta_df, supp_table1, by = "orig.ident")

# UMAP coordinates
umap_coords <- read.csv(paste0(raw_data_path, "GSE211785_EXPORT_ST_umap.txt"), sep =
                          "\t")

# Check ordering
stopifnot(identical(colnames(count_mat), meta_df$X))
stopifnot(identical(umap_coords$X, meta_df$X))

# Read in list of spatial coordinates and check what corresponds
coord_file_list <- list.files(paste0(raw_data_path, "visium_cell_locations/"), "*.csv")

# Samples in meta data
samples_in_metadata <- sort(unlist(lapply(unique(meta_df$orig.ident), function(fname) {
  strsplit(fname, "_")[[1]][1]
})))
# Samples with spatial coordinate data
samples_in_coords <- sort(unlist(lapply(coord_file_list, function(fname) {
  strsplit(fname, "_")[[1]][1]
})))
# Samples with spatial coordinate data but no meta data
cat("Samples with spatial coordinate data but no meta data:",
    samples_in_coords[which(!(samples_in_coords %in% samples_in_metadata))], "\n")
# Samples with meta data but no spatial coordinate data
cat("Samples with meta data but no spatial coordinate data:",
    samples_in_metadata[which(!(samples_in_metadata %in% samples_in_coords))], "\n")

# Number of cells per sample in meta data
meta_counts <- unlist(list(table(meta_df$orig.ident)))
# Number of cells per sample in spatial coordinate data
nrow_coords <- list()
for (fname in coord_file_list) {
  sample_name <- strsplit(fname, "_")[[1]][1]
  coords <- read.csv(paste0(raw_data_path, "visium_cell_locations/", fname))
  if ("in_tissue" %in% colnames(coords)) {
    nrow_coords[[paste0(sample_name, "_ST")]] <- sum(coords$in_tissue == 1)
  } else {
    nrow_coords[[paste0(sample_name, "_ST")]] <- nrow(coords)
  }
}
nrow_coords <- unlist(nrow_coords)
# Data frame showing number of cells per sample in spatial coordinate data v. meta data
sample_count_df <- data.frame(dplyr::bind_rows(nrow_coords, meta_counts))
rownames(sample_count_df) <- c("coord_files", "metadata")
print("Cell counts in spatial coordinate data v. metadata")
print(sample_count_df)


# Store only samples with coordinates
coord_list <- list()
count_list <- list()
meta_list <- list()
umap_list <- list()
for (fname in coord_file_list) {
  # Sample
  sample_name <- strsplit(fname, "_")[[1]][1]
  # Spatial coordinates for sample
  coords <- read.csv(paste0(raw_data_path, "visium_cell_locations/", fname))

  # Indices for cells in sample
  row_idx <- which(meta_df$orig.ident == paste0(sample_name, "_ST"))
  if (0 == length(row_idx)) {
    cat("SPATIAL FILENAME DOES NOT CORRESPOND TO METADATA:", fname, "\n")
    next
  }

  # Truncate cell IDs to remove everything after and including the underscore
  cell_id_subset <- unlist(lapply(strsplit(meta_df$X[row_idx], "_"), function(x) {
    x[1]
  }))
  # Match ordering of coordinates
  match_idx <- match(cell_id_subset, coords$barcode)
  stopifnot(identical(coords$barcode[match_idx], cell_id_subset))
  # Reorder
  coords <- coords[match_idx, ]
  # Recheck to be safe
  stopifnot(identical(coords$barcode, cell_id_subset))

  # Store
  coord_list[[sample_name]] <- coords
  count_list[[sample_name]] <- count_mat[, row_idx]
  meta_list[[sample_name]] <- meta_df[row_idx, ]
  umap_list[[sample_name]] <- umap_coords[row_idx, ]

  colnames(meta_list[[sample_name]])[1] <- "cell_id_metadata"
  colnames(umap_list[[sample_name]])[1] <- "cell_id_umap"
  colnames(coord_list[[sample_name]])[1] <- "cell_id_coords"

  # Make sure all lengths are identical
  stopifnot(1 == length(unique(c(
    length(row_idx),
    nrow(coords),
    length(unique(cell_id_subset)),
    length(cell_id_subset)
  ))))
}

# Check that the same genes are present in all samples in the same order
for (idx in 2:length(count_list)) {
  stopifnot(identical(rownames(count_list[[1]]), rownames(count_list[[idx]])))
}

# Recreate metadata by binding across samples
# and also incorporating the spatial and UMAP coordinates
meta_subset <- dplyr::bind_cols(
  dplyr::bind_rows(meta_list),
  dplyr::bind_rows(coord_list),
  dplyr::bind_rows(umap_list)
)
# Check order
stopifnot(identical(meta_subset$cell_id_metadata, meta_subset$cell_id_umap))
# Bind count matrices across samples
count_mat_subset <- do.call(cbind, count_list)
# Remove genes with all zero counts
cat("Count data dimensions:", dim(count_mat_subset), "\n")
zero_genes <- which(Matrix::rowSums(count_mat_subset) == 0)
if (0 < length(zero_genes)) {
  count_mat_subset <- count_mat_subset[-zero_genes, ]
  cat("Number of genes with all-zero counts (dropped):", length(zero_genes), "\n")
  cat("Count data dimensions:", dim(count_mat_subset), "\n")
}

# Order by largest to smallest variance in gene counts
ct_var <- sparseMatrixStats::rowVars(count_mat_subset)
count_mat_subset <- count_mat_subset[order(ct_var, decreasing = TRUE), ]

# Check order
stopifnot(identical(colnames(count_mat_subset), meta_subset$cell_id_metadata))

# See how many genes are identically zero for each sample
# Including ones without coordinates
for (samp in unique(meta_df$orig.ident)) {
  row_idx <- which(meta_df$orig.ident == samp)
  cat(samp,
      "; Has metadata: ",
      as.character(samp %in% unique(meta_subset$orig.ident)),
      "; Number of all-zero Genes: ",
      sum(Matrix::rowSums(count_mat[, row_idx]) == 0),
      "; Number of non-zero Genes: ",
      sum(Matrix::rowSums(count_mat[, row_idx]) > 0),
      "; Number of Cells: ",
      length(row_idx),
      "\n", sep=""
    )
}

# Save kidney count data and meta data as a Seurat object
seu_kid <- Seurat::CreateSeuratObject(count_mat_subset, meta.data = meta_subset)
saveRDS(seu_kid, paste0(processed_data_path, "seu_kidney.rds"))
