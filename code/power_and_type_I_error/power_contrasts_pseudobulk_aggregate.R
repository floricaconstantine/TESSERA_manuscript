#' Aggregate files for pseudobulk power results on both Poisson-Leroux and Poisson-spNNGP generated data
#' @author Florica Constantine

# Libraries
library(dplyr)

# Paths
path <- "/scratch/users/spatialseq/natgen_kidney/"
pb_data_path <- file.path(path, "power_contrasts")
processed_data_path <- file.path(path, "processed")

# List of files to read in
file_list <- base::list.files(pb_data_path)

# Aggregate files
out_list = list()
if (file.exists(file.path(processed_data_path, "power_contrasts_pb_results.rds"))) {
  out_all <- readRDS(file.path(processed_data_path, "power_contrasts_pb_results.rds"))
} else{
  summary_out_list <- list()
  for (idx in 1:length(file_list)) {
    if (0 == (idx %% 250)) {
      cat(idx, "\n")
    }
    
    # Read in output
    out_list[[idx]] <- readRDS(file.path(pb_data_path, file_list[idx]))
  }
  out_all <- dplyr::bind_rows(out_list)
  saveRDS(
    out_all,
    file.path(processed_data_path, "power_contrasts_pb_results.rds"),
    compress = "xz"
  )
}
