#' Aggregate files for timing results
#' @author Florica Constantine


# Libraries
library(dplyr)

# Paths
path <- "/scratch/users/spatialseq/natgen_kidney/"
processed_data_path <- file.path(path, "processed")
timing_data_path <- file.path(path, "timing")

# List of files to read in
file_list <- base::list.files(timing_data_path)

# Aggregate files
if (file.exists(paste0(processed_data_path, "timing_results.rds"))) {
  timing_out_all <- readRDS(paste0(processed_data_path, "timing_results.rds"))
} else {
  summary_out_list <- list()
  for (idx in 1:length(file_list)) {
    if (0 == (idx %% 250)) {
      cat(idx, "\n")
    }

    # Read in output
    out <- readRDS(paste0(timing_data_path, file_list[idx]))

    # Compute error in beta, just in case
    beta_SE <- sum((out$beta_hat - out$beta_true)^2)
    beta_norm2 <- sum((out$beta_true)^2)

    summary_out_list[[idx]] <- cbind(
      out$run_info,
      data.frame(
        data_model = out$data_model,
        gene_name = out$gene_name,
        beta_SE = beta_SE,
        beta_norm2 = beta_norm2,
        beta_rel_SE = beta_SE / beta_norm2,
        time_sec = out$alg_run_time
      )
    )
  }
  summary_out_all <- dplyr::bind_rows(summary_out_list)
  saveRDS(summary_out_all,
            paste0(processed_data_path, "timing_results.rds"))
}
