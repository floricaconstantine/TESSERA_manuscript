#' Aggregate files for single sample results that used ALL cells
#' @author Florica Constantine

# Libraries
library(dplyr)

# Paths
path <- "/scratch/users/spatialseq/natgen_kidney/"
singlesample_data_path <- file.path(path, "singlesample_allcells_output_onlyinteract")
processed_data_path <- file.path(path, "processed")

# List of files to read in
file_list <- base::list.files(singlesample_data_path)

# Aggregate files
if(file.exists(file.path(processed_data_path, "single_sample_all_cells_perfsummary_onlyinteract.rds"))){
  summary_out_all <- readRDS(file.path(processed_data_path, "single_sample_all_cells_perfsummary_onlyinteract.rds"))
}else{
  summary_out_list <- list()
  for(idx in 1:length(file_list)){
    if (0 == (idx %% 250)) {
      cat(idx, "\n")
    }

    # Read in output
    out <- readRDS(file.path(singlesample_data_path, file_list[idx]))

    # Add in information that wasn't stored
    trial_idx <- out$sim_settings$trial
    out$summary_df <- cbind(data.frame(trial=trial_idx), out$summary_df)

    summary_out_list[[idx]] <- out$summary_df
  }
  summary_out_all <- dplyr::bind_rows(summary_out_list)
  saveRDS(summary_out_all, file.path(processed_data_path, "single_sample_all_cells_perfsummary_onlyinteract.rds"))
}
