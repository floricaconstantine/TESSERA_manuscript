#' Aggregate TESSERA and GAM/GLM results from within-cell-type between-condition comparisons in a real-world kidney dataset
#' @author Florica Constantine

# Libraries
library(TESSERA)
library(dplyr)

# Paths
path <- "/scratch/users/spatialseq/natgen_kidney/"
modelfit_data_path <- file.path(path, "real_output_onlyinteract")
glm_gam_modelfit_data_path <- file.path(path, "glm_offset_real_output_onlyinteract")
processed_data_path <- file.path(path, "processed")

# Output of prep_data
prep_obj <- readRDS(file.path(processed_data_path, "prepData_lattice_onlyinteract.rds"))

# Subset to top 3000 genes by variance
seu_kid <- readRDS(file.path(processed_data_path, "seu_kidney.rds"))
dim(seu_kid)
seu_kid <- seu_kid[1:3000, ]
gene_keep <- rownames(seu_kid)
length(gene_keep)

### TESSERA ###
# List of files to read in
pattern_types <- c("Leroux", "spNNGP")
for(p_name in pattern_types){
  pattern_type <- p_name
  full_file_list <- base::list.files(modelfit_data_path, pattern = pattern_type)
  file_list <- c()
  for(f_idx in full_file_list){
    gene_name <- sub(".*gene_(.*)\\.rds", "\\1", f_idx)
    if(gene_name %in% gene_keep){
      file_list <- c(file_list, f_idx)
    }
  }
  length(file_list)
  
  if(file.exists(file.path(processed_data_path, paste0("real_results_perfsummary_", pattern_type, "_onlyinteract.rds")))
     & file.exists(file.path(processed_data_path, paste0("real_results_waldstats_", pattern_type, "_onlyinteract.rds")))){
    summary_out_all <- readRDS(file.path(processed_data_path, paste0("real_results_perfsummary_", pattern_type, "_onlyinteract.rds")))
    wald_df_all <- readRDS(file.path(processed_data_path, paste0("real_results_waldstats_", pattern_type, "_onlyinteract.rds")))
  }else{
    summary_out_list <- list()
    wald_out_list <- list()
    for(idx in 1:length(file_list)){
      if (0 == (idx %% 250)) {
        cat(idx, "\n")
      }
  
      out <- readRDS(paste0(modelfit_data_path, file_list[idx]))
      summary_out_list[[idx]] <- out$performanceSummary
  
      wald_out <- TESSERA::calc_Wald_statistics(out, prep_obj$contrast_mat)
      wald_out$contrast_description <- rownames(wald_out)
      wald_out_list[[idx]] <- wald_out
    }
    summary_out_all <- dplyr::bind_rows(summary_out_list)
    saveRDS(summary_out_all, file.path(processed_data_path, paste0("model_results_perfsummary_", pattern_type, "_onlyinteract.rds")))
    wald_df_all <- dplyr::bind_rows(wald_out_list)
    saveRDS(wald_df_all, file.path(processed_data_path, paste0("model_results_waldstats_", pattern_type, "_onlyinteract.rds")))
  }
}

### GLMs and GAMs ###
# List of files to read in
pattern_type <- ".rds"
full_file_list <- base::list.files(glm_gam_modelfit_data_path, pattern = pattern_type)
file_list <- c()
for(f_idx in full_file_list){
  gene_name <- sub(".*gene_(.*)\\.rds", "\\1", f_idx)
  if(gene_name %in% gene_keep){
    file_list <- c(file_list, f_idx)
  }
}
length(file_list)

if(file.exists(file.path(processed_data_path, paste0("glm_gam_model_results_perfsummary_onlyinteract.rds")))
   & file.exists(file.path(processed_data_path, paste0("glm_gam_model_results_waldstats_onlyinteract.rds")))){
  summary_out_all <- readRDS(file.path(processed_data_path, paste0("glm_gam_model_results_perfsummary_onlyinteract.rds")))
  wald_df_all <- readRDS(file.path(processed_data_path, paste0("glm_gam_model_results_waldstats_onlyinteract.rds")))
}else{
  summary_out_list <- list()
  wald_out_list <- list()
  for(idx in 1:length(file_list)){
    if (0 == (idx %% 250)) {
      cat(idx, "\n")
    }

    out <- readRDS(paste0(glm_gam_modelfit_data_path, file_list[idx]))
    summary_out_list[[idx]] <- out$performanceSummary
    wald_out_list[[idx]] <- out$wald_df

  }
  summary_out_all <- dplyr::bind_rows(summary_out_list)
  saveRDS(summary_out_all, file.path(processed_data_path, paste0("glm_gam_model_results_perfsummary_onlyinteract.rds")))
  wald_df_all <- dplyr::bind_rows(wald_out_list)
  saveRDS(wald_df_all, file.path(processed_data_path, paste0("glm_gam_model_results_waldstats_onlyinteract.rds")))
}

