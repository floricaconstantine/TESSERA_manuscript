#' Create pseudobulk and bulk data from single cell data
#' @author Florica Constantine
#' Dependencies in file: dplyr.

## Create pseudobulk and bulk data
make_pseudobulk_bulk <- function(count_data,
                                 meta_data,
                                 individual_colname,
                                 cell_type_colname,
                                 cell_id_colname) {
  require(dplyr)
  
  # Extract metadata as vectors aligned with the count matrix columns
  individual_vector <- meta_data[, individual_colname]
  cell_type_vector <- meta_data[, cell_type_colname]
  
  # Unique individuals
  individuals <- unique(individual_vector)
  # Unique cell types
  cell_types <- unique(cell_type_vector)
  
  # Pseudobulk
  pseudobulk_data <- matrix(0, nrow = nrow(count_data), ncol = 0)
  rownames(pseudobulk_data) <- rownames(count_data)
  pseudobulk_meta <- data.frame(cell_type_colname = character(),
                                individual_colname = numeric())
  colnames(pseudobulk_meta) <- c(cell_type_colname, individual_colname)
  
  # Bulk
  bulk_data <- matrix(0, nrow = nrow(count_data), ncol = 0)
  rownames(bulk_data) <- rownames(count_data)
  bulk_meta <- data.frame(individual_colname = numeric())
  colnames(bulk_meta) <- c(individual_colname)
  
  # Loop over individual/cell-type combinations
  for (i_idx in 1:length(individuals)) {
    # Individual
    i_subset <- which(individual_vector == individuals[i_idx])
    
    for (c_idx in 1:length(cell_types)) {
      # Cell Type
      c_subset <- which(cell_type_vector == cell_types[c_idx])
      
      # Make pseudobulk data
      subset_indices <- intersect(i_subset, c_subset)
      if (0 < length(subset_indices)) {
        pseudobulk_data <- cbind(pseudobulk_data, rowSums(as.matrix(count_data[, subset_indices])))
        tmp <- data.frame(cell_type_colname = cell_types[c_idx],
                          individual_colname = individuals[i_idx])
        colnames(tmp) <- c(cell_type_colname, individual_colname)
        pseudobulk_meta <- rbind(pseudobulk_meta, tmp)
      }
    }
    
    # Make bulk data
    bulk_data <- cbind(bulk_data, rowSums(count_data[, i_subset]))
    tmp <- data.frame(individual_colname = individuals[i_idx])
    colnames(tmp) <- c(individual_colname)
    bulk_meta <- rbind(bulk_meta, tmp)
  }
  pseudobulk_meta$pb_id <- paste0(pseudobulk_meta[, individual_colname], "_", pseudobulk_meta[, cell_type_colname])
  colnames(pseudobulk_data) <- paste0(pseudobulk_meta[, individual_colname], "_", pseudobulk_meta[, cell_type_colname])
  
  # Get rest of metadata for pseudobulk
  pseudobulk_meta$order_internal <- seq(1, nrow(pseudobulk_meta))
  tmp1 <- as.data.frame(meta_data
                        %>% dplyr::group_by_at(c(cell_type_colname, individual_colname))
                        %>% summarise_if(is.numeric, mean))
  tmp2 <- as.data.frame(meta_data[, which(colnames(meta_data) != cell_id_colname)]
                        %>% dplyr::group_by_at(c(cell_type_colname, individual_colname))
                        %>% select_if(~ !is.numeric(.x))
                        %>% distinct())
  pseudobulk_meta <- merge(
    pseudobulk_meta,
    merge(tmp1, tmp2, by = c(cell_type_colname, individual_colname)),
    by = c(cell_type_colname, individual_colname)
  )
  pseudobulk_meta <- pseudobulk_meta[order(pseudobulk_meta$order_internal), ]
  
  # Get rest of metadata for bulk
  bulk_meta$order_internal <- seq(1, nrow(bulk_meta))
  tmp1 <- as.data.frame(meta_data
                        %>% group_by_at(c(individual_colname))
                        %>% summarise_if(is.numeric, mean))
  tmp2 <- as.data.frame(meta_data[, which((colnames(meta_data) != cell_id_colname)
                                          &
                                            (colnames(meta_data) != cell_type_colname))]
                        %>% group_by_at(c(individual_colname))
                        %>% select_if(~ !is.numeric(.x))
                        %>% distinct())
  bulk_meta <- merge(bulk_meta, merge(tmp1, tmp2, by = c(individual_colname)), by =
                       c(individual_colname))
  bulk_meta <- bulk_meta[order(bulk_meta$order_internal), ]
  
  # Return pseudobulk and bulk data and associated metadata
  return(
    list(
      pseudobulk_meta = pseudobulk_meta,
      pseudobulk_data = pseudobulk_data,
      bulk_meta = bulk_meta,
      bulk_data = bulk_data
    )
  )
}
