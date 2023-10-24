library(tidyverse)
library(splitstackshape)
library(SingleCellExperiment)
library(fabricatr)

set.seed(42)

seu <- read_tsv(snakemake@input[['seurat']])
sce <- readRDS(snakemake@input[['sce']])

sce_gt <- tibble(cell_id = colnames(sce),
                 cell_type = sce$CellType)

req_cells <- as.integer(snakemake@wildcards[['cell_num']])

n_cells <- 0

if(snakemake@params$selection_type == "wout_markers"){
  sel_method <- "NoMarkerSeurat-clustering"
  while(n_cells < req_cells){
    num_clusters <- seu$cluster %>% 
      unique() %>% 
      length()
    
    if(exists("sce_subset")){
      no_to_select <- ceiling((req_cells - ncol(sce_subset)) / num_clusters)
    }else{
      no_to_select <- ceiling(req_cells / num_clusters)
    }
    
    subset <- stratified(seu, "cluster", no_to_select)
    
    n_cells <- n_cells + nrow(subset)
    oversampled <- n_cells - req_cells
    while(oversampled > 0){
      rem_from_cluster <- subset |>
        group_by(cluster) |>
        tally() |>
        filter(n == max(n)) |> 
        slice_head(n = 1) |>
        pull(cluster)
        print(paste("rem from clust", rem_from_cluster))

      rem_cell_id <- subset |>
        filter(cluster == rem_from_cluster) |>
        slice_head(n = 1) |>
        pull(cell_id)

      subset <- filter(subset, cell_id != rem_cell_id)
      if(exists("sce_subset")){
        existing <- ncol(sce_subset)
      }else{
        existing <- 0
      }
      oversampled <- existing + nrow(subset) - req_cells
    }
    
    if(exists("sce_subset")){
      sce_subset <- cbind(sce_subset, sce[, colnames(sce) %in% subset$cell_id])
    }else{
      sce_subset <- sce[, colnames(sce) %in% subset$cell_id]
    }
    n_cells <- ncol(sce_subset)
    seu <- filter(seu, !(cell_id %in% colnames(sce_subset)))
  }
}else if(snakemake@params$selection_type == "w_markers"){
  sel_method <- "MarkerSeurat-clustering"
  while(n_cells < req_cells){
    num_clusters <- seu$predicted_cell_type %>% 
      unique() %>% 
      length()
    
    if(exists("sce_subset")){
      no_to_select <- ceiling((req_cells - ncol(sce_subset)) / num_clusters)
    }else{
      no_to_select <- ceiling(req_cells / num_clusters)
    }
    
    subset <- stratified(seu, "predicted_cell_type", no_to_select)
    
    n_cells <- n_cells + nrow(subset)
    oversampled <- n_cells - req_cells
    while(oversampled > 0){
      rem_from_cluster <- subset |>
        group_by(predicted_cell_type) |>
        tally() |>
        filter(n == max(n)) |> 
        slice_head(n = 1) |>
        pull(predicted_cell_type)

      rem_cell_id <- subset |>
        filter(predicted_cell_type == rem_from_cluster) |>
        slice_head(n = 1) |>
        pull(cell_id)

      subset <- filter(subset, cell_id != rem_cell_id)
      if(exists("sce_subset")){
        existing <- ncol(sce_subset)
      }else{
        existing <- 0
      }
      oversampled <- existing + nrow(subset) - req_cells
    }
    
    if(exists("sce_subset")){
      sce_subset <- cbind(sce_subset, sce[, colnames(sce) %in% subset$cell_id])
    }else{
      sce_subset <- sce[, colnames(sce) %in% subset$cell_id]
    }
    n_cells <- ncol(sce_subset)
    seu <- filter(seu, !(cell_id %in% colnames(sce_subset)))
  }
}


# sce_subset <- sce[, colnames(sce) %in% seu_subset$cell_id]

# while(n_cells < req_cells){
#   pred_cell_types <- seu$cluster %>% 
#     unique() %>% 
#     length()
  
#   if(exists("sce_subset")){
#     no_to_select <- ceiling((req_cells - ncol(sce_subset)) / pred_cell_types)
#   }else{
#     no_to_select <- ceiling(req_cells / pred_cell_types)
#   }
  
#   subset <- stratified(seu, "cluster", no_to_select) |> 
#     pull(cell_id)
  
#   n_cells <- n_cells + length(subset)
#   oversampled <- n_cells - req_cells
#   if(oversampled > 0){
#     subset <- subset[-sample(1:length(subset), oversampled)]
#   }
  
#   if(exists("sce_subset")){
#     sce_subset <- cbind(sce_subset, sce[, colnames(sce) %in% subset])
#   }else{
#     sce_subset <- sce[, colnames(sce) %in% subset]
#   }
#   seu <- filter(seu, !(cell_id %in% colnames(sce_subset)))
# }

# For the imbalanced analysis. Ensures there is at least one cell type of both types

if(!is.null(snakemake@wildcards$similarity)){
  if(length(unique(sce_subset$CellType)) < 2){
    sce_subset <- sce_subset[, -c(1:2)]

    each_ct <- filter(sce_gt, !(cell_id %in% colnames(sce_subset))) |>
      group_by(cell_type) |>
      slice_head(n = 1)
    
    sce_subset <- cbind(sce_subset, sce[, colnames(sce) %in% each_ct$cell_id])
  }
}

gt <- tibble(cell_id = colnames(sce_subset)) |>
  left_join(sce_gt) |>
  mutate(method = sel_method,
         params = paste0("knn-", snakemake@wildcards[['neighbors']], "-res-", 
                         snakemake@wildcards[['res']], '-cell-num-', 
                         snakemake@wildcards[['cell_num']]))

if(nrow(gt) != as.integer(snakemake@wildcards$cell_num)){
  stop("Not the correct number of cells selected")
}

saveRDS(sce_subset, snakemake@output[['sce']])
write_tsv(gt, snakemake@output[['ground_truth']])
