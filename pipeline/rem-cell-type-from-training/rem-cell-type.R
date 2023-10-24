suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(data.table)
  library(caret)
  library(SingleCellExperiment)
  library(fabricatr)
})
source("pipeline/whatsthatcell-helpers.R")
set.seed(1)

markers <- read_yaml(snakemake@input$markers)
if(grepl('entropy', snakemake@wildcards$strat)){
  criterion <- 'entropy'
}else{
  criterion <- 'maxp'
}

# Create marker file with only markers for missing cell type
missing_cell_type_markers <- markers

if(is.null(snakemake@params$rem_cell_type_list)){
  cell_type_to_rem <- snakemake@wildcards$rem_cell_type
}else{
  cell_type_to_rem <- snakemake@params$rem_cell_type_list
}

marker_index_rem <- which(names(markers$cell_types) %in% cell_type_to_rem)
missing_cell_type_markers$cell_types <- markers$cell_types[marker_index_rem]

sce <- readRDS(snakemake@input$sce)

# Create PCA embedding and filtered expression
features <- create_features(sce)

### REMOVED CELL TYPE #####
# Create initial training set with removed cell type
selected_cells_rem_cell_type <- get_training_type_rem(features$expression, 
                                                      snakemake@wildcards$initial, 
                                                      markers,
                                                      cell_type_to_rem,
                                                      as.integer(snakemake@wildcards$num))

if(nrow(selected_cells_rem_cell_type) != as.integer(snakemake@wildcards$num)){
  stop("Number of cells selected is not equal to the number of cells requested")
}

# Adds the cell type labels
missing_cell_type_PCA <- left_join(features$PCA, selected_cells_rem_cell_type)

# Run active learning with removed cell type
missing_cell_type_uncertainty <- rem_cell_type_AL_wrapper(missing_cell_type_PCA,
                                                          snakemake@wildcards$AL_alg,
                                                          snakemake@wildcards$strat,
                                                          rand = 0,
                                                          criterion = criterion) |>
  bind_rows() |> 
  as_tibble() |> 
  mutate(comp = "missing cell types",
         num_missing_cells = 0) |> 
  left_join(select(features$expression, X1, gt_cell_type), by = c("cell_id" = "X1")) # adds gt cell type



### KEPT CELL TYPE #####
# Find all the cells of the type removed above that would otherwise have been selected 
# These will then be included in the training data
selected_cells_kept_df <- get_training_type_kept(features$expression, cell_type_to_rem, 
                                                 snakemake@wildcards$initial, 
                                                 missing_cell_type_markers, 
                                                 selected_cells_rem_cell_type)

num_cells <- lapply(selected_cells_kept_df, nrow) |>
  unlist()
if(any(num_cells != as.integer(snakemake@wildcards$num))){
  stop("Number of cells selected is not equal to the number of cells requested")
}
# Run active learning with kept cell type
# Runs three times with 1, 2 and 3 cells of the removed type
kept_cells_uncertainty <- lapply(selected_cells_kept_df, function(x){
  df_PCA <- left_join(features$PCA, 
                      select(x, -num_missing_cells), 
                      by = "X1")
  
  entr <- rem_cell_type_AL_wrapper(df_PCA,
                                   snakemake@wildcards$AL_alg,
                                   snakemake@wildcards$strat,
                                   rand = 0,
                                   criterion = criterion)
  
  bind_rows(entr) |> 
    as_tibble() |> 
    mutate(comp = "kept cell types",
           num_missing_cells = unique(x$num_missing_cells))
}) |> bind_rows() |> 
  left_join(select(features$expression, X1, gt_cell_type), by = c("cell_id" = "X1"))


# Combine
uncertainties <- bind_rows(
  missing_cell_type_uncertainty,
  kept_cells_uncertainty
  ) |> 
  mutate(params = paste0("AL_alg-", snakemake@wildcards$AL_alg, "-strat-",
                         snakemake@wildcards$strat, "-init-", 
                         snakemake@wildcards$initial, "-rem_celltype-",
                         paste0(cell_type_to_rem, collapse = ", "), "-seed-", snakemake@wildcards$s))

write_tsv(uncertainties, snakemake@output$tsv)
