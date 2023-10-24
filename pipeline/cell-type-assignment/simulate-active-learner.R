library(caret)
library(data.table)
library(yaml)
library(SingleCellExperiment)
library(tidyverse)
library(fabricatr)
source("pipeline/whatsthatcell-helpers.R")
set.seed(42)

### [ PARAMS ] #####
max_AL_iterations <- as.integer(snakemake@params[['max_cell_num']]) / 10

if(snakemake@wildcards[['AL_type']] == 'Active-Learning_maxp'){
  criterion <- "maxp"
}else if(snakemake@wildcards[['AL_type']] == 'Active-Learning_entropy'){
  criterion <- "entropy"
}

### [ LOAD & PROCESS DATA ] #####
markers <- read_yaml(snakemake@input[['markers']])

sce <- readRDS(snakemake@input[['expression']])
if(any(grepl("CD16_32", rownames(sce)))){
  rownames(sce)[grep("CD16_32", rownames(sce))] <- "CD16-32"
}

df_expression <- load_scs(sce)
df_expression$cell_type <- NA
df_expression$gt_cell_type <- sce$CellType
df_expression$iteration <- NA
df_expression$corrupted <- NA

all_cell_types <- unique(sce$CellType)

## Corrupt a fraction of labels
if(as.numeric(snakemake@wildcards[['corrupt']]) != 0){
  df_expression <- corrupt_labels(df_expression, all_cell_types, snakemake@wildcards[['corrupt']])
}

total_initial_cells <- 20

## Assign one of each cell type for classifiers to work in imbalanced setting
if(!is.null(snakemake@wildcards[['similarity']])){
  cell_types <- unique(df_expression$gt_cell_type)

  one_of_each_cells <- lapply(cell_types, function(x){
    select(df_expression, X1, gt_cell_type) |> 
      filter(gt_cell_type == x) |> 
      slice_head(n = 1) |> 
      pull(X1)
  }) |> unlist()

  one_of_each_idx <- which(df_expression$X1 %in% one_of_each_cells)
  df_expression$cell_type[one_of_each_idx] <- df_expression$gt_cell_type[one_of_each_idx]
  df_expression$iteration[one_of_each_idx] <- 0

  total_initial_cells <- total_initial_cells - length(one_of_each_idx)
}

### [ INITIAL CELL TYPE ASSIGNMENT ] #####
if(snakemake@wildcards[['initial']] == 'ranking'){
  df_expression <- cell_ranking_wrapper(df_expression, markers)
}else if(snakemake@wildcards[['initial']] == 'random'){
  random_cell_idx <- sample(1:nrow(df_expression), total_initial_cells)
  df_expression$cell_type[random_cell_idx] <- df_expression$gt_cell_type[random_cell_idx]
  df_expression$iteration[random_cell_idx] <- 0
}

### [ ACTIVE LEARNING CELL TYPE ASSIGNMENT ] #####
# create a list to save all entropies into
entropies <- list()

# Remove genes with 0 expression
df_expression_non_zero <- select(df_expression, -where(is.character), -iteration, -cell_type, -corrupted)
df_expression_non_zero <- df_expression_non_zero[,colSums(df_expression_non_zero) > 0]

# Calculate PCA embedding
df_PCA <- as.matrix(df_expression_non_zero) |> 
  prcomp(center = TRUE, scale. = TRUE)

df_PCA <- df_PCA$x |> 
  as.data.frame()

df_PCA <- bind_cols(
  tibble(X1 = df_expression$X1),
  df_PCA[,1:min(20, ncol(df_PCA))], 
  tibble(cell_type = df_expression$cell_type,
         gt_cell_type = df_expression$gt_cell_type,
         iteration = df_expression$iteration)
)

table(df_PCA$cell_type)

for(i in 1:max_AL_iterations){
  AL <- active_learning_wrapper(select(df_PCA, -gt_cell_type, -iteration), 
                                snakemake@wildcards[['AL_alg']], 
                                snakemake@wildcards[['strat']], 
                                i, 
                                entropies, 
                                as.numeric(snakemake@wildcards[['rand']]),
                                criterion)
  if(length(AL$new_cells) != 10){
    stop("Error: Active learning did not select the right number of cells") # this is to prevent issues with quantile selection (when quantiles are equal to max no cells were selected)
  }

  entropies[[length(entropies) + 1]] <- AL$criterion_table
  
  # What index do the selected cells correspond to?
  to_assign_index <- match(AL$new_cells, df_PCA$X1)
  
  # Get ground truth labels based on the index
  df_PCA$cell_type[to_assign_index] <- df_PCA$gt_cell_type[to_assign_index]
  print(table(df_PCA$cell_type))
  df_PCA$iteration[to_assign_index] <- i
  
  not_annotated <- filter(df_PCA, is.na(cell_type)) %>% 
    nrow()
  print(i)
  if(not_annotated < 10){
    break
  }
}


### [ SAVE OUTPUT ] #####
entropies %>% 
  bind_rows() %>% 
  as_tibble() %>% 
  mutate(method = paste0("Active-Learning-groundTruth-initial_sel-", 
                         snakemake@wildcards[['initial']], '-strategy-',
                         snakemake@wildcards[['strat']], '-AL_alg-', 
                         snakemake@wildcards[['AL_alg']], '-randomCells-', 
                         snakemake@wildcards[['rand']], '-corrupted-', 
                         snakemake@wildcards[['corrupt']], '-seed-',
                         snakemake@wildcards[['s']])) %>% 
  write_tsv(snakemake@output[['entropy']])

original_cell_types <- tibble(cell_id = colnames(sce),
                              cell_type = sce$CellType)

df_PCA %>% 
  select(X1, cell_type, iteration) %>%
  dplyr::rename("cell_id" = "X1", 
                "corrupted_cell_type" = "cell_type") %>% # this is the field iteratively filled in by AL, thus may be corrupted
  left_join(original_cell_types, by = 'cell_id') %>% 
  filter(!is.na(corrupted_cell_type)) %>% 
  mutate(method = paste0("Active-Learning-groundTruth-", 
                         snakemake@wildcards[['initial']], '-strategy-',
                         snakemake@wildcards[['strat']], '-AL_alg-', 
                         snakemake@wildcards[['AL_alg']], '-randomCells-', 
                         snakemake@wildcards[['rand']], '-corrupted-', 
                         snakemake@wildcards[['corrupt']], 'seed-',
                         snakemake@wildcards[['s']])) %>% 
  write_tsv(snakemake@output[['assignments']])