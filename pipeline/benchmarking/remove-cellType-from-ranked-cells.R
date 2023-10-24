library(tidyverse)
library(yaml)
library(data.table)
library(caret)
library(SingleCellExperiment)
library(fabricatr)
source("pipeline/whatsthatcell-helpers.R")

full_markers <- read_yaml("markers/CyTOF.yml")#snakemake@input[['markers']])

cell_type_to_rem <- "CLP"
marker_index_rem <- which(names(full_markers$cell_types) == cell_type_to_rem)

markers <- full_markers
markers$cell_types <- markers$cell_types[-marker_index_rem]
unique_markers <- unique(unlist(markers$cell_types))

sce <- readRDS("data/CyTOF/CyTOF-train.rds")#snakemake@input[['expression']])
if(TRUE %in% grepl("CD16_32", rownames(sce))){
  rownames(sce)[grep("CD16_32", rownames(sce))] <- "CD16-32"
}

df_expression <- load_scs(sce)
df_expression$cell_type <- NA
df_expression$iteration <- NA
df_expression$gt_cell_type <- sce$CellType

### Create ground truth tibble
ground_truth <- tibble(cell_id = rownames(colData(sce)),
                       cell_type = sce$CellType)

selected_cells <- 0
needed_cells <- 10
discarded_cells <- c()

ranked_cells <- tibble(
  X1 = character(),
  cell_type = character()
)

ranking_expression <- df_expression
while(selected_cells < needed_cells){
  left_to_rank_num <- needed_cells - selected_cells
  
  ranking_expression <- cell_ranking_wrapper(ranking_expression, 
                                             markers, 
                                             left_to_rank_num) |> 
    filter(cell_type != cell_type_to_rem  | is.na(cell_type))
  
  # Save selected cells minus cells to remove
  ranked_cells <- bind_rows(
    ranked_cells,
    filter(ranking_expression, !is.na(cell_type)) |> 
    select(X1, cell_type)
  )
  
  # Remove already selected cells so that new cells are selected with new ranking
  ranking_expression <- filter(ranking_expression, is.na(cell_type))
  
  # How many cells have been selected total?
  selected_cells <- nrow(ranked_cells)
}

## [ ACTIVE LEARNING ] ###
df_expression$corrupted <- NA
df_PCA <- select(df_expression, -c(X1, cell_type, iteration, gt_cell_type, corrupted)) |> 
  as.matrix() |> 
  prcomp(center = TRUE, scale. = TRUE)

df_PCA <- df_PCA$x |> 
  as.data.frame()

df_PCA <- bind_cols(
  tibble(X1 = df_expression$X1),
  df_PCA[,1:min(20, ncol(df_PCA))], 
  tibble(gt_cell_type = df_expression$gt_cell_type,
         iteration = df_expression$iteration)
) |> left_join(ranked_cells)

rem_cell_type_AL_wrapper <- function(df_PCA, iter = 6){
  entropies <- list()
  for(i in 1:iter){
    AL <- active_learning_wrapper(select(df_PCA, -gt_cell_type, -iteration), 
                                  'multinom', #snakemake@wildcards[['AL_alg']], 
                                  'highest_entropy', #snakemake@wildcards[['strat']], 
                                  i, 
                                  entropies, 
                                  0, #as.numeric(snakemake@wildcards[['rand']]),
                                  'entropy')#criterion)
    
    entropies[[length(entropies) + 1]] <- AL$criterion_table
    
    # What index do the selected cells correspond to?
    to_assign_index <- match(AL$new_cells, df_PCA$X1)
    
    # Get ground truth labels based on the index
    df_PCA$cell_type[to_assign_index] <- df_PCA$gt_cell_type[to_assign_index]
    df_PCA$iteration[to_assign_index] <- i
    
    not_annotated <- filter(df_PCA, is.na(cell_type)) %>% 
      nrow()
    print(i)
    if(not_annotated < 10){
      break
    }
  }
  
  entropies
}

rem_entropies <- rem_cell_type_AL_wrapper(df_PCA)

missing_cell_type_uncertainty <- bind_rows(rem_entropies) |> 
  as_tibble() |> 
  mutate(comp = "missing cell types",
         num_missing_cells = 0) |> 
  left_join(select(df_expression, X1, gt_cell_type), by = c("cell_id" = "X1"))



### [ REPEAT WITH ADDED CELL TYPE ] #####
missing_cell_type_markers <- full_markers
missing_cell_type_markers$cell_types <- full_markers$cell_types[marker_index_rem]
missing_df_expression <- filter(df_expression, gt_cell_type == cell_type_to_rem)

missing_df_expression <- cell_ranking_wrapper(missing_df_expression, 
                                              missing_cell_type_markers, 
                                              3) |> 
  filter(!is.na(cell_type)) |> 
  select(X1, cell_type)

include_ranked_cell <- 
  lapply(1:3, function(x){
  bind_rows(
    slice_head(ranked_cells, n = (nrow(ranked_cells) - x)),
    slice_head(missing_df_expression, n = x)
  ) |> mutate(num_missing_cells = x)
})


kept_cells_uncertainty <- lapply(include_ranked_cell, function(x){
  df_PCA$cell_type <- NULL
  df_PCA <- left_join(df_PCA, 
                      select(x, -num_missing_cells), 
                      by = "X1")
  
  entr <- rem_cell_type_AL_wrapper(df_PCA)
  
  bind_rows(entr) |> 
    as_tibble() |> 
    mutate(comp = "kept cell types",
           num_missing_cells = unique(x$num_missing_cells))
}) |> bind_rows() |> 
  left_join(select(df_expression, X1, gt_cell_type), by = c("cell_id" = "X1"))


bind_rows(missing_cell_type_uncertainty,
          kept_cells_uncertainty) |> 
  ggplot(aes(x = as.character(no_cells_annotated), y = criterion_val, fill = as.character(num_missing_cells))) + 
  geom_boxplot() + 
  facet_wrap(~gt_cell_type)



### OLD #####
remove_marker_AL <- function(markers, remove_index, df_expression, ground_truth){
  subset_markers <- markers$cell_types[-remove_index]
  
  # Get list of cell types
  all_cell_types <- names(subset_markers)
  
  iterations <- ((nrow(df_expression) - 2*length(all_cell_types)) / length(all_cell_types))
  iterations <- floor(iterations) - 1
  
  for(i in 1:iterations){
    # Get initial set of cells based on their marker expression ranking
    df_expression <- cell_ranking_wrapper(df_expression, markers, random_selection = 0)

    if(all(all_cell_types %in% unique(df_expression$cell_type))){
      break
    }
  }
  
  # If any cell types have been assigned to the removed cell type make them NA again
  df_expression$cell_type[which(df_expression$cell_type == names(markers$cell_types[remove_index]))] <- NA
  
  for(i in 1:20){
    AL <- active_learning_wrapper(df_expression, unique_markers, i)
    df_expression <- AL$expression
  }
  
  df_expression %>% 
    filter(!is.na(cell_type)) %>% 
    select(X1, iteration, cell_type, gt_cell_type) %>% 
    mutate(removed_markers = names(markers$cell_types[remove_index])) %>%
    mutate(seed = snakemake@wildcards[['seed']])
}


set.seed(as.integer(snakemake@wildcards[['seed']]))

cell_subset <- df_expression$X1[createDataPartition(df_expression$gt_cell_type, p = 0.3)$Resample1]

df <- filter(df_expression, X1 %in% cell_subset)
gt <- filter(ground_truth, cell_id %in% cell_subset)

AL_removed <- lapply(1:length(markers$cell_types), function(x){
  remove_marker_AL(markers, x, df, ground_truth)
}) %>% 
  bind_rows()

write_tsv(AL_removed, snakemake@output[['tsv']])

# AL_removed %>%
#   bind_rows() %>% 
#   group_by(removed_markers, iteration, cell_type) %>% 
#   tally() %>% 
#   ggplot(aes(x = iteration, y = n, fill = cell_type)) +
#   geom_bar(stat = 'identity') +
#   scale_fill_manual(values = cell_type_colours()) +
#   facet_wrap(~removed_markers) +
#   whatsthatcell_theme()
