library(tidyverse)
library(yaml)
library(data.table)
library(caret)
library(SingleCellExperiment)
source("pipeline/whatsthatcell-helpers.R")


markers <- read_yaml(snakemake@input[['markers']])
unique_markers <- unique(unlist(markers$cell_types))

sce <- readRDS(snakemake@input[['expression']])
df_expression <- load_scs(sce)
df_expression$cell_type <- NA
df_expression$iteration <- NA
df_expression$gt_cell_type <- sce$CellType

### Create ground truth tibble
ground_truth <- tibble(cell_id = rownames(colData(sce)),
                       cell_type = sce$CellType)


remove_marker_AL <- function(markers, df_expression, ground_truth){
  subset_markers <- markers$cell_types
  
  # Get list of cell types
  all_cell_types <- names(subset_markers)
  
  iterations <- ((nrow(df_expression) - 2*length(all_cell_types)) / length(all_cell_types))
  iterations <- floor(iterations) - 1
  
  for(i in 1:iterations){
    # Get initial set of cells based on their marker expression ranking
    df_expression <- cell_ranking_wrapper(df_expression, markers)

    if(all(all_cell_types %in% unique(df_expression$cell_type))){
      break
    }
  }
  
  for(i in 1:20){
    AL <- active_learning_wrapper(df_expression, unique_markers, i, random_selection = 0)
    df_expression <- AL$expression
  }
  
  df_expression %>% 
    filter(!is.na(cell_type)) %>% 
    select(X1, iteration, cell_type, gt_cell_type) %>% 
    mutate(removed_markers = "all-kept") %>%
    mutate(seed = snakemake@wildcards[['seed']])
}


set.seed(as.integer(snakemake@wildcards[['seed']]))

cell_subset <- df_expression$X1[createDataPartition(df_expression$gt_cell_type, p = 0.3)$Resample1]

df <- filter(df_expression, X1 %in% cell_subset)
gt <- filter(ground_truth, cell_id %in% cell_subset)
#save.image('debug')
AL_removed <- remove_marker_AL(markers, df, gt)

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
