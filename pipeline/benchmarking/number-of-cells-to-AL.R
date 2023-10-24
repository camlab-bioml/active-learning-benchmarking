
library(tidyverse)
library(yaml)
source("pipeline/whatsthatcell-helpers.R")


markers <- read_yaml("markers/scRNASeq.yml")#snakemake@input[['markers']])
unique_markers <- unique(unlist(markers$cell_types))

sce <- readRDS("data/scRNASeq/scRNASeq-train.rds")#snakemake@input[['expression']])
df_expression <- load_scs(sce)
df_expression$cell_type <- NA
df_expression$iteration <- NA
df_expression$gt_cell_type <- sce$cell_type

# Get list of cell types
all_cell_types <- ground_truth$cell_type %>% unique()


get_preranked_cell_no <- function(df, prop, all_cell_types, seed){
  set.seed(seed)
  subset <- createDataPartition(df$gt_cell_type, p = prop)$Resample1
  df <- df[subset, ]
  
  iterations <- ((nrow(df) - 2*length(all_cell_types)) / length(all_cell_types))
  iterations <- floor(iterations) - 1

  ### [ RANKED CELL TYPE ASSIGNMENT ] #####
  for(i in 1:iterations){
    # Get initial set of cells based on their marker expression ranking
    ranked_cells <- select_initial_cells(df, markers$cell_types)
    
    # What index do the selected cells correspond to?
    to_assign_index <- match(ranked_cells, df$X1)
   
    # Get ground truth labels based on the index
    assignment <- df$gt_cell_type[to_assign_index]
    
    df$cell_type[to_assign_index] <- assignment
    df$iteration[to_assign_index] <- i
    
    if(all(all_cell_types %in% unique(df$cell_type))){
      break
    }
  }
  cells_to_completion <- df %>% 
    filter(!is.na(cell_type))%>% 
    nrow()
  
  frequencies <- df$gt_cell_type %>% 
    table %>% 
    as.data.frame %>% 
    arrange(Freq)
  
  list(a = tibble(nrow(df), cells_to_completion, prop, seed,
         frequencies$.[1], frequencies$Freq[1]),
       b = filter(df, !is.na(cell_type)))
}

  
a <- lapply(seq(0.1, 1, by = 0.1), function(i){
  lapply(1:8, function(j){
    get_preranked_cell_no(df_expression, i, all_cell_types, j)
  })
})

lapply(1:length(a), function(x) a[x][[1]])

b <- bind_rows(a[[1]])
colnames(b) <- c("total_cells_in_subset", "cells_to_completion", 
                 "proportion_of_cells_selected", "seed", 
                 "cell_type_of_lowest_frequency", "lowest_frequency")

a %>% 
  ggplot(aes(x = as.character(proportion_of_cells_selected), y = cells_to_completion)) +
  geom_point()


a %>% 
  ggplot(aes(x = total_cells_in_subset, y = cells_to_completion)) +
  geom_point()


a %>% 
  ggplot(aes(x = lowest_frequency, y = cells_to_completion)) +
  geom_point()


  