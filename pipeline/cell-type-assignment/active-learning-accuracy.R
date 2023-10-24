suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(yaml)
  library(data.table)
  library(fabricatr)
  library(yardstick)
})
source("pipeline/whatsthatcell-helpers.R")

markers <- read_yaml(snakemake@input$markers)
sce <- readRDS(snakemake@input$sce)

if(any(grepl("CD16_32", rownames(sce)))){
  rownames(sce)[grep("CD16_32", rownames(sce))] <- "CD16-32"
}

if(grepl("maxp", snakemake@wildcards$strat)){
  criterion <- "maxp"
}else if(grepl("entropy", snakemake@wildcards$strat)){
  criterion <- "entropy"
}

df_expression <- load_scs(sce)
df_expression$cell_type <- NA
df_expression$gt_cell_type <- sce$CellType
df_expression$non_corrupted_gt_cell_type <- sce$CellType
df_expression$iteration <- NA
df_expression$corrupted <- NA

all_cell_types <- unique(sce$CellType)

## Corrupt a fraction of labels
if(as.numeric(snakemake@wildcards[['corrupt']]) != 0){
  df_expression <- corrupt_labels(df_expression, all_cell_types, snakemake@wildcards[['corrupt']])
}

### [ INITIAL CELL TYPE ASSIGNMENT ] #####
if(snakemake@wildcards[['initial']] == 'ranking'){
  df_expression <- cell_ranking_wrapper(df_expression, markers)
}else if(snakemake@wildcards[['initial']] == 'random'){
  random_cell_idx <- sample(1:nrow(df_expression), 20)
  df_expression$cell_type[random_cell_idx] <- df_expression$gt_cell_type[random_cell_idx]
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
         non_corrupted_gt_cell_type = df_expression$non_corrupted_gt_cell_type,
         iteration = df_expression$iteration)
)

iters <- 20

f1_scores <- tibble(
  iteration = 1:20,
  f1 = NA
)

for(i in 1:iters){
  print(i)
  annotated_cells <- df_PCA %>% 
    filter(!is.na(cell_type)) %>% 
    filter(cell_type != "Skipped", cell_type != "Unclear")
  print(dim(annotated_cells))
  
  left_cells <- df_PCA %>% 
    filter(is.na(cell_type))
  
  ModelFit <- fit_AL_classifier(select(annotated_cells, -gt_cell_type, -iteration, -non_corrupted_gt_cell_type),
                                snakemake@wildcards$AL_alg)
  
  ### Calculate F1-score
  predicted_scores <- predict(ModelFit, 
                              select(left_cells, -X1, -cell_type, -gt_cell_type, -iteration, -non_corrupted_gt_cell_type),
                              type = "raw")
  
  preds <- tibble(
    cell_id = left_cells$X1,
    annotated_cell_type = left_cells$non_corrupted_gt_cell_type,
    predicted_cell_type = predicted_scores
  )
  
  cell_types <- unique(union(preds$predicted_cell_type, preds$annotated_cell_type))
  preds$annotated_cell_type <- factor(preds$annotated_cell_type, levels = cell_types)
  preds$predicted_cell_type <- factor(preds$predicted_cell_type, levels = cell_types)
  
  f1 <- tryCatch(f_meas(preds, annotated_cell_type, predicted_cell_type), error=function(e) NULL) |> 
    pull(.estimate)
  f1_scores$f1[i] <- f1
  
  
  # Continue AL - predict probabilities
  predicted_scores <- predict(ModelFit, 
                              select(left_cells, -X1, -cell_type, -gt_cell_type, -iteration, -non_corrupted_gt_cell_type),
                              type = "prob")
  
  # Get next set of cells
  sel_cells <- entropy_maxp_cell_selection(criterion, predicted_scores, left_cells, 
                                           snakemake@wildcards$strat, 10, 0, annotated_cells)
  
  
  # What index do the selected cells correspond to?
  to_assign_index <- match(sel_cells[[1]], df_PCA$X1)
  
  # Get ground truth labels based on the index
  df_PCA$cell_type[to_assign_index] <- df_PCA$gt_cell_type[to_assign_index]
}


f1_scores |> 
  mutate(cohort = snakemake@wildcards$modality,
         params = paste0("Init-", snakemake@wildcards$initial,
                         "-strat-", snakemake@wildcards$strat,
                         "-AL_alg-", snakemake@wildcards$AL_alg,
                         "-rand-", snakemake@wildcards$rand,
                         "-corrupt-", snakemake@wildcards$corrupt,
                         "-seed-", snakemake@wildcards$s)) |> 
  write_tsv(snakemake@output$tsv)



