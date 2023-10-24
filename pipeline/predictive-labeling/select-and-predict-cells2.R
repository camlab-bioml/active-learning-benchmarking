suppressPackageStartupMessages({
  library(caret)
  library(SingleCellExperiment)
  library(yaml)
  library(data.table)
  library(tibble)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggalluvial)
})
source("pipeline/whatsthatcell-helpers.R")

### [ READ IN DATA ] #####
markers <- read_yaml(snakemake@input[['markers']]) %>% 
  unlist() %>% 
  unique()

sce <- readRDS(snakemake@input[['sce']])
if(any(grepl("CD16_32", rownames(sce)))){
  rownames(sce)[grep("CD16_32", rownames(sce))] <- "CD16-32"
}

# Get selected set of cells to subset to for training
selected_cells <- read_tsv(snakemake@input[['selected_cells']])


### [ PROCESS DATA ] #####
# create dataframe & add labels
df_expression <- load_scs(sce)
df_expression <- df_expression[, c(TRUE, colSums(df_expression[,2:ncol(df_expression)]) > 0)]

df_PCA <- select(df_expression, -X1) |> 
  as.matrix() |> 
  prcomp(center = TRUE, scale. = TRUE)

df_PCA <- df_PCA$x |> 
  as.data.frame()

df_PCA <- bind_cols(
  tibble(X1 = df_expression$X1),
  df_PCA[,1:20], 
  tibble(cell_type = sce$CellType)
)

# Separate into labelled and unlabeled 
labelled <- df_PCA %>% 
  filter(X1 %in% selected_cells$cell_id)

unlabeled <- df_PCA %>% 
  filter(!(X1 %in% selected_cells$cell_id))

### [ TRAIN AND PREDICT ] #####
# Train LR on labelled subset
multiNomModelFit <- train(cell_type ~ .,
                          data = select(labelled, -X1),
                          method = snakemake@wildcards[['pred_alg']],
                          trace = FALSE,
                          trControl = trainControl(method = "boot", number = 100))

# Predict probabilities
unlabeled_pred <- predict(multiNomModelFit,
                          select(unlabeled, -X1, -cell_type), 
                          type = "prob")

## Calculate entropy and max probability
entropy <- apply(unlabeled_pred, 1, calculate_entropy)
unlabeled_pred$entropy <- entropy

# Add predicted labels
unlabeled_pred$pred_cell_type <- predict(multiNomModelFit, select(unlabeled, -X1, -cell_type))

# Add cell ID
unlabeled_pred$cell_id <- unlabeled$X1

# Save full prediction table
bind_rows(
  mutate(selected_cells, entropy = NA, labeling = "gt") |> 
    select(cell_id, entropy, cell_type, labeling),
  
  select(unlabeled_pred, cell_id, entropy, pred_cell_type) |> 
    dplyr::rename("cell_type" = "pred_cell_type") |> 
    mutate(labeling = "pred") |> 
    select(cell_id, entropy, cell_type, labeling)
) |>
  write_tsv(snakemake@output[['full_pred']])

