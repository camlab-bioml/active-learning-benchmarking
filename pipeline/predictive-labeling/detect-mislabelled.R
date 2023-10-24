# Can predictive labeling find mislabelled cell types?
library(data.table)
library(readr)
library(tibble)
library(dplyr)
library(caret)
library(ggalluvial)
library(tidyr)
source("pipeline/whatsthatcell-helpers.R")
set.seed(42)

sce <- readRDS(snakemake@input$sce)

# Create features
features <- create_features(sce)
train_features <- features$PCA |> 
  select(-iteration)

# Get training set
train_features <- train_features[sample(1:nrow(train_features), 250), ]
train_features$cell_type <- train_features$gt_cell_type

# Corrupt 10% of cells to random label
all_cell_types <- unique(features$PCA$gt_cell_type)
corrupt_idx <- sample(1:nrow(train_features), 25)

for(i in corrupt_idx){
  gt_cell_type <- train_features$cell_type[i]
  list_corrupt_ct <- all_cell_types[!all_cell_types %in% gt_cell_type]
  
  train_features$cell_type[i] <- sample(list_corrupt_ct, 1)
}

# Train initial model
multiNomModelFit <- train(cell_type ~ .,
                          data = select(train_features, -X1, -gt_cell_type),
                          method = snakemake@wildcards[['pred_alg']],
                          trace = FALSE)


# Predict probabilities
unlabeled_pred <- predict(multiNomModelFit,
                          select(train_features, -X1, -gt_cell_type, -cell_type), 
                          type = "prob")
unlabeled_pred$pred_type <- predict(multiNomModelFit,
                                    select(train_features, -X1, -gt_cell_type, -cell_type))

unlabeled_pred$cell_id <- train_features$X1
unlabeled_pred$corr_cell_type <- train_features$cell_type
unlabeled_pred$gt_cell_type <- train_features$gt_cell_type

unlabeled_pred$params <- paste0("mod-", snakemake@wildcards$modality,
                                "-predLabAlg-", snakemake@wildcards[['pred_alg']],
                                "-seed-", snakemake@wildcards[['seed']])

write_tsv(unlabeled_pred, snakemake@output$pred)

# 
# unlabeled_pred |> 
#   select(cell_id, gt_cell_type, corr_cell_type, pred_type) |> 
#   filter(gt_cell_type != corr_cell_type) |> 
#   pivot_longer(-cell_id, values_to = "cell_type", names_to = "category") |> 
#   mutate(category = factor(category, levels = c("gt_cell_type", "corr_cell_type", "pred_type"))) |> 
#   ggplot(aes(x = category, alluvium = cell_id, stratum = cell_type, fill = cell_type)) +
#   geom_stratum() +
#   stat_flow() +
#   whatsthatcell_theme()



