# How accurate can a random forest be & what do the probabilities look like when trained on lots of data?
suppressPackageStartupMessages({
  library(caret)
  library(data.table)
  library(patchwork)
  library(tidyverse)
})
source("pipeline/whatsthatcell-helpers.R")

# Create PCA embedding and filtered expression
sce <- readRDS(snakemake@input$sce)
features <- create_features(sce)

gt <- tibble(cell_id = colnames(sce),
             cell_type = sce$cell_type)

# Join cell types and remove one cell type
pc_features <- left_join(features$PCA, gt, by = c("X1" = "cell_id")) |> 
  filter(cell_type != snakemake@wildcards$rem)

tctrl <- trainControl(method='cv', 
                      number=5, 
                      repeats=3)


tunegrid <- expand.grid(mtry=c(1:20))

tunegrid
set.seed(1)
# First cross validation
ModelFit <- train(cell_type ~ ., 
                  data = select(pc_features, -X1, -gt_cell_type, -iteration),
                  method = 'rf',
                  tuneGrid=tunegrid, 
                  trControl = tctrl)


# Now re-train model on full dataset with best mtry
tctrl <- trainControl(method='none')
tunegrid <- data.frame(mtry = ModelFit$bestTune$mtry)
ModelFit <- train(cell_type ~ .,
                  data = select(pc_features, -X1, -gt_cell_type, -iteration),
                  method = 'rf',
                  tuneGrid=tunegrid, 
                  trControl = tctrl)

### Test on held out test set
test <- readRDS(snakemake@input$test)
test_features <- create_features(test)

predicted_scores <- predict(ModelFit, 
                            select(test_features$PCA, -gt_cell_type, -iteration),
                            type = "prob")

predicted_scores$entropy <- apply(predicted_scores, 1, calculate_entropy)
predicted_scores$cell_id <- test_features$PCA$X1
predicted_scores$cell_type <- test_features$PCA$gt_cell_type

probs <-  predicted_scores |> 
  pivot_longer(-c(cell_id, cell_type, entropy), 
               names_to =  "predicted", values_to = "prob") |> 
  ggplot(aes(x = predicted, y = prob)) +
  geom_boxplot() +
  facet_wrap(~cell_type, nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

entropy <- predicted_scores |> 
  ggplot(aes(x = cell_type, y = entropy)) +
  geom_boxplot() +
  labs(x = "ground truth cell type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

pdf(snakemake@output$probs_entropy, width = 10, height = 6)
  probs / entropy + plot_layout(heights = c(2, 1))
dev.off()
