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
                          trace = FALSE)

# Predict probabilities
unlabeled_pred <- predict(multiNomModelFit,
                          select(unlabeled, -X1, -cell_type), 
                          type = "prob")

## Calculate entropy and max probability
entropy <- apply(unlabeled_pred, 1, calculate_entropy)

max_p_idx <- apply(unlabeled_pred, 1, which.max)
max_p <- lapply(seq_len(nrow(unlabeled_pred)), function(x){
  unlabeled_pred[x, max_p_idx[x]]
}) %>% unlist()

unlabeled_pred$entropy <- entropy
unlabeled_pred$max_p <- max_p

# Add predicted labels
unlabeled_pred$pred_cell_type <- predict(multiNomModelFit, select(unlabeled, -X1, -cell_type))

# Add cell ID
unlabeled_pred$cell_id <- unlabeled$X1

# Save full prediction table
bind_rows(
  mutate(selected_cells, entropy = NA, max_p = NA, labeling = "gt") |> 
    select(cell_id, entropy, max_p, cell_type, labeling),
  
  select(unlabeled_pred, cell_id, entropy, max_p, pred_cell_type) |> 
    dplyr::rename("cell_type" = "pred_cell_type") |> 
    mutate(labeling = "pred") |> 
    select(cell_id, entropy, max_p, cell_type, labeling)
) |>
  write_tsv(snakemake@output[['full_pred']])

# Remove probabilities for each cell type
unlabeled_pred <- select(unlabeled_pred, cell_id, pred_cell_type, entropy, max_p)

# Calculate number of cells to select
additional_cell_num <- round(nrow(selected_cells) * 0.20)

entropy_cells <- unlabeled_pred %>% 
  arrange(entropy) %>% 
  slice_head(n = additional_cell_num) %>% 
  select(cell_id, pred_cell_type) %>% 
  mutate(selection_strat = "entropy") %>% 
  as_tibble()

maxp_cells <- unlabeled_pred %>% 
  arrange(-max_p) %>% 
  slice_head(n = additional_cell_num) %>% 
  select(cell_id, pred_cell_type) %>% 
  mutate(selection_strat = "maxp") %>% 
  as_tibble()

high_prob_cells <- unlabeled_pred %>%
  filter(max_p > 0.99) %>%
  select(cell_id, pred_cell_type) %>% 
  mutate(selection_strat = "99p") %>% 
  as_tibble()

entropy_maxp <- bind_rows(entropy_cells, maxp_cells, high_prob_cells)

pdf(snakemake@output[['pdf']])
  entropy_maxp %>% 
    left_join(select(unlabeled, X1, cell_type), by = c("cell_id" = "X1")) %>% 
    dplyr::rename("ground_truth" = "cell_type") %>% 
    pivot_longer(c(pred_cell_type, ground_truth), names_to = "type", values_to = "cell_type") %>% 
    ggplot(aes(x = type, stratum = cell_type, fill = cell_type, alluvium = cell_id)) +
    geom_stratum() +
    geom_flow() +
    scale_fill_manual(values = cell_type_colours(snakemake@wildcards[['modality']])) +
    whatsthatcell_theme() +
    facet_wrap(~selection_strat, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()

entropy_cells <- select(entropy_cells, -selection_strat) %>%
  dplyr::rename("cell_type" = "pred_cell_type")
maxp_cells <- select(maxp_cells, -selection_strat) %>%
  dplyr::rename("cell_type" = "pred_cell_type")
highp_cells <- select(high_prob_cells, -selection_strat) %>% 
  dplyr::rename("cell_type" = "pred_cell_type")

if(snakemake@wildcards[['selection_procedure']] == "Active-Learning_entropy" |
   snakemake@wildcards[['selection_procedure']] == "Active-Learning_maxp"){
  
  entropy_cells$corrupted_cell_type <- NA
  entropy_cells$iteration <- 'predicted'
  entropy_cells$method <- paste0(unique(selected_cells$method), "-predLabelSelection-entropy-AL_alg-", snakemake@wildcards[['AL_alg']])
  entropy_cells$cell_num <- additional_cell_num
  
  maxp_cells$corrupted_cell_type <- NA
  maxp_cells$iteration <- 'predicted'
  maxp_cells$method <- paste0(unique(selected_cells$method), "-predLabelSelection-maxp-AL_alg-", snakemake@wildcards[['AL_alg']])
  entropy_cells$cell_num <- additional_cell_num

  highp_cells$corrupted_cell_type <- NA
  highp_cells$iteration <- 'predicted'
  highp_cells$method <- paste0(unique(selected_cells$method), "-predLabelSelection-99p-AL_alg-", snakemake@wildcards[['AL_alg']])
  highp_cells$cell_num <- nrow(highp_cells)

  selected_cells$iteration <- as.character(selected_cells$iteration)
}else if(snakemake@wildcards[['selection_procedure']] == "random"){
  entropy_cells$params <- paste0(unique(selected_cells$params), "-predLabelSelection-entropy-AL_alg-", snakemake@wildcards[['AL_alg']])
  maxp_cells$params <- paste0(unique(selected_cells$params), "-predLabelSelection-maxp-AL_alg-", snakemake@wildcards[['AL_alg']])
  highp_cells$params <- paste0(unique(selected_cells$params), "-predLabelSelection-99p-AL_alg-", snakemake@wildcards[['AL_alg']])

}else if(grepl("Seurat-clustering", snakemake@wildcards[['selection_procedure']])){
  entropy_cells$method <- snakemake@wildcards[['selection_procedure']]
  maxp_cells$method <- snakemake@wildcards[['selection_procedure']]
  highp_cells$method <- snakemake@wildcards[['selection_procedure']]

  entropy_cells$params <- paste0(unique(selected_cells$params), "-predLabelSelection-maxp-AL_alg-", snakemake@wildcards[['AL_alg']])
  maxp_cells$params <- paste0(unique(selected_cells$params), "-predLabelSelection-entropy-AL_alg-", snakemake@wildcards[['AL_alg']])
  highp_cells$params <- paste0(unique(selected_cells$params), "-predLabelSelection-99p-AL_alg-", snakemake@wildcards[['AL_alg']])
}

if(nrow(entropy_cells) > 0){
  entropy_cells <- bind_rows(selected_cells, entropy_cells)
}else{
  entropy_cells <- selected_cells
}
if(nrow(maxp_cells) > 0){
  maxp_cells <- bind_rows(selected_cells, maxp_cells)
}else{
  maxp_cells <- selected_cells
}
if(nrow(highp_cells) > 0){
  highp_cells <- bind_rows(selected_cells, highp_cells)
}else{
  highp_cells <- selected_cells
}

write_tsv(entropy_cells, snakemake@output[['entropy_cells']])
write_tsv(maxp_cells, snakemake@output[['maxp_cells']])
write_tsv(highp_cells, snakemake@output[['high_prob_cells']])

