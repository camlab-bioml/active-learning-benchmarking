suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(yaml)
  library(data.table)
})
source("pipeline/whatsthatcell-helpers.R")

get_AL_probs <- function(df, AL_method){
  annotated_cells <- df %>% 
    filter(!is.na(cell_type)) %>% 
    filter(cell_type != "Skipped", cell_type != "Unclear")
  
  left_cells <- df %>% 
    filter(is.na(cell_type))
  
  
  bootstrap_reps <- 1000
  tctrl <- trainControl(method = "boot", number = bootstrap_reps)
  ModelFit <- train(cell_type ~ ., 
                    data = select(annotated_cells, -X1),
                    method = AL_method,
                    trace = FALSE,
                    trControl = tctrl)
  
  predicted_scores <- predict(ModelFit, 
                              select(left_cells, -X1, -cell_type),
                              type = "prob") |> 
    mutate(cell_id = left_cells$X1)
  predicted_scores
}


if(is.null(snakemake@params$rem_cell_type_list)){
  cell_type_to_rem <- snakemake@wildcards$rem_cell_type
}else{
  cell_type_to_rem <- snakemake@params$rem_cell_type_list
}

markers <- read_yaml(snakemake@input$markers)


# Create marker file with only markers for missing cell type
missing_cell_type_markers <- markers
marker_index_rem <- which(names(markers$cell_types) %in% cell_type_to_rem)
missing_cell_type_markers$cell_types <- markers$cell_types[marker_index_rem]

sce <- readRDS(snakemake@input$sce)
gt <- tibble(cell_id = colnames(sce),
             gt_cell_type = sce$CellType)

# Create PCA embedding and filtered expression
features <- create_features(sce)

### REMOVED CELL TYPE #####
# Create initial training set with removed cell type
selected_cells_rem_cell_type <- get_training_type_rem(features$expression, 
                                                      snakemake@wildcards$initial, 
                                                      markers,
                                                      cell_type_to_rem,
                                                      needed_cells = as.integer(snakemake@wildcards$num))

# Adds the cell type labels
missing_cell_type_PCA <- left_join(features$PCA, selected_cells_rem_cell_type)

rem_cell_type_probs <- get_AL_probs(
  select(missing_cell_type_PCA, -gt_cell_type, -iteration),
  snakemake@wildcards$AL_alg
)

### KEPT CELL TYPE #####
# Find all the cells of the type removed above that would otherwise have been selected 
# These will then be included in the training data
selected_cells_kept_df <- get_training_type_kept(features$expression, cell_type_to_rem, 
                                                 snakemake@wildcards$initial, 
                                                 missing_cell_type_markers, 
                                                 selected_cells_rem_cell_type)

# Run active learning with kept cell type
# Runs three times with 1, 2 and 3 cells of the removed type
kept_cells_uncertainty <- lapply(selected_cells_kept_df, function(x){
  df_PCA <- left_join(features$PCA, 
                      select(x, -num_missing_cells), 
                      by = "X1")
  
  kept_cell_type_probs <- get_AL_probs(
    select(df_PCA, -gt_cell_type, -iteration),
    snakemake@wildcards$AL_alg
  )
  
  kept_cell_type_probs |> 
    mutate(num_missing_cells = unique(x$num_missing_cells))
}) |> bind_rows()


rownames(rem_cell_type_probs) <- NULL
rem_cell_type_probs$num_missing_cells <- 0
rownames(kept_cells_uncertainty) <- NULL

probs <- bind_rows(
  pivot_longer(rem_cell_type_probs, -c(cell_id, num_missing_cells), 
               names_to = "predicted", values_to = "prob"),
  pivot_longer(kept_cells_uncertainty, -c(cell_id, num_missing_cells), 
               names_to = "predicted", values_to = "prob")
)

pdf(snakemake@output$pdf, height = 8, width = 12)
    probs |> 
    left_join(gt) |> 
    ggplot(aes(x = predicted, y = prob, fill = as.character(num_missing_cells))) +
    geom_boxplot() +
    labs(fill = "Num cells\nin training\ndataset",
        y = "Probability", x = "Predicted cell type",
        title = paste0("Removed cell type: ", cell_type_to_rem),
        subtitle = "Faceted by ground truth cell type") +
    facet_wrap(~gt_cell_type) +
    whatsthatcell_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()
