library(tidyverse)
library(yardstick)

### [ LOAD GROUND TRUTH ] #####
scRNAseq <- readRDS("data/scRNASeq/scRNASeq-test.rds")
scRNAseq_gt <- tibble(cell_id = colnames(scRNAseq),
                      annotated_cell_type = scRNAseq$CellType)

CyTOF <- readRDS("data/CyTOF/CyTOF-test.rds")
CyTOF_gt <- tibble(cell_id = colnames(CyTOF),
                   annotated_cell_type = CyTOF$CellType)


### [ FUNCTIONS ] #####
acc_wrap <- function(tt) {
  cell_types <- unique(union(tt$predicted_cell_type, tt$annotated_cell_type))
  
  tt$annotated_cell_type <- factor(tt$annotated_cell_type, levels = cell_types)
  tt$predicted_cell_type <- factor(tt$predicted_cell_type, levels = cell_types)
  
  bind_rows(
    tryCatch({kap(tt, annotated_cell_type, predicted_cell_type)}, error=function(e) NULL),
    tryCatch({bal_accuracy(tt, annotated_cell_type, predicted_cell_type)}, error=function(e) NULL),
    tryCatch({mcc(tt, annotated_cell_type, predicted_cell_type)}, error=function(e) NULL),
    tryCatch({sensitivity(tt, annotated_cell_type, predicted_cell_type)}, error=function(e) NULL),
    tryCatch({specificity(tt, annotated_cell_type, predicted_cell_type)}, error=function(e) NULL)
  )
}

### [ LOAD PREDICTIONS ] #####
scRNAseq_cl <- dir("output/v1/rare-subtype-benchmarking/", full.names = TRUE, pattern = 'scRNA')
scRNAseq_clusters <- dir("output/v1/cluster-and-interpret/", full.names = TRUE, pattern = 'Seurat')
CyTOF_cl <- dir('output/v1/rare-subtype-benchmarking/', full.names = TRUE, pattern = 'CyTOF')

scRNA_cl <- lapply(scRNAseq_cl, read_tsv) %>% 
  bind_rows() %>% 
  left_join(scRNAseq_gt)
scRNAseq_clusters <- lapply(scRNAseq_clusters, read_tsv) %>% 
  bind_rows()

CyTOF_cl <- lapply(CyTOF_cl, read_tsv) %>% 
  bind_rows() %>% 
  left_join(CyTOF_gt)

### [ CELL LEVEL ACCURACIES ] ####
cell_level_predictions <- bind_rows(scRNA_cl, CyTOF_cl)

accuracies <- cell_level_predictions %>% 
  group_by(modality, training_annotator, selection_procedure, prediction_params) %>% 
  acc_wrap()

#### NEED TO INCLUDE ALGORITHM IN PREDICTION SO THAT IT IS SAVED AND I CAN GROUP BY/FACET HERE

accuracies %>% 
  #filter(modality == "scRNASeq") %>% 
  #mutate(method = paste(modality, "-", prediction_params))
  ggplot(aes(x = selection_procedure, y = .estimate)) +
  geom_point() +
  facet_grid(prediction_params~.metric)


 accuracies %>% 
  filter(modality == "CyTOF") %>% 
  #mutate(method = paste(modality, "-", prediction_params))
  ggplot(aes(x = selection_procedure, y = .estimate)) +
  geom_boxplot() +
  facet_grid(prediction_params~.metric)

  


### [ CELL TYPE BARPLOTS ] #####
cell_level_predictions %>% 
  filter(modality == "scRNASeq") %>% 
  group_by(modality, training_annotator, selection_procedure, 
           prediction_params, predicted_cell_type) %>% 
  tally() %>% 
  ggplot(aes(x = modality, y = n, fill = predicted_cell_type)) +
  geom_bar(stat = 'identity')

cell_level_predictions %>% 
  filter(modality == "CyTOF") %>% 
  group_by(modality, training_annotator, selection_procedure, 
           prediction_params, predicted_cell_type) %>% 
  tally() %>% 
  ggplot(aes(x = modality, y = n, fill = predicted_cell_type)) +
  geom_bar(stat = 'identity')
  
### WHY ARE SOME PREDICTED CELL TYPES == NA/?