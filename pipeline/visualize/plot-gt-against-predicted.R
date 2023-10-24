suppressPackageStartupMessages({
  library(tidyverse)
  library(ggalluvial)
  library(SingleCellExperiment)
})
source("pipeline/whatsthatcell-helpers.R")

sce <- readRDS(snakemake@input$sce)
if(is.null(sce$CellType)){
  sce$CellType <- sce$cell_type
}
if(snakemake@wildcards[['modality']] == "snRNASeq"){
  colnames(sce) <- gsub("-", "_", colnames(sce))
}
gt <- tibble(cell_id = colnames(sce),
             ground_truth = sce$CellType)

predictions <- lapply(snakemake@input$predicted, function(x){
    df <- read_tsv(x) |> 
        separate(prediction_params, c('m1', 'm2', 'rm_knn', 'knn',
                                    'rm_res', 'res', 'rm_cell', 'cell_num', 'rm_sel',
                                    'random', 'rm_cor', 'corrupted', 'rm_init', 'initial', 'rm_s', 'seed'), sep = '-') %>% 
        select(-starts_with('rm')) %>% 
        unite('method', c(m1, m2), sep = '-')
    df
}) |> bind_rows()

predictions <- left_join(predictions, gt, by = "cell_id")
predictions$cell_num <- factor(predictions$cell_num, levels = c("100", "250", "500"))

pdf(snakemake@output$pdf, width = 13, height = 15)
predictions |>
    pivot_longer(c(predicted_cell_type, ground_truth), names_to = 'predicted_or_gt', values_to = 'cell_type') |> 
    ggplot(aes(x = predicted_or_gt, stratum = cell_type, fill = cell_type, alluvium = cell_id)) +
    geom_stratum() +
    geom_flow() +
    scale_fill_manual(values = cell_type_colours(snakemake@wildcards$modality)) +
    whatsthatcell_theme() +
    facet_grid(cell_num + method ~ seed) +
    labs(x = "", y = "cell") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()