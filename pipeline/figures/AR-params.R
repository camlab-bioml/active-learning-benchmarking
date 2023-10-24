suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
source("pipeline/whatsthatcell-helpers.R")


cytof_acc <- read_tsv(snakemake@input$cytof_acc) |> 
  mutate(cohort = "CyTOF")
scrna_acc <- read_tsv(snakemake@input$scrna_acc) |> 
  mutate(cohort = "scRNASeq")
snrna_acc <- read_tsv(snakemake@input$snrna_acc) |> 
  mutate(cohort = "snRNASeq")
scrna_lung_acc <- read_tsv(snakemake@input$scrna_lung) |> 
  mutate(cohort = "scRNALung")
liverAtlas_acc <- read_tsv(snakemake@input$liver) |> 
  mutate(cohort = "liverAtlas")
tabulaVasc_acc <- read_tsv(snakemake@input$vasc) |> 
  mutate(cohort = "tabulaVasc")

acc <- bind_rows(cytof_acc, scrna_acc, snrna_acc,
                 scrna_lung_acc, liverAtlas_acc, tabulaVasc_acc)



plot_knn_res <- function(acc, fill, cohort){
  acc <- filter(acc, selection_procedure == "MarkerSeurat_clustering" |
                  selection_procedure == "NoMarkerSeurat_clustering") |> 
    mutate(selection_procedure = case_when(selection_procedure == "MarkerSeurat_clustering" ~ "AR Marker",
                                           selection_procedure == "NoMarkerSeurat_clustering" ~ "AR No Marker")) |> 
    filter(.metric == "f_meas")
  
  if(fill == "res"){
    acc |>
      ggplot(aes(x = as.character(cell_num), y = .estimate, fill = as.character(res))) +
      geom_boxplot() +
      labs(x = "Number of cells", y = "F1 score", fill = "Clustering\nresolution",
           title = cohort) +
      facet_grid(method ~ selection_procedure + knn) +
      whatsthatcell_theme()
  }else if(fill == "knn"){
    acc |> 
      ggplot(aes(x = as.character(cell_num), y = .estimate, fill = as.character(knn))) +
      geom_boxplot() +
      labs(x = "Number of cells", y = "F1 score", fill = "Number of\nnearest\nneighbours",
           title = cohort) +
      facet_grid(method ~ selection_procedure + res) +
      whatsthatcell_theme()
  }
}


pdf(snakemake@output$res, 
    height = 23, width = 15)
  (plot_knn_res(cytof_acc, "res", "CyTOF - Bone marrow") /
    plot_knn_res(scrna_acc, "res", "scRNASeq - Breast cancer cell lines") /
    plot_knn_res(snrna_acc, "res", "snRNASeq - Pancreas cancer")) |
  (plot_knn_res(scrna_lung_acc, "res", "scRNASeq - Lung cancer cell lines") /
     plot_knn_res(liverAtlas_acc, "res", "scRNASeq - Liver") /
     plot_knn_res(tabulaVasc_acc, "res", "scRNASeq - Vasculature")) +
  plot_layout(guides = "collect")
dev.off()


pdf(snakemake@output$knn, 
    height = 23, width = 15)
  (plot_knn_res(cytof_acc, "knn", "CyTOF - Bone marrow") /
      plot_knn_res(scrna_acc, "knn", "scRNASeq - Breast cancer cell lines") /
      plot_knn_res(snrna_acc, "knn", "snRNASeq - Pancreas cancer")) |
    (plot_knn_res(scrna_lung_acc, "knn", "scRNASeq - Lung cancer cell lines") /
       plot_knn_res(liverAtlas_acc, "knn", "scRNASeq - Liver") /
       plot_knn_res(tabulaVasc_acc, "knn", "scRNASeq - Vasculature")) +
    plot_layout(guides = "collect")
dev.off()



