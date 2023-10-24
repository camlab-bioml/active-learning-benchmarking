library(tidyverse)
library(ComplexHeatmap)

cytof_sim <- read_tsv(snakemake@input$cytof_sim)
scrna_sim <- read_tsv(snakemake@input$scrna_sim)
snrna_sim <- read_tsv(snakemake@input$snrna_sim)

plot_dist_mat <- function(df){
  pivot_wider(df, names_from = "cell_type2", values_from = "cosine_similarity") |> 
    select(-c(cohort, seed)) |> 
    as.data.frame() |> 
    column_to_rownames("cell_type1") |> 
    as.matrix() |> 
    Heatmap(name = "Cosine\nsimilarity")
}

pdf(snakemake@input$scrna, height = 5, width = 5.5)
  plot_dist_mat(scrna_sim)
dev.off()

pdf(snakemake@input$snrna, height = 5, width = 5.5)
  plot_dist_mat(snrna_sim)
dev.off()

pdf(snakemake@input$cytof, height = 5, width = 5.5)
  plot_dist_mat(cytof_sim)
dev.off()
