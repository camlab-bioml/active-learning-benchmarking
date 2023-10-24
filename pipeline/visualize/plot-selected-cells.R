suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  library(ComplexHeatmap)
  library(viridis)
  library(yaml)
})
source("pipeline/whatsthatcell-helpers.R")
files <- strsplit(snakemake@input[['sel_cells']], split = " ")
selected_cells <- lapply(files, function(x){
    f <- read_tsv(x)
    f$seed <- as.vector(str_match(x, "seed-[0-9]"))
    f
}) |> bind_rows()

if(snakemake@wildcards[['selection']] == 'Active-Learning_entropy' | snakemake@wildcards[['selection']] == 'Active-Learning_maxp'){
    selected_cells <- dplyr::rename(selected_cells, "params" = "method") |>
        mutate(params = gsub("seed-[0-9]", "", params),
               params = paste0(params, "-", cell_num))
}

head(selected_cells)

counts <- selected_cells |>
    group_by(params, seed, cell_type) |>
    tally() |>
    mutate(selection = snakemake@wildcards[['selection']])

# Selection specific separation of params
if(snakemake@wildcards[['selection']] == "random"){
    counts <- counts |>
        separate(params, into = c("rm_num", "cell_num", "rm_corr", "rm"), sep = "-")
}else if(snakemake@wildcards[['selection']] == "Seurat-clustering"){
    counts <- counts |>
        separate(params, into = c("knn", 'knn_val', 'res', 'res_val', 'num', 'num2', 'cell_num'), sep = '-')
}

if(snakemake@wildcards[['selection']] == "Seurat-clustering" | snakemake@wildcards[['selection']] == "random"){
    pdf(snakemake@output[['barplot']], height = 7, width = 10)
        counts |> 
            ggplot(aes(x = reorder(cell_num, as.integer(cell_num)), y = n, fill = cell_type)) + 
            geom_bar(stat = 'identity', position = 'fill') + 
            labs(x = "Number of cells") + 
            scale_fill_manual(values = cell_type_colours(snakemake@wildcards[['modality']])) +
            whatsthatcell_theme() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
            facet_wrap(~seed, nrow = 2)
    dev.off()
}else if(snakemake@wildcards[['selection']] == 'Active-Learning_entropy' | snakemake@wildcards[['selection']] == 'Active-Learning_maxp'){
    counts <- counts |>
        separate(params, c('rm_a', 'rm_l', 'rm_gt', 'initial', 'rm_strat', 'strat', 'rm_al', 'al', 'rm_rand', 'rand', 'rm_cor', 'cor', 'cell_num'), sep = '-')
    
    p <- ggplot(counts, aes(x = reorder(cell_num, as.integer(cell_num)), y = n, fill = cell_type)) + 
            geom_bar(stat = 'identity', position = 'fill') + 
            labs(x = "Number of cells") + 
            scale_fill_manual(values = cell_type_colours(snakemake@wildcards[['modality']])) +
            whatsthatcell_theme() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    
    if(snakemake@params[['cor_or_rand']] == "cor"){
        pdf(snakemake@output[['barplot']], height = 12, width = 12)
            print(p + facet_grid(initial + cor ~ seed))
        dev.off()
    }else if(snakemake@params[['cor_or_rand']] == 'rand'){
        pdf(snakemake@output[['barplot']], height = 12, width = 12)
            print(p + facet_grid(initial + rand ~ seed))
        dev.off()
    }
}