suppressPackageStartupMessages({
  library(yaml)
  library(SingleCellExperiment)
  devtools::load_all('/leader')
  library(Seurat)
})

source("pipeline/whatsthatcell-helpers.R")
set.seed(1)

sce <- readRDS(snakemake@input$sce)
markers <- get_markers(snakemake@input$markers)

sce_gt <- tibble(cell_id = colnames(sce),
                 cell_type = sce$CellType)

seu <- CreateSeuratObject(counts = assays(sce)$logcounts)
seu <- FindVariableFeatures(seu, selection.method = 'vst', nfeatures = 2000)
all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)

AR <- adaptive_reweighting(seu, markers$cell_types, n_requested_cells = 100, mode = "marker_based")

select(AR, cell_id) |> 
    left_join(sce_gt, by = "cell_id") |>
    mutate(params = "Default") |> 
    write_csv(snakemake@output$selected)