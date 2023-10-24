library(tidyverse)
library(SingleCellExperiment)
library(scater)
library(annotables)
library(patchwork)
library(Matrix)
set.seed(1)

sces <- lapply(snakemake@input$mat, readRDS)
sce <- do.call('cbind', sces)

sce$CellType <- str_split_fixed(colnames(sce), "_", 2)[,1]

## More QC
grch38 <- select(grch38, symbol, chr) |> 
    filter(!grepl("GL|KI|CHR_", chr)) |>
    filter(symbol %in% rownames(sce)) |>
    distinct() |> 
    group_by(symbol) |>
    slice_head(n = 1)

gene_metadata <- tibble(symbol = rownames(sce)) |>
    left_join(grch38)

if(!all(gene_metadata$symbol == rownames(sce))){
    stop("Error: gene order is not the same")
}

rowData(sce)$chr <- gene_metadata$chr

sce <- logNormCounts(sce)

head(rownames(sce))
mt_genes <- which(rowData(sce)$chr == "MT")
ribo_genes <- grepl("^RP[LS]", rownames(sce))
feature_ctrls <- list(mito = rownames(sce)[mt_genes],
                      ribo = rownames(sce)[ribo_genes])
lapply(feature_ctrls, head)

sce <- addPerCellQC(sce, subsets = feature_ctrls)

pdf(snakemake@output$detected)
  plotColData(sce, y = "subsets_mito_percent")
dev.off()


# run PCA and UMAP
sce <- runPCA(sce, n_dimred = 20)
sce <- runUMAP(sce, n_dimred = 20, ncomponents = 2)

dim_reductions <- (plotReducedDim(sce, dimred = 'UMAP', colour_by = "detected") | 
  plotReducedDim(sce, dimred = 'UMAP', colour_by = 'CellType')) /
  (plotReducedDim(sce, dimred = 'PCA', colour_by = 'detected') |
  plotReducedDim(sce, dimred = 'PCA', colour_by = 'CellType'))

png(snakemake@output$dim_red, width = 1500, height = 1500, res = 150)
  dim_reductions
dev.off()

saveRDS(sce, snakemake@output$sce)