library(DropletUtils)
library(Matrix)
library(tidyverse)
library(scater)
library(annotables)
library(SingleCellExperiment)

matrix <- readMM(snakemake@input$mat)
barcodes <- read_tsv(snakemake@input$barcodes, col_names = FALSE)
features <- read_tsv(snakemake@input$features, col_names = FALSE)

rownames(matrix) <- features$X1
colnames(matrix) <- barcodes$X1

grch38_filtered <- grch38 |> 
  filter(!grepl("GL|KI|CHR_", chr) & biotype == "protein_coding") |>
  filter(!grepl("^RP[L|S]|^MALAT1|^HSP|^FOS|^JUN", symbol))
gene_map <- select(grch38_filtered, ensgene, symbol) |> deframe()

# Remove any genes that don't have a hugo name
matrix <- matrix[rownames(matrix) %in% names(gene_map), ]

# change rownames to hugo
rownames(matrix) <- gene_map[rownames(matrix)]

sel_cells <- grepl(snakemake@wildcards$cell_line, colnames(matrix))
subset_mat <- matrix[, sel_cells]

subset_mat <- as.matrix(subset_mat) |>
  as.data.frame() |>
  rownames_to_column("gene") |> 
  pivot_longer(-gene, names_to = "cell_id", values_to = "counts")

subset_mat <- subset_mat |>
  group_by(gene, cell_id) %>% 
  summarise(counts = sum(counts)) |> 
  ungroup()

subset_mat <- pivot_wider(subset_mat, values_from = "counts", names_from = "cell_id") |>
  as.data.frame() |>
  column_to_rownames("gene") |>
  t()

subset_mat <- as(subset_mat, "sparseMatrix")

sce <- SingleCellExperiment(assays = list(counts = t(subset_mat)))

saveRDS(sce, snakemake@output$sce)

