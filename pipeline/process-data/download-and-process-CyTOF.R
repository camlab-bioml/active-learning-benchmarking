suppressPackageStartupMessages({
  library(HDCytoData)
  library(tidyverse)
  library(SingleCellExperiment)
  library(caret)
})

# Download dataset
CyTOF <- Samusik_all_SE()

# remove non-celltype markers & unassigned cells
assigned <- CyTOF[rowData(CyTOF)$population_id != "unassigned", 
                  colData(CyTOF)$marker_class == "type"]

# Convert to single cell experiment & normalize
expression <- assay(assigned)
cell_types <- as.character(rowData(assigned)$population_id)

# Normalize
expression <- asinh(expression / 5)

sce <- SingleCellExperiment(list(logcounts = t(expression)),
                            colData = DataFrame(cell_type = cell_types))
sce$cell_type <- cell_types
colnames(sce) <- paste0("Levine_32_", seq(1:ncol(sce)))

subset_types <- c("CLP", "MPP", "CMP", "GMP", "B-cell Frac A-C (pro-B cells)",
                  "IgM- IgD- B-cells", "IgD- IgMpos B cells")

sce <- sce[, sce$cell_type %in% subset_types]

set.seed(42)
subset_sce <- sce[, createDataPartition(sce$cell_type, p = 0.1)$Resample1]

# cell_types_remove <- subset_sce$cell_type %>% 
#   table() %>% 
#   as_tibble %>% 
#   filter(n < 15) %>% 
#   pull(".")
# 
# subset_sce <- subset_sce[ , !(subset_sce$cell_type %in% cell_types_remove)]

saveRDS(subset_sce, snakemake@output[['Levine_CyTOF']])
