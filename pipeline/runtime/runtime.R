suppressPackageStartupMessages({
  library(tidyverse)
  library(SingleCellExperiment)
  library(yaml)
  library(data.table)
  library(caret)
  library(Seurat)
  devtools::load_all("/leader")
})
source("pipeline/whatsthatcell-helpers.R")
set.seed(1)

markers <- read_yaml(snakemake@input$markers)
sce <- readRDS(snakemake@input$sce)

df_expression <- load_scs(sce)
df_expression$cell_type <- NA
df_expression$gt_cell_type <- sce$CellType
df_expression$iteration <- NA
df_expression$corrupted <- NA

total_initial_cells <- 20

random_cell_idx <- sample(1:nrow(df_expression), total_initial_cells)
df_expression$cell_type[random_cell_idx] <- df_expression$gt_cell_type[random_cell_idx]
df_expression$iteration[random_cell_idx] <- 0

entropies <- list()

# Remove genes with 0 expression
df_expression_non_zero <- select(df_expression, -where(is.character), -iteration, -cell_type, -corrupted)
df_expression_non_zero <- df_expression_non_zero[,colSums(df_expression_non_zero) > 0]

# Calculate PCA embedding
df_PCA <- as.matrix(df_expression_non_zero) |> 
  prcomp(center = TRUE, scale. = TRUE)

df_PCA <- df_PCA$x |> 
  as.data.frame()

df_PCA <- bind_cols(
  tibble(X1 = df_expression$X1),
  df_PCA[,1:min(20, ncol(df_PCA))], 
  tibble(cell_type = df_expression$cell_type,
         gt_cell_type = df_expression$gt_cell_type,
         iteration = df_expression$iteration)
)

# Entropy 
AL_entropy_timing <- lapply(c("rf", "multinom"), function(al_alg){
  lapply(c("highest_entropy", "0.75_quant_entropy", "0.95_quant_entropy"), function(strat){
    start <- Sys.time()
    AL <- active_learning_wrapper(select(df_PCA, -gt_cell_type, -iteration), 
                                  al_alg,
                                  strat,
                                  1, 
                                  entropies, 
                                  0,
                                  'entropy')
    end <- Sys.time()
    diff <- end - start
    
    tibble(AL_alg = al_alg, method = strat, time = as.numeric(diff))
  })
}) |> bind_rows()


# Maxp
AL_maxp_timing <- lapply(c("rf", "multinom"), function(al_alg){
  lapply(c('0.05_quant_maxp', '0.25_quant_maxp', 'lowest_maxp'), function(strat){
    start <- Sys.time()
    AL <- active_learning_wrapper(select(df_PCA, -gt_cell_type, -iteration), 
                                  al_alg,
                                  strat,
                                  1, 
                                  entropies, 
                                  0,
                                  'maxp')
    end <- Sys.time()
    diff <- end - start
    
    tibble(AL_alg = al_alg, method = strat, time = as.numeric(diff))
  })
}) |> bind_rows()


# AR
seu <- CreateSeuratObject(counts = assays(sce)$logcounts)
seu <- FindVariableFeatures(seu, selection.method = 'vst', nfeatures = 2000)

all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)

# Cluster based
start <- Sys.time()
AR <- adaptive_reweighting(seu, n_requested_cells = 100, mode = "cluster_based")
end <- Sys.time()

AR_Nomarker <- tibble(AL_alg = NA, method = "AR No Marker", time = as.numeric(end - start))

# Marker based
markers <- get_markers(snakemake@input$markers)
start <- Sys.time()
AR <- adaptive_reweighting(seu, markers, n_requested_cells = 100, mode = "marker_based")
end <- Sys.time()

AR_marker <- tibble(AL_alg = NA, method = "AR Marker", time = as.numeric(end - start))


bind_rows(AL_maxp_timing, 
          AL_entropy_timing,
          AR_marker, AR_Nomarker) |> 
  mutate(method = case_when(method == "0.05_quant_maxp" ~ "AL 0.05-maxp",
                            method == "0.25_quant_maxp" ~ "AL 0.25-maxp",
                            method == "lowest_maxp" ~ "AL Lowest maxp",
                            method == "0.95_quant_entropy" ~ "AL 0.95-entropy",
                            method == "0.75_quant_entropy" ~ "AL 0.75-entropy",
                            method == "highest_entropy" ~ "AL Highest entropy",
                            TRUE ~ method),
           modality = snakemake@wildcards$modality,
           seed = snakemake@wildcards$s) |> 
  write_tsv(snakemake@output$timing)


