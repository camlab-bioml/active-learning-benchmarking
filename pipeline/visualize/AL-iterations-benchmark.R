library(tidyverse)
library(yardstick)
source("pipeline/whatsthatcell-helpers.R")

### [ LOAD GROUND TRUTH ] #####
scRNAseq <- readRDS(snakemake@input[['sce']])
sce <- readRDS("data/CyTOF/CyTOF-test.rds")
expression_gt <- tibble(cell_id = colnames(sce),
                        annotated_cell_type = sce$CellType)

AL_files <- list.files("output/v3/rare-subtype-benchmarking/",
                 pattern = "CyTOF-Active-Learning", full.names = TRUE)

AL_files <- snakemake@input[['assignments']]

AL <- lapply(AL_files, read_tsv) %>% 
  bind_rows()  %>% 
  separate(prediction_params, c('m1', 'm2', 'rm_it', 'rm_na', 'rm_knn', 'rm_na2', 
                     'rm_res', 'rm_na3', 'rm_cell_num', 'cell_num', 'rm_rand', 'rand', 
                     'rm_corr', 'corrupted'), sep = '-') %>% 
  select(-starts_with('rm')) %>% 
  unite(method, c(m1, m2), sep = '-') %>% 
  mutate(selection_procedure = gsub("Active-Learning-strategy-", "", selection_procedure))


AL <- left_join(AL, expression_gt)

acc <- AL %>% 
  filter(predicted_cell_type != "unassigned") %>% 
  group_by(method, cell_num, rand, corrupted, selection_procedure) %>% 
  acc_wrap()

pdf(snakemake@output[['pdf']], width = 12, height = 4)
  acc %>% 
    mutate(cell_num = factor(as.integer(cell_num))) %>% 
    filter(method == "CyTOF-LDA") %>% 
    ggplot(aes(x = cell_num, y = .estimate, color = selection_procedure, group = selection_procedure)) +
    geom_point() +
    geom_line(aes(group = selection_procedure)) +
    scale_x_discrete(breaks = acc$cell_num) +
    ylim(0,1) +
    theme_bw() +
    facet_grid(rand + corrupted~.metric)
dev.off()
