library(tidyverse)
library(RColorBrewer)
library(patchwork)
source("pipeline/whatsthatcell-helpers.R")

files <- snakemake@input[['f']]

## Create plots

if(snakemake@wildcards[['AL_type']] == 'Active-Learning_entropy'){
  ylab = "Entropy"
}else if(snakemake@wildcards[['AL_type']] == 'Active-Learning_maxp'){
  ylab = "Maxp"
}
plots <- lapply(files, function(x){
  vals <- read_tsv(x)
  col_numbers <- vals$no_cells_annotated %>% 
    unique %>% 
    length()
  colors <- colorRampPalette(brewer.pal(9, "Blues"))(col_numbers)

  create_entropy_boxplot(vals, ylab, colors)
})


pdf(snakemake@output[['box']], height = 9, width = 14)
  plots[[1]] / plots[[2]] / plots[[3]]
dev.off()
  
