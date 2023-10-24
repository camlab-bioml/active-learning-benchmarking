library(tidyverse)

training <- read_tsv(snakemake@input[['assignment']]) %>% 
  filter(!is.na(iteration)) |>
  arrange(iteration)

subset <- training[1:as.integer(snakemake@wildcards[['subset_val']]),]
subset$cell_num <- snakemake@wildcards[['subset_val']]

write_tsv(subset, snakemake@output[['split']])