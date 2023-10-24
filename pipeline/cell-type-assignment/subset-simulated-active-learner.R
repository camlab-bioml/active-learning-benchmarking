library(tidyverse)


training <- read_tsv(snakemake@input[['assignment']]) %>% 
  arrange(iteration)

ranked <- filter(training, iteration == 0)

subset <- as.integer(snakemake@wildcards[['subset_val']])
if(subset != 0){
  AL <- filter(training, iteration != 0) %>% 
    filter(iteration <= subset)

  df <- bind_rows(ranked, AL)
}else{
  df <- ranked
}



write_tsv(df, snakemake@output[['split']])