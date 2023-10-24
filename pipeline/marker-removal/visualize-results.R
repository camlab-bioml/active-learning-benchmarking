## benchmark marker corruption
suppressPackageStartupMessages({
  library(tidyverse)
})
source("pipeline/whatsthatcell-helpers.R")

acc <- lapply(snakemake@input$accs, function(x){
  read_tsv(x) |> 
    mutate(cohort = case_when(grepl("CyTOF", basename(x)) ~ "CyTOF",
                              grepl("liverAtlas", basename(x)) ~ "liverAtlas",
                              grepl("scRNALung", basename(x)) ~ "scRNALung",
                              grepl("scRNASeq", basename(x)) ~ "scRNASeq",
                              grepl("snRNASeq", basename(x)) ~ "snRNASeq",
                              grepl("tabulaVasc", basename(x)) ~ "tabulaVasc"))
}) |> bind_rows()

acc <- filter(acc, .metric == "sensitivity")

median_acc <- acc |> 
  group_by(rem_percentage,method, cohort, selection_procedure) |> 
  summarize(median_acc = median(na.omit(.estimate)))

average_method_line <- group_by(median_acc, rem_percentage, cohort) |> 
  summarize(mean_line = median(median_acc))

pdf(snakemake@output$fig, height = 5, width = 8)
  ggplot(NULL, aes(x = rem_percentage)) +
    geom_point(aes(y = median_acc, colour = method, alpha = 0.5), 
               data = median_acc) +
    geom_line(aes(y = median_acc, colour = method, alpha = 0.5), 
              data = median_acc) +
    geom_smooth(aes(y = median_acc), data = median_acc) +
    labs(x = "Percentage of markers corrupted", y = "Median sensitivity", 
         colour = "Cell type prediction method") +
    facet_wrap(~ cohort, scales = "free_y") +
    whatsthatcell_theme() +
    guides(alpha = 'none')
dev.off()

