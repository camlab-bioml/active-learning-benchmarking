suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
})
source("pipeline/whatsthatcell-helpers.R")

entropies <- lapply(snakemake@input$entropy_files, function(x){
  read_tsv(x) |> 
    filter(no_cells_annotated == 20) |> 
    mutate(rem_cell_type = gsub('AL_alg-(.*?)-strat-(.*?)-init-(.*?)-rem_celltype-',
                                "", params),
           rem_cell_type = gsub("-seed-[0-9]", "", rem_cell_type),
           AL_alg = gsub("AL_alg-", "", params),
           AL_alg = gsub("-strat-.*", "", AL_alg),
           strat = gsub("AL_alg-(.*?)-strat-", "", params),
           strat = gsub("-init-.*", "", strat),
           initial = gsub("(.*?)-init-", "", params),
           initial = gsub("-rem_celltype.*", "", initial),
           seed = gsub(".*-seed-", "", params)) |> 
    select(-c(criterion, no_cells_annotated, params, comp))
}) |> bind_rows()

png(snakemake@output$entropies, height = 600, width = 800)
  mutate(entropies, gt_rem = rem_cell_type == gt_cell_type) |> 
    ggplot(aes(x = as.factor(num_missing_cells), y = criterion_val, fill = gt_rem)) +
    geom_boxplot(outlier.colour = "grey60", outlier.size = 0.8) +
    scale_fill_manual(values = c("#ffa600", "#9763ff")) +
    labs(x = "Number of missing cells", y = "Scaled entropy/maximum probability", 
         fill = "Cell type\nremoved?", title = snakemake@wildcards$modality) +
    facet_grid(strat ~ AL_alg + initial) +
    whatsthatcell_theme()
dev.off()


mutate(entropies, gt_rem = rem_cell_type == gt_cell_type) |> 
  filter(strat == "highest_entropy" & initial == "random") |> 
  mutate(gt_rem = case_when(gt_rem == TRUE ~ "Cell type removed",
                            gt_rem == FALSE ~ "Cell type not removed")) |>
  group_by(num_missing_cells, gt_rem, AL_alg) |>
  summarize(criterion_val = median(criterion_val)) |>
  write_tsv(snakemake@output$example)
