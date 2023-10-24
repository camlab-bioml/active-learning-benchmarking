suppressPackageStartupMessages({
  library(tidyverse)
})
source("pipeline/whatsthatcell-helpers.R")

similarity <- read_tsv(snakemake@input$similarity)

# files <- list.files("output/v6/results/rem_cell_type",
#                     pattern = glob2rx("Init-sel-random-rem-*-scRNASeq-highest_entropy-ALAlg-multinom-cell_num-20-seed-0.tsv"),
#                     full.names = TRUE)


entropy <- lapply(snakemake@input$entropies, read_tsv) |> 
  bind_rows() |> 
  filter(num_missing_cells == 0 | num_missing_cells == 3) |> 
  select(-criterion, -comp)
#save.image('debug-dist-rem.Rdata')
entr <- entropy |> 
  mutate(al = case_when(grepl('AL_alg-multinom', params) ~ "multinom",
                        grepl('AL_alg-rf', params) ~ 'rf'),
         strat = case_when(grepl('strat-highest_entropy', params) ~ 'highest_entropy',
                           grepl('strat-lowest_maxp', params) ~ 'lowest_maxp'),
         init = case_when(grepl('init-random', params) ~ "random",
                          grepl('init-ranking', params) ~ "ranking"),
         ct = str_extract(params, 'rem_celltype.*'),
         s = str_extract(params, '-seed-.*')) |>
  mutate(ct = gsub("rem_celltype-", "", ct),
         ct = gsub("-seed-[0-9]", "", ct),
         s = gsub("-seed-", "", s)) |>
  select(-params) |>
  pivot_wider(names_from = "num_missing_cells", values_from = "criterion_val") |> 
  mutate(entropy_diff = `0` - `3`) |> 
  select(-c(`0`, `3`))


# entr |> 
#   left_join(similarity, by = c("rem_cell_type" = "cell_type1", 
#                                "gt_cell_type" = "cell_type2")) |> 
#   ggplot(aes(x = mean_cosine_similarity, y = entropy_diff, color = gt_cell_type)) +
#   geom_point() +
#   labs(x = "Mean cosine similarity across seeds from removed cell type",
#        y = "Entropy difference\n(entropy w/ 0 cells - entropy w/ 3 cells)",
#        title = "Faceted by removed cell type",
#        color = "Ground truth cell type") +
#   scale_color_manual(values = cell_type_colours("scRNASeq")) +
#   facet_wrap(~rem_cell_type, nrow = 1) +
#   whatsthatcell_theme() +
#   theme(legend.position = "bottom")


closest_cell_type <- similarity |> 
  filter(cell_type1 != cell_type2) |>
  group_by(cell_type1) |> 
  arrange(-mean_cosine_similarity) |> 
  slice_head(n = 1)

pdf(snakemake@output$pdf, height = 6, width = 8)
  entr |> 
    left_join(closest_cell_type, by = c("gt_cell_type" = "cell_type1")) |> 
    group_by(gt_cell_type) |> 
    mutate(mean_entropy_diff = mean(na.omit(entropy_diff))) |> 
    ungroup() |> 
    ggplot(aes(x = mean_cosine_similarity, y = mean_entropy_diff, color = gt_cell_type)) +
    geom_point() +
    labs(x = "Mean cosine similarity across seeds from closest cell type",
         y = "Mean entropy difference\n(entropy w/ 0 cells - entropy w/ 3 cells)",
         color = "Ground truth cell type") +
    scale_color_manual(values = cell_type_colours(snakemake@wildcards$modality)) +
    whatsthatcell_theme()
dev.off()

