suppressPackageStartupMessages({
  library(tidyverse)
})
source("pipeline/whatsthatcell-helpers.R")

uncertainties <- lapply(snakemake@input$tsv, read_tsv) |> 
  bind_rows() |> 
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
  select(-params)

cell_selections <- uncertainties |> 
  group_by(ct, no_cells_annotated, num_missing_cells, s)

if(snakemake@wildcards$strat == 'highest_entropy'){
  cell_selections <- cell_selections |>
    arrange(-criterion_val) |> 
    slice_head(n = 10)
  y_lab <- "Entropy"
}else if(snakemake@wildcards$strat == 'lowest_maxp'){
  cell_selections <- cell_selections |>
    arrange(criterion_val) |> 
    slice_head(n = 10)
  y_lab <- "Maxp"
}

old_selected_cells_plot <- cell_selections |> 
    mutate(sel_of_rem_type = ct == gt_cell_type) |> 
    group_by(ct, num_missing_cells, no_cells_annotated, s) |> 
    summarise(prop = sum(sel_of_rem_type) / 10) |> 
    ungroup() |> 
    group_by(ct, num_missing_cells, no_cells_annotated) |> 
    summarise(mean_prop = mean(prop), sd_prop = sd(prop)) |> 
    ggplot(aes(x = no_cells_annotated,
                colour = as.character(num_missing_cells), group = num_missing_cells)) +
    geom_ribbon(aes(ymin = mean_prop - sd_prop, ymax = mean_prop + sd_prop, 
                    fill = as.character(num_missing_cells)), alpha = .2, linetype='blank') +
    geom_line(aes(y = mean_prop), size = 2) +
    geom_point(aes(y = mean_prop), size = 4) +
    scale_x_continuous(breaks = c(20, 30, 40)) +
    labs(x = "Number of cells annotated", 
        y = "Mean proportion of cells assigned to the removed type",
        color = "Number of cells\nin training set",
        fill = "Number of cells\nin training set") +
    facet_wrap(~ct) +
    whatsthatcell_theme()

pdf(snakemake@output$plot, height = 4, width = 35)
  uncertainties |> 
    filter(no_cells_annotated == as.integer(snakemake@wildcards$num)) |> 
    ggplot(aes(x = gt_cell_type, y = criterion_val, fill = as.character(num_missing_cells))) + 
    geom_boxplot() + 
    labs(y = y_lab, fill = "Number of cells of\nthe type removed\nin ranked dataset",
        title = paste("Facet = cell type removed"), 
        x = "Ground truth cell type") +
    facet_wrap(~ct, nrow = 1) + 
    whatsthatcell_theme() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()
