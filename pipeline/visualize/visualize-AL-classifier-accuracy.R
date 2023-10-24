suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
source("pipeline/whatsthatcell-helpers.R")

maxp_acc <- lapply(snakemake@input$maxp_accs, read_tsv) |> 
  bind_rows()
entr_acc <- lapply(snakemake@input$entr_accs, read_tsv) |> 
  bind_rows()

plot_f1 <- function(acc, sel_cohort, cohort_title){
  filter(acc, cohort == sel_cohort) |> 
    separate(params, c("rm_init", "init", "rm_strat", "strat", "rm_al", "AL", "rm_rand",
                       "rand", "rm_cor", "corr", "rm_s", "seed"), "-") |> 
    select(-starts_with("rm")) |> 
    mutate(seed_cor = paste0(seed, corr, sep = "_")) |> 
    group_by(iteration, corr, init, AL, strat) |> 
    mutate(mean_f1 = mean(f1), sd_f1 = sd(f1), 
           lower_f1 = mean_f1 - sd_f1, high_f1 = mean_f1 + sd_f1) |> 
    ggplot(aes(x = iteration, color = as.character(corr), group = seed_cor)) +
    geom_ribbon(aes(ymin = lower_f1, ymax = high_f1, 
                    fill = as.character(corr)), 
                alpha = 0.05, linetype = 'blank') +
    geom_point(aes(y = mean_f1)) +
    geom_line(aes(y = mean_f1)) +
    scale_color_manual(values = c("#ffa600", "#9763ff")) +
    scale_fill_manual(values = c("#ffa600", "#9763ff")) +
    labs(x = "Training iteration", y = "Mean f1 score", 
         color = "Proportion of\nlabels corrupted",
         fill = "Proportion of\nlabels corrupted",
         title = cohort_title) +
    facet_grid(init + AL ~ strat) +
    whatsthatcell_theme()
}

cytof <- bind_rows(maxp_acc, entr_acc) |>
  plot_f1("CyTOF", "CyTOF - Bone marrow")
scrnaseq <- bind_rows(maxp_acc, entr_acc) |>
  plot_f1("scRNASeq", "scRNASeq - Breast cancer cell lines")
snrnaseq <- bind_rows(maxp_acc, entr_acc) |>
  plot_f1("snRNASeq", "snRNASeq - Pancreas cancer")
scrnalung <- bind_rows(maxp_acc, entr_acc) |>
  plot_f1("scRNALung", "scRNASeq - Lung cancer cell lines")
liverAtlas <- bind_rows(maxp_acc, entr_acc) |>
  plot_f1("liverAtlas", "scRNASeq - Liver")
tabulaVasc <- bind_rows(maxp_acc, entr_acc) |>
  plot_f1("tabulaVasc", "scRNASeq - Vasculature")

pdf(snakemake@output$pdf1, width = 11, height = 14)
  cytof / scrnaseq / snrnaseq + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A')
dev.off()

pdf(snakemake@output$pdf2, width = 11, height = 14)
  scrnalung / liverAtlas / tabulaVasc + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A')
dev.off()

