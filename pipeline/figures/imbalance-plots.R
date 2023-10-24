library(tidyverse)
library(scales)
library(patchwork)

sc_acc <- read_tsv(snakemake@input$scrna_acc)
sn_acc <- read_tsv(snakemake@input$snrna_acc)
cy_acc <- read_tsv(snakemake@input$cytof_acc)

source("pipeline/whatsthatcell-helpers.R")


imbalance_score <- bind_rows(sc_acc, sn_acc, cy_acc) |> 
  dplyr::filter(.metric == "bal_accuracy") |> 
  pivot_wider(names_from = "similarity", values_from = ".estimate") |> 
  mutate(similar = (`imbalanced-similar` - `balanced-similar`) / `balanced-similar`,
         different =  (`imbalanced-different` - `balanced-different`) / `balanced-different`) |> 
  filter(!is.na(similar)) |> 
  filter(!is.na(different)) |> 
  pivot_longer(c(similar, different), names_to = "comp", values_to = "diff")
  
al_imbalance <- filter(imbalance_score, !grepl("Seurat|random", strat))
AR_rand_imbalance <- filter(imbalance_score, grepl("Seurat|random", strat))

comb_imb_score <- bind_rows(al_imbalance,
          mutate(AR_rand_imbalance, al = "rf", init = "random"),
          mutate(AR_rand_imbalance, al = "rf", init = "ranking"),
          mutate(AR_rand_imbalance, al = "multinom", init = "random"),
          mutate(AR_rand_imbalance, al = "multinom", init = "ranking")) |> 
  mutate(strat = case_when(strat == "highest_entropy" ~ "AL Highest entropy",
                           strat == "0.75_quant_entropy" ~ "AL 0.75-entropy",
                           strat == "0.95_quant_entropy" ~ "AL 0.95-entropy",
                           strat == "lowest_maxp" ~ "AL Lowest maxp",
                           strat == "0.05_quant_maxp" ~ "AL 0.05-maxp",
                           strat == "0.25_quant_maxp" ~ "AL 0.25-maxp",
                           strat == "MarkerSeurat_clustering" ~ "AR Marker",
                           strat == "NoMarkerSeurat_clustering" ~ "AR No Marker",
                           strat == "random" ~ "Random"),
         al = case_when(al == "rf" ~ "RF",
                        al == "multinom" ~ "LR"))
  
scrnaseq <- filter(comb_imb_score, modality == "scRNASeq") |> 
  ggplot(aes(x = comp, y = diff, fill = strat)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = sel_met_cols) +
  labs(x = "Cell type similarity", 
       title = "scRNASeq", fill = "Selection procedure") +
  facet_wrap(~ al + init + method, scales = "free_y") +
  whatsthatcell_theme() +
  theme(axis.title.y = element_blank())

snrnaseq <- filter(comb_imb_score, modality == "snRNASeq") |> 
  ggplot(aes(x = comp, y = diff, fill = strat)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  scale_y_continuous(breaks = c(-0.5, 0, 0.5)) +
  scale_fill_manual(values = sel_met_cols) +
  labs(x = "Cell type similarity", 
       y = "Change in balanced accuracy\n(imbalanced - balanced) / balanced",
       title = "snRNASeq", fill = "Selection procedure") +
  facet_wrap(~ al + init + method, scales = "free_y", ncol = 3) +
  whatsthatcell_theme()

cytof <- filter(comb_imb_score, modality == "CyTOF") |> 
  ggplot(aes(x = comp, y = diff, fill = strat)) +
  geom_boxplot() +
  geom_hline(yintercept = 0) +
  scale_y_continuous(breaks = c(-0.4, -0.2, 0)) +
  scale_fill_manual(values = sel_met_cols) +
  labs(x = "Cell type similarity", 
       y = "Change in balanced accuracy\n(imbalanced - balanced) / balanced",
       title = "CyTOF", fill = "Selection procedure") +
  facet_wrap(~ al + init + method, scales = "free_y", ncol = 4) +
  whatsthatcell_theme() +
  theme(axis.title.y = element_blank())

plotting_order <- filter(comb_imb_score, modality == "snRNASeq" & al == "RF" & init == "ranking") |> 
  group_by(strat) |> 
  summarize(mean_estimate = mean(diff)) |> 
  arrange(mean_estimate) |> 
  pull(strat)

pdf(snakemake@output$main_fig, height = 3, width = 7)
  filter(comb_imb_score, modality == "snRNASeq" & al == "RF" & init == "ranking") |> 
    mutate(strat = factor(strat, levels = plotting_order)) |> 
    ggplot(aes(x = comp, y = diff, fill = strat)) +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    scale_y_continuous(breaks = c(-0.5, 0, 0.5)) +
    scale_fill_manual(values = sel_met_cols) +
    labs(x = "Cell type similarity", 
         y = "Change in balanced accuracy\n(imbalanced - balanced) / balanced",
         fill = "Selection procedure") +
    facet_wrap(~ method) +
    whatsthatcell_theme()
dev.off()


pdf(snakemake@output$sup_fig, height = 18, width = 9)
  (scrnaseq / snrnaseq / cytof) + 
    plot_annotation(tag_levels = "A") + 
  plot_layout(guides = "collect", heights = c(2, 2, 1))
dev.off()

