suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(patchwork)
})
source("pipeline/whatsthatcell-helpers.R")

# Read in imbalanced accuracies
sc_acc <- read_tsv(snakemake@input$scrna_acc)
sn_acc <- read_tsv(snakemake@input$snrna_acc)
cy_acc <- read_tsv(snakemake@input$cytof_acc)
lung_acc <- read_tsv(snakemake@input$scrna_lung_acc)
liver_acc <- read_tsv(snakemake@input$liver_acc)
vasc_acc <- read_tsv(snakemake@input$vasc_acc)

# Calculate improvement
imbalance_score <- bind_rows(sc_acc, sn_acc, cy_acc, lung_acc, liver_acc, vasc_acc) |> 
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
    facet_wrap(~ method, nrow = 1) +
    whatsthatcell_theme() +
    theme(legend.position = "bottom")
dev.off()



plot_supplementary_imb <- function(imb_score, sel_cohort, sel_method,
                                   plot_title = NULL, axis_title = FALSE){
  p <- filter(imb_score, modality == sel_cohort, method == sel_method) |> 
    ggplot(aes(x = comp, y = diff, fill = strat)) +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    scale_fill_manual(values = sel_met_cols) +
    labs(x = "Cell type similarity", fill = "Selection procedure") +
    facet_grid(method ~ al + init, scales = "free_y") +
    whatsthatcell_theme() +
    theme(axis.title.y = element_blank())
  
  if(!is.null(plot_title)){
    p <- p + 
      labs(title = plot_title)
  }
  
  if(!axis_title){
    p <- p + 
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank())
  }else{
    p <- p + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  }
  
  p
}

scrna_supp <- plot_supplementary_imb(comb_imb_score, "scRNASeq", "Random-Forest", 
                                     "scRNASeq - Breast cancer cell lines") /
  plot_supplementary_imb(comb_imb_score, "scRNASeq", "SCN-labels") /
  plot_supplementary_imb(comb_imb_score, "scRNASeq", "SVM-rejection") /
  plot_supplementary_imb(comb_imb_score, "scRNASeq", "scmap-clusters") /
  plot_supplementary_imb(comb_imb_score, "scRNASeq", "scmap-sc") /
  plot_supplementary_imb(comb_imb_score, "scRNASeq", "singleR-labels", NULL, TRUE)
scrna_lung_supp <- plot_supplementary_imb(comb_imb_score, "scRNALung", "Random-Forest", 
                                     "scRNASeq - Lung cancer cell lines") /
  plot_supplementary_imb(comb_imb_score, "scRNALung", "SCN-labels") /
  plot_supplementary_imb(comb_imb_score, "scRNALung", "SVM-rejection") /
  plot_supplementary_imb(comb_imb_score, "scRNALung", "scmap-clusters") /
  plot_supplementary_imb(comb_imb_score, "scRNALung", "scmap-sc") /
  plot_supplementary_imb(comb_imb_score, "scRNALung", "singleR-labels", NULL, TRUE)


pdf(snakemake@output$supp_1, height = 13, width = 12)
  (scrna_supp | scrna_lung_supp) + 
    plot_layout(guides = "collect")
dev.off()


vasc_supp <- plot_supplementary_imb(comb_imb_score, "tabulaVasc", "Random-Forest", 
                                    "scRNASeq - Vasculature") /
  plot_supplementary_imb(comb_imb_score, "tabulaVasc", "SCN-labels") /
  plot_supplementary_imb(comb_imb_score, "tabulaVasc", "SVM-rejection") /
  plot_supplementary_imb(comb_imb_score, "tabulaVasc", "scmap-clusters") /
  plot_supplementary_imb(comb_imb_score, "tabulaVasc", "scmap-sc") /
  plot_supplementary_imb(comb_imb_score, "tabulaVasc", "singleR-labels", NULL, TRUE)

liver_supp <- plot_supplementary_imb(comb_imb_score, "liverAtlas", "Random-Forest", 
                                    "scRNASeq - Liver") /
  plot_supplementary_imb(comb_imb_score, "liverAtlas", "SCN-labels") /
  plot_supplementary_imb(comb_imb_score, "liverAtlas", "SVM-rejection") /
  plot_supplementary_imb(comb_imb_score, "liverAtlas", "scmap-clusters") /
  plot_supplementary_imb(comb_imb_score, "liverAtlas", "scmap-sc") /
  plot_supplementary_imb(comb_imb_score, "liverAtlas", "singleR-labels", NULL, TRUE)

pdf(snakemake@output$supp_2, height = 13, width = 12)
  (vasc_supp | liver_supp) + 
    plot_layout(guides = "collect")
dev.off()


cytof_supp <- plot_supplementary_imb(comb_imb_score, "CyTOF", "Random-Forest", 
                                     "CyTOF - Bone marrow") /
  plot_supplementary_imb(comb_imb_score, "CyTOF", "CyTOF-LDA", NULL, TRUE)

scrna_lung_supp <- plot_supplementary_imb(comb_imb_score, "scRNALung", "Random-Forest", 
                                          "scRNALung - Lung cancer cell lines") /
  plot_supplementary_imb(comb_imb_score, "scRNALung", "SCN-labels") /
  plot_supplementary_imb(comb_imb_score, "scRNALung", "SVM-rejection") /
  plot_supplementary_imb(comb_imb_score, "scRNALung", "scmap-clusters") /
  plot_supplementary_imb(comb_imb_score, "scRNALung", "scmap-sc") /
  plot_supplementary_imb(comb_imb_score, "scRNALung", "singleR-labels", NULL, TRUE)

pdf(snakemake@output$supp_3, height = 16, width = 8)
  (cytof_supp / scrna_lung_supp) +
    plot_layout(guides = "collect", nrow = 3, ncol = 1, heights = c(0.3, 0.3, 3))
dev.off()

  