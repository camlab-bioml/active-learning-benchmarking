# Imbalance2 result summary
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})
source("pipeline/whatsthatcell-helpers.R")

f <- list.files("output/v8/results/imbalance2/acc", full.names = TRUE)

acc <- lapply(f, read_tsv) |> 
  bind_rows() |> 
  select(-c(rand, corrupted, .estimator)) |> 
  pivot_wider(names_from = balance, values_from = .estimate) |> 
  mutate(improvement =  (imbalanced - balanced) / balanced,
         strat = case_when(strat == "highest_entropy" ~ "AL Highest entropy",
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

plot_acc <- function(acc, metric, x_axis = FALSE){
  p <- acc |> 
    filter(.metric == metric) |> 
    mutate(.metric = case_when(.metric == "bal_accuracy" ~ "Balanced accuracy",
                               .metric == "f_meas" ~ "F1-score",
                               .metric == "kap" ~ "Kappa",
                               .metric == "mcc" ~ "Matthews correlation coefficient",
                               .metric == "sensitivity" ~ "Sensitivity")) |> 
    ggplot(aes(x = method, y = improvement, fill = strat)) +
    geom_hline(yintercept = 0, color = "black") +
    geom_boxplot() +
    labs(title = metric) +
    scale_fill_manual(values = sel_met_cols) +
    facet_wrap(~modality, scales = "free", nrow = 1) +
    whatsthatcell_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  if(x_axis){
    p <- p + labs(x = "Cell type assignment method")
  }else{
    p <- p + theme(axis.text.x = element_blank(),
                   axis.ticks = element_blank(),
                   axis.title.x = element_blank())
  }
  
  p
}

plotting_order <- acc |> 
  filter(improvement != Inf) |> 
  group_by(strat) |> 
  summarize(mean_estimate = mean(na.omit(improvement))) |> 
  arrange(mean_estimate) |> 
  pull(strat)

acc <- mutate(acc, strat = factor(strat, levels = plotting_order))


acc1 <- filter(acc, modality %in% c("CyTOF", "scRNASeq", "snRNASeq"))
pdf("output/v8/paper-figures/imbalance2-fig1.pdf", height = 12, width = 14)
  (plot_acc(acc1, "bal_accuracy") /
      plot_acc(acc1, "f_meas") /
      plot_acc(acc1, "kap") /
      plot_acc(acc1, "mcc") /
      plot_acc(acc1, "sensitivity", TRUE)) +
    plot_layout(guides = 'collect') &
    labs(fill = "Selection procedure")
dev.off()


acc2 <- filter(acc, modality %in% c("scRNALung", "tabulaLiver", "tabulaVasc"))
pdf("output/v8/paper-figures/imbalance2-fig2.pdf", height = 12, width = 14)
  (plot_acc(acc2, "bal_accuracy") /
      plot_acc(acc2, "f_meas") /
      plot_acc(acc2, "kap") /
      plot_acc(acc2, "mcc") /
      plot_acc(acc2, "sensitivity", TRUE)) +
    plot_layout(guides = 'collect') &
    labs(fill = "Selection procedure")
dev.off()


### Number of times at top

# average over seed, al params (ranking/random & RF/LR), AR params
mean_acc <- acc |> 
  filter(improvement != Inf) |> 
  group_by(method, strat, modality, .metric) |> 
  summarize(mean_improvement = mean(na.omit(improvement))) |> 
  arrange(-mean_improvement) |> 
  ungroup() |> 
  mutate(.metric = case_when(.metric == "bal_accuracy" ~ "Balanced accuracy",
                             .metric == "f_meas" ~ "F1-score",
                             .metric == "kap" ~ "Kappa",
                             .metric == "mcc" ~ "Matthews correlation coefficient",
                             .metric == "sensitivity" ~ "Sensitivity"))

# out of all methods & modalities, how often is each metric number 1?
top_1 <- group_by(mean_acc, method, modality, .metric) |> 
  slice_head(n = 1) |> 
  ungroup() |> 
  group_by(strat, modality, .metric) |> 
  tally() |> 
  ungroup() |> 
  pivot_wider(names_from = .metric, values_from = n, values_fill = 0) |> 
  pivot_longer(-c(strat, modality), names_to = ".metric", values_to = "n")

sel_procedures <- unique(mean_acc$strat)
top1_missing <- sel_procedures[!sel_procedures %in% top_1$strat]

top1_missing <- expand.grid(top1_missing, unique(top_1$modality), unique(top_1$.metric))
colnames(top1_missing) <- c("strat", "modality", ".metric")
top1_missing$n <- 0

top_1 <- bind_rows(top_1, top1_missing)
  
top_3 <- group_by(mean_acc, method, modality, .metric) |> 
  slice_head(n = 3) |> 
  ungroup() |> 
  group_by(strat, modality, .metric) |> 
  tally() |> 
  ungroup() |> 
  pivot_wider(names_from = .metric, values_from = n, values_fill = 0) |> 
  pivot_longer(-c(strat, modality), names_to = ".metric", values_to = "n")

top3_missing <- sel_procedures[!sel_procedures %in% top_3$strat]

top3_missing <- expand.grid(top3_missing, unique(top_3$modality), unique(top_3$.metric))
colnames(top3_missing) <- c("strat", "modality", ".metric")
top3_missing$n <- 0

top_3 <- bind_rows(top_3, top3_missing)

top1_order <- top_1 |> group_by(strat) |> 
  summarize(mean_n = mean(n)) |> 
  arrange(mean_n) |> 
  pull(strat)

top1_p <- mutate(top_1, strat = factor(strat, levels = top1_order)) |> 
  ggplot(aes(x = strat, y = n)) +
  geom_boxplot() +
  labs(y = "Number of times selection\nmethod is in the best performing",
       title = "Top performing") +
  facet_grid(~.metric) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank())

top3_order <- top_3 |> group_by(strat) |> 
  summarize(mean_n = mean(n)) |> 
  arrange(mean_n) |> 
  pull(strat)
  
top3_p <- mutate(top_3, strat = factor(strat, levels = top3_order)) |> 
  ggplot(aes(x = strat, y = n)) +
  geom_boxplot() +
  labs(x = "Selection procedure",
       y = "Number of times selection\nmethod is in the top 3 best performing",
       title = "Top 3 performing") +
  facet_grid(~.metric) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

pdf("output/v8/paper-figures/imbalance2-figure-barplot.pdf", height = 6, width = 12)
  top1_p / top3_p
dev.off()

