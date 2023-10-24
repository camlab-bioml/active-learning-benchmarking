suppressPackageStartupMessages({
  library(tidyverse)
  library(ggalluvial)
  library(patchwork)
  library(scales)
})
source("pipeline/whatsthatcell-helpers.R")

### PREDICTIVE LABELLING ACCURACY
pred_lab_acc <- read_tsv(snakemake@input$acc)

# SUPPLEMENTAL FILE
pdf(snakemake@output$supp, height = 10, width = 12)
  pred_lab_acc |> 
    filter(.metric == "f_meas") |> 
    mutate(pred_sel = gsub("top", "", pred_sel),
           pred_sel = paste0(pred_sel, "%")) |> 
    mutate(pred_sel = factor(pred_sel, c("10%", "50%", "100%"))) |> 
    mutate(pred = case_when(pred == "multinom" ~ "LR",
                            pred == "rf" ~ "RF"),
           selection_procedure = case_when(selection_procedure == "MarkerSeurat-clustering" ~ "AR Marker",
                                           selection_procedure == "NoMarkerSeurat-clustering" ~ "AR No Marker",
                                           selection_procedure == "Active-Learning_entropy" ~ "AL Highest entropy",
                                           selection_procedure == "Active-Learning_maxp" ~ "AL Lowest maxp",
                                           selection_procedure == "random" ~ "Random")) |> 
    ggplot(aes(x = pred_sel, y = .estimate, fill = pred)) +
    geom_boxplot() +
    labs(x = "Percentage of most confidently labelled cells selected", 
         fill = "Self-\ntraining\nalgorithm", y = "F-1 score") +
    scale_fill_manual(values = al_colours()) +
    facet_grid(selection_procedure~mod+cell_num) +
    whatsthatcell_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

# PANEL A, MAIN FIGURE
pred_lab_acc_pub <- pred_lab_acc |> 
  filter(.metric == "f_meas" & selection_procedure == "MarkerSeurat-clustering") |> 
  mutate(pred_sel = gsub("top", "", pred_sel),
         pred_sel = paste0(pred_sel, "%")) |> 
  mutate(pred_sel = factor(pred_sel, c("10%", "50%", "100%"))) |> 
  mutate(pred = case_when(pred == "multinom" ~ "LR",
                          pred == "rf" ~ "RF"))
  
plot_pred_lab_acc <- function(acc, cohort, x_lab = FALSE){
  if(x_lab){
    x_lab <- "Percentage of most confidently\nlabelled cells selected"
  }else{
    x_lab <- ""
  }
  
  filter(acc, mod == cohort) |> 
    ggplot(aes(x = pred_sel, y = .estimate, fill = pred)) +
    geom_boxplot() +
    labs(x = x_lab, fill = "Self-\ntraining\nalgorithm",
         title = cohort, y = "F1-score") +
    scale_fill_manual(values = al_colours()) +
    facet_wrap(~cell_num, nrow = 1) +
    whatsthatcell_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}

pred_lab_acc_plot <- plot_pred_lab_acc(pred_lab_acc_pub, "CyTOF") |
  plot_pred_lab_acc(pred_lab_acc_pub, "scRNASeq", TRUE) & theme(axis.title.y = element_blank()) |
  plot_pred_lab_acc(pred_lab_acc_pub, "snRNASeq") & theme(axis.title.y = element_blank())



## Benchmarking with predictive labelling data
scrna <- read_tsv(snakemake@input$scrna) |> 
  mutate(AL_alg = sub(".*-ALAlg-", "",  selection_procedure),
         selection_procedure = sub("-ALAlg-.*", "", selection_procedure),
         cohort = "scRNASeq")
snrna <- read_tsv(snakemake@input$snrna) |> 
  mutate(AL_alg = sub(".*-ALAlg-", "",  selection_procedure),
         selection_procedure = sub("-ALAlg-.*", "", selection_procedure),
         cohort = "snRNASeq")

cytof <- read_tsv(snakemake@input$cytof) |> 
  mutate(AL_alg = sub(".*-ALAlg-", "",  selection_procedure),
         selection_procedure = sub("-ALAlg-.*", "", selection_procedure),
         cohort = "CyTOF")

acc <- bind_rows(scrna, snrna, cytof)

cell_nums <- as.character(sort(unique(acc$cell_num)))
acc$cell_num <- factor(acc$cell_num, levels = cell_nums)


### Main figure - accuracies
baseline_acc <- acc |> 
  select(-rand, -corrupted, -.estimator, -pred_lab_alg) |> 
  filter(cell_selection == "baseline") |> 
  unite("method", c(method, knn, res, cell_num, initial, seed, 
                    selection_procedure, .metric, AL_alg, cohort), sep = ",")

pred_lab_acc <- acc |> 
  select(-rand, -corrupted, -.estimator) |> 
  filter(cell_selection != "baseline") |> 
  unite("method", c(method, knn, res, cell_num, initial, seed, selection_procedure, 
                    .metric, AL_alg, cohort), sep = ',') |> 
  pivot_wider(names_from = pred_lab_alg, values_from = .estimate)

acc_gap <- left_join(pred_lab_acc, 
                     select(baseline_acc, -cell_selection), 
                     by = "method") |> 
  pivot_longer(c(multinom, rf), values_to = ".estimate_pred_lab", names_to = "pred_labeller") |> 
  mutate(gap = ((.estimate_pred_lab / .estimate) * 100) - 100) |> 
  separate(method, c("method", "knn", "res", "cell_num", "initial", "seed", 
                     "selection_procedure", ".metric", "AL_alg", "cohort"), sep = ",")


# Which is better rf or LR?
lr_vs_rf <- acc_gap |> 
  mutate(pred_labeller = case_when(pred_labeller == "multinom" ~ "LR",
                                   pred_labeller == "rf" ~ "RF"),
         cell_selection = gsub("top", "", cell_selection),
         cell_selection = paste0(cell_selection, "%"),
         cell_selection = factor(cell_selection, levels = c("10%", "50%", "100%"))) |> 
  filter(selection_procedure == "MarkerSeurat-clustering" & cell_num == 100, .metric == "f_meas") |> 
  ggplot(aes(x = method, y = gap, fill = pred_labeller)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#DA94D4", "#7EA3CC")) +
  labs(x = "Cell type assignment method", y = "% change in F1-score",
       fill = "Self-\ntraining\nalgorithm") +
  facet_wrap(~cohort + cell_selection, scales = "free", nrow = 1) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Compare improvement to original accuracy by selection procedure
comp_dataset <- acc_gap |> 
  filter(.metric == "f_meas", cell_num == 100) |> 
  filter(initial == "random" | is.na(initial) | initial == "NA") |> 
  filter(AL_alg == 'rf' | is.na(AL_alg) | AL_alg == "NA") |> 
  filter(pred_labeller == "multinom" & cell_selection == "top50") |> 
  mutate(selection_procedure = case_when(selection_procedure == "random" ~ "Random",
                                         selection_procedure == "MarkerSeurat-clustering" ~ "AR Marker",
                                         selection_procedure == "NoMarkerSeurat-clustering" ~ "AR No Marker",
                                         selection_procedure == "0.05-maxp-AL" ~ "AL 0.05-maxp",
                                         selection_procedure == "0.25-maxp-AL" ~ "AL 0.25-maxp",
                                         selection_procedure == "0.75-entropy-AL" ~ "AL 0.75-entropy",
                                         selection_procedure == "0.95-entropy-AL" ~ "AL 0.95-entropy",
                                         selection_procedure == "highest-entropy-AL" ~ "AL Highest entropy",
                                         selection_procedure == "lowest-maxp-AL" ~ "AL Lowest maxp"))

plot_comp <- function(df, ncol, ylab = "", xlab = ""){
  if(ylab != ""){
    ylab <- "Baseline F1-score"
  }
  if(xlab != ""){
    xlab <- "% Self-training improvement"
  }
  df |> 
    ggplot(aes(x = gap, y = .estimate, colour = selection_procedure)) +
    geom_point() +
    scale_color_manual(values = sel_met_cols) +
    labs(x = xlab,
         y = ylab,
         colour = "Selection procedure") +
    whatsthatcell_theme() +
    facet_wrap(~method, ncol = ncol)
}

cytof_gap <- comp_dataset |> 
  filter(cohort == "CyTOF") |> 
  plot_comp(1, ylab = "Label") +
  labs(title = "CyTOF")

scrnaseq_gap <- comp_dataset |> 
  filter(cohort == "scRNASeq") |> 
  plot_comp(2, xlab = "label") +
  labs(title = "scRNASeq")

snrnaseq_gap <- comp_dataset |> 
  filter(cohort == "snRNASeq") |> 
  plot_comp(2) +
  labs(title = "snRNASeq")

gap_comb <- (cytof_gap | scrnaseq_gap | snrnaseq_gap) + 
  plot_layout(guides = "collect", widths = c(0.65, 1, 1))

# Detecting mislabelled cells
cytof <- snakemake@input$mislabelled_cytof
scrna <- snakemake@input$mislabelled_scrna
snrna <- snakemake@input$mislabelled_snrna

mislabelled_pred <- lapply(c(cytof, scrna, snrna), function(x){
  df <- read_tsv(x)
  probs <- select(df, -c(cell_id, pred_type, corr_cell_type, gt_cell_type, params))
  
  df$entropy <- apply(probs, 1, calculate_entropy)
  df$entropy <- df$entropy / log(ncol(probs), 2)
  
  select(df, cell_id, entropy, pred_type, corr_cell_type, gt_cell_type, params)
}) |> bind_rows() |> 
  separate(params, c("rm_mod", "modality", "rm_pred", "predAlg", "rm_seed", "seed")) |> 
  select(-starts_with("rm"))

mislabelled_pred_plot <- mislabelled_pred |> 
  mutate(cell_is_corrupt = corr_cell_type != gt_cell_type,
         predAlg = case_when(predAlg == "multinom" ~ "LR",
                             predAlg == "rf" ~ "RF")) |> 
  ggplot(aes(x = gt_cell_type, y = entropy, fill = cell_is_corrupt)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#8B80F9", "#F03560")) +
  labs(x = "Ground truth cell type", y = "Scaled entropy", fill = "Cell corrupted\nduring training") +
  facet_grid(predAlg~modality, scales = "free") +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

pdf(snakemake@output$main, height = 14, width = 12)
  (wrap_elements(full = pred_lab_acc_plot + plot_layout(guides = "collect"))) /
    wrap_elements(full = lr_vs_rf & labs(title = "")) /
    wrap_elements(full = gap_comb) /
    mislabelled_pred_plot +
    plot_layout(heights = c(1, 1.1, 1.3, 1.4)) +
    plot_annotation(tag_levels = "A")
dev.off()


sup_acc_gap <- acc_gap |> 
  mutate(selection_procedure = case_when(selection_procedure == "random" ~ "Random",
                                         selection_procedure == "MarkerSeurat-clustering" ~ "AR Marker",
                                         selection_procedure == "NoMarkerSeurat-clustering" ~ "AR No Marker",
                                         selection_procedure == "0.05-maxp-AL" ~ "AL 0.05 maxp",
                                         selection_procedure == "0.25-maxp-AL" ~ "AL 0.25 maxp",
                                         selection_procedure == "0.75-entropy-AL" ~ "AL 0.75 entropy",
                                         selection_procedure == "0.95-entropy-AL" ~ "AL 0.95 entropy",
                                         selection_procedure == "highest-entropy-AL" ~ "AL highest entropy",
                                         selection_procedure == "lowest-maxp-AL" ~ "AL lowest maxp")) |> 
  mutate(pred_labeller = case_when(pred_labeller == "multinom" ~ "LR",
                                   pred_labeller == "rf" ~ "RF"))
  
plot_sup_gap <- function(df, sel_cohort){
  filter(df, cohort == sel_cohort & .metric == "f_meas") |> 
    mutate(cell_selection = gsub("top", "", cell_selection),
           cell_selection = paste0(cell_selection, "%"),
           cell_selection = factor(cell_selection, levels = c("10%", "50%", "100%"))) |> 
    ggplot(aes(x = method, y = gap, fill = pred_labeller)) +
    geom_boxplot() +
    scale_fill_manual(values = c("#DA94D4", "#7EA3CC")) +
    labs(x = "Cell type assignment method", y = "% change in F1-score",
         fill = "Self-\ntraining\nalgorithm") +
    facet_grid(cell_selection + cell_num~selection_procedure, scales = "free") +
    whatsthatcell_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}

pdf(snakemake@output$sup_cytof, width = 12, height = 8)
  plot_sup_gap(sup_acc_gap, "CyTOF")
dev.off()

pdf(snakemake@output$sup_scrna, width = 12, height = 8)
  plot_sup_gap(sup_acc_gap, "scRNASeq")
dev.off()

pdf(snakemake@output$sup_snrna, width = 12, height = 8)
  plot_sup_gap(sup_acc_gap, "snRNASeq")
dev.off()

# ### As function of number of cells predictively labelled
# test <- bind_rows(
#   select(pred_lab_acc, method, cell_selection, multinom) |> 
#     dplyr::rename(".estimate" = "multinom"),
#   baseline_acc
# ) |> 
#   filter(cell_selection != "top200") |> 
#   mutate(cell_selection = factor(cell_selection, c("baseline", "top10", "top50", "top100"))) |> 
#   #filter(method == "Random-Forest,10,0.4,100,NA,0,MarkerSeurat-clustering,kap,NA,scRNASeq") |>
#   separate(method, c("alg", "knn", "res", "cell_num", "initial", "seed", 
#                      "selection_procedure", ".metric", "AL_alg", "cohort"), sep = ",",
#            remove = FALSE) |> 
#   filter(.metric == "sensitivity") |> 
#   filter(initial == "NA" | is.na(initial) | initial == "random") |> 
#   filter(AL_alg == "NA" | is.na(AL_alg) | AL_alg == "multinom")
  
# test |> 
#   filter(alg == "Random-Forest", selection_procedure == "random") |> 
#   ggplot(aes(x = cell_selection, y = .estimate, group = method, colour = cell_num)) +
#   geom_point() +
#   geom_line() +
#   facet_grid(selection_procedure ~ cohort + alg) +
#   whatsthatcell_theme()

# left_join(pred_lab_acc, 
#           select(baseline_acc, -cell_selection), 
#           by = "method") |> 
#   #pivot_wider(names_from = "cell_selection", values_from = "multinom")
#   pivot_longer()

# med_gap <- acc_gap |> 
#   filter(.metric == "f_meas") |> 
#   filter(initial == "NA" | is.na(initial) | initial == "random") |> 
#   filter(AL_alg == "NA" | is.na(AL_alg) | AL_alg == "multinom") |> 
#   group_by(method, cell_num, selection_procedure, cohort, cell_selection, pred_labeller) |> 
#   summarize(median_gap = median(na.omit(gap)))

# med_gap |> 
#   filter(pred_labeller == "multinom" & method ) |> 
#   ggplot(aes(x = cell_selection, y = median_gap, colour = cell_num)) +
#   geom_line() +
#   facet_grid(selection_procedure ~ cohort + method) +
#   whatsthatcell_theme()


# med_gap |> 
#   filter(cell_selection != "top200") |> 
#   mutate(cell_selection = factor(cell_selection, levels = c("top10", "top50", "top100"))) |> 
#   filter(pred_labeller == "multinom" & method == "CyTOF-LDA" & cell_num == 100, cohort == "CyTOF" &
#          selection_procedure == "0.05-maxp-AL") |> 
#   ggplot(aes(x = cell_selection, y = median_gap)) +
#   geom_point() 
#   #geom_line()# +
#   # facet_grid(selection_procedure ~ cohort + method) +
#   # whatsthatcell_theme()

