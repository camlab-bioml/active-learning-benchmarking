suppressPackageStartupMessages({
  library(tidyverse)
  library(ggalluvial)
  library(patchwork)
  library(scales)
})
source("pipeline/whatsthatcell-helpers.R")

### PREDICTIVE LABELLING ACCURACY
pred_lab_acc <- read_tsv("output/v8/new/pred-labeling-accuracy.tsv")#snakemake@input$acc)

# SUPPLEMENTAL FILE
pdf(snakemake@output$supp1, height = 14, width = 10)
  pred_lab_acc |> 
    filter(.metric == "f_meas") |> 
    filter(mod == "scRNASeq" | mod == "scRNALung" | mod == "CyTOF") |> 
    mutate(pred_sel = gsub("top", "", pred_sel),
           pred_sel = paste0(pred_sel, "%")) |> 
    mutate(pred_sel = factor(pred_sel, c("10%", "50%", "100%"))) |> 
    mutate(pred = case_when(pred == "multinom" ~ "LR",
                            pred == "rf" ~ "RF"),
           selection_procedure = case_when(selection_procedure == "MarkerSeurat-clustering" ~ "AR Marker",
                                           selection_procedure == "NoMarkerSeurat-clustering" ~ "AR No Marker",
                                           selection_procedure == "Active-Learning_entropy" ~ "AL Highest entropy",
                                           selection_procedure == "Active-Learning_maxp" ~ "AL Lowest maxp",
                                           selection_procedure == "random" ~ "Random"),
           mod = case_when(mod == "tabulaVasc" ~ "scRNASeq\nVasculature",
                           mod == "snRNASeq" ~ "snRNASeq\nPancreas cancer",
                           mod == "scRNASeq" ~ "scRNASeq\nBreast cancer cell lines",
                           mod == "scRNALung" ~ "scRNASeq\nLung cancer cell lines",
                           mod == "liverAtlas" ~ "scRNASeq\nLiver",
                           mod == "CyTOF" ~ "CyTOF\nBone marrow")) |> 
    ggplot(aes(x = pred_sel, y = .estimate, fill = pred)) +
    geom_boxplot() +
    labs(x = "Percentage of most confidently labelled cells selected", 
         fill = "Self-\ntraining\nalgorithm", y = "F-1 score") +
    scale_fill_manual(values = al_colours()) +
    facet_grid(mod+cell_num~selection_procedure) +
    whatsthatcell_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

pdf(snakemake@output$supp2, height = 14, width = 10)
  pred_lab_acc |> 
    filter(.metric == "f_meas") |> 
    filter(mod == "snRNASeq" | mod == "liverAtlas" | mod == "tabulaVasc") |> 
    mutate(pred_sel = gsub("top", "", pred_sel),
           pred_sel = paste0(pred_sel, "%")) |> 
    mutate(pred_sel = factor(pred_sel, c("10%", "50%", "100%"))) |> 
    mutate(pred = case_when(pred == "multinom" ~ "LR",
                            pred == "rf" ~ "RF"),
           selection_procedure = case_when(selection_procedure == "MarkerSeurat-clustering" ~ "AR Marker",
                                           selection_procedure == "NoMarkerSeurat-clustering" ~ "AR No Marker",
                                           selection_procedure == "Active-Learning_entropy" ~ "AL Highest entropy",
                                           selection_procedure == "Active-Learning_maxp" ~ "AL Lowest maxp",
                                           selection_procedure == "random" ~ "Random"),
           mod = case_when(mod == "tabulaVasc" ~ "scRNASeq\nVasculature",
                           mod == "snRNASeq" ~ "snRNASeq\nPancreas cancer",
                           mod == "scRNASeq" ~ "scRNASeq\nBreast cancer cell lines",
                           mod == "scRNALung" ~ "scRNASeq\nLung cancer cell lines",
                           mod == "liverAtlas" ~ "scRNASeq\nLiver",
                           mod == "CyTOF" ~ "CyTOF\nBone marrow")) |> 
    ggplot(aes(x = pred_sel, y = .estimate, fill = pred)) +
    geom_boxplot() +
    labs(x = "Percentage of most confidently labelled cells selected", 
         fill = "Self-\ntraining\nalgorithm", y = "F-1 score") +
    scale_fill_manual(values = al_colours()) +
    facet_grid(mod+cell_num~selection_procedure) +
    whatsthatcell_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

# PANEL A, SUPPLEMENTAL
pred_lab_acc_pub <- pred_lab_acc |> 
  filter(.metric == "f_meas" & strat == "highest_entropy" &
           initial == "ranking" & al == "rf") |> 
  mutate(pred_sel = gsub("top", "", pred_sel),
         pred_sel = paste0(pred_sel, "%")) |> 
  mutate(pred_sel = factor(pred_sel, c("10%", "50%", "100%"))) |> 
  mutate(pred = case_when(pred == "multinom" ~ "LR",
                          pred == "rf" ~ "RF"))

# PANEL A, MAIN FIGURE
summary_pred_by_percent_included <- pred_lab_acc_pub |> 
  group_by(mod, initial, selection_procedure, strat, al, pred, pred_sel) |> 
  summarize(mean_estimate = mean(na.omit(.estimate))) |> 
  mutate(mod = case_when(mod == "tabulaVasc" ~ "scRNASeq\nVasculature",
                         mod == "snRNASeq" ~ "snRNASeq\nPancreas cancer",
                         mod == "scRNASeq" ~ "scRNASeq\nBreast cancer cell lines",
                         mod == "scRNALung" ~ "scRNASeq\nLung cancer cell lines",
                         mod == "liverAtlas" ~ "scRNASeq\nLiver",
                         mod == "CyTOF" ~ "CyTOF\nBone marrow")) |> 
  ggplot(aes(x = pred_sel, y = mean_estimate, colour = pred, group = pred)) +
  geom_point() + 
  geom_line() +
  scale_colour_manual(values = al_colours()) +
  labs(x = "Percentage of most confidently labeled cells", y = "Mean F1-score",
       colour = "Self-training algorithm") +
  facet_wrap(~mod, nrow = 1) +
  whatsthatcell_theme() +
  theme(legend.position = "bottom")


## Benchmarking with predictive labelling data
scrna <- read_tsv("output/v8/new/pred2/benchmark-predictive-labeling-scRNASeq.tsv") |> #snakemake@input$scrna) |> 
  mutate(AL_alg = sub(".*-ALAlg-", "",  selection_procedure),
         selection_procedure = sub("-ALAlg-.*", "", selection_procedure),
         cohort = "scRNASeq\nBreast cancer cell lines")
snrna <- read_tsv("output/v8/new/pred2/benchmark-predictive-labeling-snRNASeq.tsv") |> #snakemake@input$snrna) |> 
  mutate(AL_alg = sub(".*-ALAlg-", "",  selection_procedure),
         selection_procedure = sub("-ALAlg-.*", "", selection_procedure),
         cohort = "snRNASeq\nPancreas cancer")
cytof <- read_tsv("output/v8/new/pred2/benchmark-predictive-labeling-CyTOF.tsv") |> #snakemake@input$cytof) |> 
  mutate(AL_alg = sub(".*-ALAlg-", "",  selection_procedure),
         selection_procedure = sub("-ALAlg-.*", "", selection_procedure),
         cohort = "CyTOF\nBone marrow")
scrnaLung <- read_tsv("output/v8/new/pred2/benchmark-predictive-labeling-scRNALung.tsv") |> 
  mutate(AL_alg = sub(".*-ALAlg-", "",  selection_procedure),
         selection_procedure = sub("-ALAlg-.*", "", selection_procedure),
         cohort = "scRNASeq\nLung cancer cell lines")
liverAtlas <- read_tsv("output/v8/new/pred2/benchmark-predictive-labeling-tabulaVasc.tsv") |> 
  mutate(AL_alg = sub(".*-ALAlg-", "",  selection_procedure),
         selection_procedure = sub("-ALAlg-.*", "", selection_procedure),
         cohort = "scRNASeq\nLiver")
tabulaVasc <- read_tsv("output/v8/new/pred2/benchmark-predictive-labeling-tabulaVasc.tsv") |> 
  mutate(AL_alg = sub(".*-ALAlg-", "",  selection_procedure),
         selection_procedure = sub("-ALAlg-.*", "", selection_procedure),
         cohort = "scRNASeq\nVasculature")

acc <- bind_rows(scrna, snrna, cytof, scrnaLung, liverAtlas, tabulaVasc)

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
  filter(selection_procedure == "highest-entropy-AL" & cell_num == 100, .metric == "f_meas") |> 
  ggplot(aes(x = method, y = gap, fill = pred_labeller)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#DA94D4", "#7EA3CC")) +
  labs(x = "Cell type assignment method", y = "% change in F1-score",
       fill = "Self-training algorithm") +
  facet_wrap(~cohort + cell_selection, scales = "free", nrow = 2) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        legend.position = "bottom")

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
  filter(cohort == "CyTOF\nBone marrow") |> 
  plot_comp(1, ylab = "Label") +
  labs(title = "CyTOF - Bone marrow")
scrnaseq_gap <- comp_dataset |> 
  filter(cohort == "scRNASeq\nBreast cancer cell lines") |> 
  plot_comp(2, xlab = "label") +
  labs(title = "scRNASeq - Breast cancer cell lines")
snrnaseq_gap <- comp_dataset |> 
  filter(cohort == "snRNASeq\nPancreas cancer") |> 
  plot_comp(2) +
  labs(title = "snRNASeq - Pancreas cancer")
scrnalung_gap <- comp_dataset |> 
  filter(cohort == "scRNASeq\nLung cancer cell lines") |> 
  plot_comp(2, ylab = "Label") + 
  labs(title = "scRNASeq - Lung cancer cell lines")
liverAtlas_gap <- comp_dataset |> 
  filter(cohort == "scRNASeq\nLiver") |> 
  plot_comp(2, ylab = "Label") + 
  labs(title = "scRNASeq - Liver")
tabulaVasc_gap <- comp_dataset |> 
  filter(cohort == "scRNASeq\nVasculature") |> 
  plot_comp(2, xlab = "label") +
  labs(title = "scRNASeq - Vasculature")

gap_comb <- ((cytof_gap | scrnaseq_gap | liverAtlas_gap)) /
  ((scrnalung_gap | tabulaVasc_gap | snrnaseq_gap)) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Detecting mislabelled cells
cytof <- list.files("output/v8/identify_mislabelled/CyTOF", full.names = TRUE) #snakemake@input$mislabelled_cytof
scrna <- list.files("output/v8/identify_mislabelled/scRNASeq", full.names = TRUE) #snakemake@input$mislabelled_scrna
snrna <- list.files("output/v8/identify_mislabelled/snRNASeq", full.names = TRUE) #snakemake@input$mislabelled_snrna
scrnaLung <- list.files("output/v8/identify_mislabelled/scRNALung", full.names = TRUE) # snakemake@input$mislabelled_scrnalung
scrnaVasc <- list.files("output/v8/identify_mislabelled/tabulaVasc", full.names = TRUE) # snakemake@input$mislabelled_scrnalung
scrnaLiver <- list.files("output/v8/identify_mislabelled/liverAtlas//", full.names = TRUE) # snakemake@input$mislabelled_scrnalung


mislabelled_pred <- lapply(c(cytof, scrna, snrna, scrnaLung, scrnaVasc, scrnaLiver), function(x){
  df <- read_tsv(x)
  probs <- select(df, -c(cell_id, pred_type, corr_cell_type, gt_cell_type, params))
  
  df$entropy <- apply(probs, 1, calculate_entropy)
  df$entropy <- df$entropy / log(ncol(probs), 2)
  
  select(df, cell_id, entropy, pred_type, corr_cell_type, gt_cell_type, params)
}) |> bind_rows() |> 
  separate(params, c("rm_mod", "modality", "rm_pred", "predAlg", "rm_seed", "seed")) |> 
  select(-starts_with("rm")) 

mislabelled_pred_plot1 <- mislabelled_pred |> 
  filter(modality == "CyTOF" | modality == "tabulaVasc" | modality == "liverAtlas") |>
  mutate(cell_is_corrupt = corr_cell_type != gt_cell_type,
         predAlg = case_when(predAlg == "multinom" ~ "LR",
                             predAlg == "rf" ~ "RF"),
         modality = case_when(modality == "CyTOF" ~ "CyTOF - Bone marrow",
                              modality == "liverAtlas" ~ "scRNASeq - Liver",
                              modality == "tabulaVasc" ~ "scRNASeq - Vasculature")) |> 
  ggplot(aes(x = gt_cell_type, y = entropy, fill = cell_is_corrupt)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#8B80F9", "#F03560")) +
  labs(x = "Ground truth cell type", y = "Scaled entropy", fill = "Cell corrupted during training") +
  facet_grid(predAlg~modality, scales = "free") +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank())

mislabelled_pred_plot2 <- mislabelled_pred |> 
  filter(modality == "scRNASeq" | modality == "snRNASeq" | modality == "scRNALung") |> 
  mutate(cell_is_corrupt = corr_cell_type != gt_cell_type,
         predAlg = case_when(predAlg == "multinom" ~ "LR",
                             predAlg == "rf" ~ "RF"),
         modality = case_when(modality == "scRNASeq" ~ "scRNASeq - Breast cancer cell lines",
                              modality == "snRNASeq" ~ "snRNASeq - Pancreas cancer",
                              modality == "scRNALung" ~ "scRNASeq - Lung cancer cell lines")) |> 
  ggplot(aes(x = gt_cell_type, y = entropy, fill = cell_is_corrupt)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#8B80F9", "#F03560")) +
  labs(x = "Ground truth cell type", y = "Scaled entropy", fill = "Cell corrupted during training") +
  facet_grid(predAlg~modality, scales = "free") +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


mislabelled_pred_plot <- mislabelled_pred_plot1 / mislabelled_pred_plot2 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

pdf("output/v8/paper-figures/pred-labelling-test.pdf", height = 27, width = 16)
  (summary_pred_by_percent_included /
    lr_vs_rf /
    wrap_elements(full = gap_comb + plot_layout(guides = "collect")) /
    wrap_elements(full = mislabelled_pred_plot + plot_layout(guides = "collec"))) +
    plot_layout(heights = c(0.45, 1.6, 3.7, 3.6)) +
    plot_annotation(tag_levels = "A")
dev.off()

pdf("output/v8/paper-figures/pred-labelling-test.pdf", #snakemake@output$main, 
    height = 20, width = 12)
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
                                   pred_labeller == "rf" ~ "RF")) |> 
  mutate(cohort = gsub("\n", " - ", cohort))
  
plot_sup_gap <- function(df, sel_cohort){
  filter(df, cohort == sel_cohort & .metric == "f_meas") |> 
    mutate(cell_selection = gsub("top", "", cell_selection),
           cell_selection = paste0(cell_selection, "%"),
           cell_selection = factor(cell_selection, levels = c("10%", "50%", "100%"))) |> 
    ggplot(aes(x = method, y = gap, fill = pred_labeller)) +
    geom_boxplot() +
    scale_fill_manual(values = c("#DA94D4", "#7EA3CC")) +
    labs(x = "Cell type assignment method", y = "% change in F1-score",
         fill = "Self-\ntraining\nalgorithm", title = sel_cohort) +
    facet_grid(cell_selection + cell_num~selection_procedure, scales = "free") +
    whatsthatcell_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}

pdf("output/v8/paper-figures/Supp-CyTOF-f1-improvement.pdf", #snakemake@output$sup_cytof, 
    width = 12, height = 8)
  plot_sup_gap(sup_acc_gap, "CyTOF - Bone marrow")
dev.off()

pdf("output/v8/paper-figures/Supp-scRNASeq-f1-improvement.pdf", #snakemake@output$sup_scrna, 
    width = 15, height = 9)
  plot_sup_gap(sup_acc_gap, "scRNASeq - Breast cancer cell lines")
dev.off()

pdf("output/v8/paper-figures/Supp-snRNASeq-f1-improvement.pdf", #snakemake@output$sup_snrna,
    width = 15, height = 9)
  plot_sup_gap(sup_acc_gap, "snRNASeq - Pancreas cancer")
dev.off()


pdf("output/v8/paper-figures/Supp-scRNALung-f1-improvement.pdf", #snakemake@output$sup_scrna_lung,
    width = 15, height = 9)
  plot_sup_gap(sup_acc_gap, "scRNASeq - Lung cancer cell lines")
dev.off()

pdf("output/v8/paper-figures/Supp-liverAtlas-f1-improvement.pdf", #snakemake@output$sup_liverAtlas,
    width = 15, height = 9)
  plot_sup_gap(sup_acc_gap, "scRNASeq - Liver")
dev.off()

pdf("output/v8/paper-figures/Supp-tabulaVasc-f1-improvement.pdf", #snakemake@output$sup_tabulaVasc,
    width = 15, height = 9)
  plot_sup_gap(sup_acc_gap, "scRNASeq - Vasculature")
dev.off()



# Initial cell number annotation dependence
pdf("output/v8/paper-figures/supp-pred-lab-num-cells.pdf", height = 5, width = 12)
  acc_gap |> 
    filter(.metric == "f_meas" & selection_procedure == "highest-entropy-AL") |> 
    mutate(gap = abs(gap),
           cell_selection = case_when(cell_selection == "top10" ~ "10%",
                                      cell_selection == "top50" ~ "50%",
                                      cell_selection == "top100" ~ "100%"),
           cell_selection = factor(cell_selection, levels = c("10%", "50%", "100%"))) |> 
    ggplot(aes(x = method, y = gap, fill = as.factor(cell_num))) +
    geom_boxplot() +
    labs(x = "Cell type assignment method", fill = "Number of cells\nin self-learning\ntraining set",
         y = "Absolute change in F1-score for self-learning over baseline") +
    facet_wrap(~cell_selection) +
    whatsthatcell_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()

