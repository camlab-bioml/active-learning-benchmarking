suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(patchwork)
  library(magick)
  library(ggpubr)
})
devtools::load_all("/ggplot2")
source("pipeline/whatsthatcell-helpers.R")

## Part A: Proportion of cells initially selected
get_proportion_rep <- function(sce_path, rand_files, rank_files, cohort){
  sce <- readRDS(sce_path)
  
  cell_types <- unique(sce$cell_type)
  
  if(is.null(cell_types)){
    cell_types <- unique(sce$CellType)
  }
  
  random <- lapply(rand_files, function(x){
    sel_cells <- read_tsv(x) |>
      filter(iteration == 0) |> 
      pull(cell_type)
    
    sum(cell_types %in% sel_cells) / length(cell_types)
  }) |> unlist()
  
  ranking <- lapply(rank_files, function(x){
    sel_cells <- read_tsv(x) |> 
      filter(iteration == 0) |> 
      pull(cell_type)
    
    sum(cell_types %in% sel_cells) / length(cell_types)
  }) |> unlist()
  
  tibble(prop_sel = c(random, ranking), 
         type = c(rep("Random", length(random)),
                  rep("Ranked", length(ranking))),
         cohort = cohort)
}

props <- bind_rows(
  get_proportion_rep(snakemake@input$cytof,
                     snakemake@input$cytof_AL_rand,
                     snakemake@input$cytof_AL_rank,
                     "CyTOF\nBone marrow"),
  get_proportion_rep(snakemake@input$scrna,
                     snakemake@input$scrna_AL_rand,
                     snakemake@input$scrna_AL_rank,
                     "scRNASeq\nBreast cancer cell lines"),
  get_proportion_rep(snakemake@input$snrna,
                     snakemake@input$snrna_AL_rand,
                     snakemake@input$snrna_AL_rank,
                     "snRNASeq\nPancreas cancer"),
  get_proportion_rep(snakemake@input$snrna_lung,
                     snakemake@input$snrna_lung_AL_rand,
                     snakemake@input$snrna_lung_AL_rank,
                     "scRNASeq\nLung cancer cell lines"),
  get_proportion_rep(snakemake@input$liverAtlas,
                     snakemake@input$liverAtlas_AL_rand,
                     snakemake@input$liverAtlas_AL_rank,
                     "scRNASeq\nLiver"),
  get_proportion_rep(snakemake@input$tabulaVasc,
                     snakemake@input$tabulaVasc_AL_rand,
                     snakemake@input$tabulaVasc_AL_rank,
                     "scRNASeq\nVasculature")
)

prop <- props |> 
  ggplot(aes(x = prop_sel, y = type, fill = type)) +
  geom_boxplot() +
  labs(x = "Proportion of cell types selected in initial training set of 20 cells",
       y = " \nInitial selection\nmethod") +
  scale_fill_manual(values = c("#8F3985", "#98DFEA")) +
  facet_wrap(~cohort, nrow = 2) +
  whatsthatcell_theme() + 
  theme(legend.position = "none")

pdf(snakemake@output$rank_vs_rand_cell_num, height = 2.5, width = 7)
  prop
dev.off()


## Part B: heatmaps
acc <- lapply(snakemake@input$accs, function(x){
  df <- read_tsv(x) |> 
    mutate(cohort = case_when(grepl("CyTOF", basename(x)) ~ "CyTOF",
                              grepl("snRNASeq", basename(x)) ~ "snRNASeq",
                              grepl("scRNASeq", basename(x)) ~ "scRNASeq",
                              grepl("scRNALung", basename(x)) ~ "scRNALung",
                              grepl("liverAtlas", basename(x)) ~ "liverAtlas",
                              grepl("tabulaVasc", basename(x)) ~ "tabulaVasc"))
  df
}) |> bind_rows() |>
  mutate(selection_procedure = case_when(selection_procedure == "NoMarkerSeurat_clustering" ~ "AR NoMarker",
                                         selection_procedure == "MarkerSeurat_clustering" ~ "AR Marker",
                                         selection_procedure == "0.95-entropy-AL" ~ "AL 0.95-entropy",
                                         selection_procedure == "highest-entropy-AL" ~ "AL highest-entropy",
                                         selection_procedure == "lowest-maxp-AL" ~ "AL lowest-maxp",
                                         selection_procedure == "0.05-maxp-AL" ~ "AL 0.05-maxp",
                                         selection_procedure == "0.25-maxp-AL" ~ "AL 0.25-maxp",
                                         selection_procedure == "0.75-entropy-AL" ~ "AL 0.75-entropy",
                                         TRUE ~ selection_procedure))


create_heatmap <- function(acc, sel_cohort, comp, title, legend = FALSE){
  subset_acc <- filter(acc, corrupted == 0) |> 
    filter(rand == 0 | is.na(rand)) |> 
    filter(.metric == 'f_meas' & cohort == sel_cohort) |> 
    filter(selection_procedure != "random" & !grepl("-AR", selection_procedure))
  
  if(comp == "lr_vs_rf"){
    improvement <- subset_acc |> 
      pivot_wider(names_from = "AL_alg", values_from = ".estimate") |> 
      mutate(improvement = (rf - multinom) / multinom) |> 
      filter(improvement != Inf & !is.na(improvement)) |> 
      filter(initial == "random") |> 
      select(-c(multinom, rf, .estimator, .metric, cohort, rand, corrupted, knn, res, initial))
  }else if(comp == "random_vs_ranked"){
    improvement <- subset_acc |> 
      pivot_wider(names_from = "initial", values_from = ".estimate") |> 
      mutate(improvement = (ranking - random) / random) |> 
      filter(improvement != Inf & !is.na(improvement)) |> 
      select(-c(ranking, random, .estimator, .metric, cohort, rand, corrupted, knn, res, AL_alg))
  }
  
  mat <- improvement |> 
    group_by(cell_num, selection_procedure) |> 
    summarize(mean_improvement = mean(improvement)) |> 
    ungroup() |> 
    mutate(mean_improvement = case_when(mean_improvement > 0.2 ~ 0.2,
                                        mean_improvement < -0.2 ~ -0.2,
                                        TRUE ~ mean_improvement)) |> 
    pivot_wider(names_from = selection_procedure, values_from = mean_improvement)
  
  ha <- HeatmapAnnotation(`Cell number` = gsub("_[0-9]", "", mat$cell_num),
                          which = "row",
                          show_legend = legend,
                          col = list(`Cell number` = c(`100` = "#FF99C9", 
                                                       `250` = "#A2C7E5", 
                                                       `500` = "#C1BDDB")))
  
  mat |> 
    as.data.frame() |> 
    column_to_rownames("cell_num") |> 
    as.matrix() |> 
    Heatmap(right_annotation = ha,
            column_title = title,
            show_heatmap_legend = legend,
            row_order = c("100", "250", "500"),
            col = colorRamp2(seq(-0.2, 0.2, length = 3), c("blue", "#EEEEEE", "red")),
            name = "Improvement\nscore")
}

### Ranking vs random
row1 <- create_heatmap(acc, "scRNASeq", "random_vs_ranked", "scRNASeq\nBreast cancer cell lines") + 
  create_heatmap(acc, "snRNASeq", "random_vs_ranked", "snRNASeq\nPancreas cancer") + 
  create_heatmap(acc, "CyTOF", "random_vs_ranked", "CyTOF\nBone marrow", TRUE)

row2 <- create_heatmap(acc, "scRNALung", "random_vs_ranked", "scRNASeq\nLung cancer cell lines") + 
  create_heatmap(acc, "liverAtlas", "random_vs_ranked", "scRNASeq\nLiver") +
  create_heatmap(acc, "tabulaVasc", "random_vs_ranked", "scRNASeq\nVasculature", TRUE)

row1 <- draw(row1, column_title = "Selecting initial cells based on marker expression")
row2 <- draw(row2)

pdf(snakemake@output$row_1_ranking_vs_random, height = 3.5, width = 7.5)
  row1
dev.off()

pdf(snakemake@output$row_2_ranking_vs_random, height = 3.4, width = 7.5)
  row2
dev.off()


### Logistic regression vs random forest
lr_vs_rf_row1 <- create_heatmap(acc, "scRNASeq", "lr_vs_rf", "scRNASeq\nBreast cancer cell lines") + 
  create_heatmap(acc, "snRNASeq", "lr_vs_rf", "snRNASeq\nPancreas cancer") +
  create_heatmap(acc, "CyTOF", "lr_vs_rf", "CyTOF\nBone marrow", TRUE)
lr_vs_rf_row1 <- draw(lr_vs_rf_row1,
                      column_title = "F1-score improvement by random forest compared to logistic regression")

lr_vs_rf_row2 <- create_heatmap(acc, "scRNALung", "lr_vs_rf", "scRNASeq\nLung cancer cell lines") + 
  create_heatmap(acc, "liverAtlas", "lr_vs_rf", "scRNASeq\nLiver") +
  create_heatmap(acc, "tabulaVasc", "lr_vs_rf", "scRNASeq\nVasculature", TRUE)
lr_vs_rf_row2 <- draw(lr_vs_rf_row2)

pdf(snakemake@output$lr_rf_hm_1, height = 3.5, width = 7.5)
  lr_vs_rf_row1
dev.off()

pdf(snakemake@output$lr_rf_hm_2, height = 3.4, width = 7.5)
  lr_vs_rf_row2
dev.off()



### Supplementary figures
# LR vs RF
rf_vs_rf_random_1 <- filter(acc, corrupted == 0) |> 
  filter(rand == 0 | is.na(rand)) |> 
  filter(selection_procedure != "random" & !grepl("-AR", selection_procedure)) |> 
  filter(initial == "random") |> 
  filter(!grepl("-AR", selection_procedure)) |> 
  filter(cohort == "CyTOF" | cohort == "scRNASeq" | cohort == "snRNASeq") |> 
  mutate(AL_alg = case_when(AL_alg == "rf" ~ "RF",
                            AL_alg == "multinom" ~ "LR"),
         cohort = case_when(cohort == "CyTOF" ~ "CyTOF\nBone marrow",
                            cohort == "scRNASeq" ~ "scRNASeq\nBreast cancer cell lines",
                            cohort == "snRNASeq" ~ "snRNASeq\nPancreas cancer")) |> 
  ggplot(aes(x = selection_procedure, y = .estimate, fill = AL_alg)) +
  geom_boxplot() +
  scale_fill_manual(values = al_colours()) +
  labs(x = "Selection procedure", fill = "Active learning\nalgorithm", 
       title = "Performance of selection methods comparing random forest and logistic regression active learning strategies",
       subtitle = "Shown are results obtained when randomly selecting the initial set of cells") +
  facet_grid(.metric ~ cohort + cell_num) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

rf_vs_rf_random_2 <- filter(acc, corrupted == 0) |> 
  filter(rand == 0 | is.na(rand)) |> 
  filter(selection_procedure != "random" & !grepl("-AR", selection_procedure)) |> 
  filter(initial == "random") |> 
  filter(!grepl("-AR", selection_procedure)) |> 
  filter(cohort == "scRNALung" | cohort == "liverAtlas" | cohort == "tabulaVasc") |> 
  mutate(AL_alg = case_when(AL_alg == "rf" ~ "RF",
                            AL_alg == "multinom" ~ "LR"),
         cohort = case_when(cohort == "scRNALung" ~ "scRNASeq\nLung cancer cell lines",
                            cohort == "liverAtlas" ~ "scRNASeq\nLiver",
                            cohort == "tabulaVasc" ~ "scRNASeq\nVasculature")) |> 
  ggplot(aes(x = selection_procedure, y = .estimate, fill = AL_alg)) +
  geom_boxplot() +
  scale_fill_manual(values = al_colours()) +
  labs(x = "Selection procedure", fill = "Active learning\nalgorithm", 
       title = "Performance of selection methods comparing random forest and logistic regression active learning strategies",
       subtitle = "Shown are results obtained when randomly selecting the initial set of cells") +
  facet_grid(.metric ~ cohort + cell_num) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

rf_vs_rf_ranking_1 <- filter(acc, corrupted == 0) |> 
  filter(rand == 0 | is.na(rand)) |> 
  filter(selection_procedure != "random" & !grepl("-AR", selection_procedure)) |> 
  filter(initial == "ranking") |> 
  filter(!grepl("-AR", selection_procedure)) |> 
  filter(cohort == "CyTOF" | cohort == "scRNASeq" | cohort == "snRNASeq") |> 
  mutate(AL_alg = case_when(AL_alg == "rf" ~ "RF",
                            AL_alg == "multinom" ~ "LR"),
         cohort = case_when(cohort == "CyTOF" ~ "CyTOF\nBone marrow",
                            cohort == "scRNASeq" ~ "scRNASeq\nBreast cancer cell lines",
                            cohort == "snRNASeq" ~ "snRNASeq\nPancreas cancer")) |> 
  ggplot(aes(x = selection_procedure, y = .estimate, fill = AL_alg)) +
  geom_boxplot() +
  scale_fill_manual(values = al_colours()) +
  labs(x = "Selection procedure", fill = "Active learning\nalgorithm", 
       title = "Performance of selection methods comparing random forest and logistic regression active learning strategies",
       subtitle = "Shown are results obtained when ranking the initial set of cells") +
  facet_grid(.metric ~ cohort + cell_num) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust =1))


rf_vs_rf_ranking_2 <- filter(acc, corrupted == 0) |> 
  filter(rand == 0 | is.na(rand)) |> 
  filter(selection_procedure != "random" & !grepl("-AR", selection_procedure)) |> 
  filter(initial == "ranking") |> 
  filter(!grepl("-AR", selection_procedure)) |> 
  filter(cohort == "scRNALung" | cohort == "liverAtlas" | cohort == "tabulaVasc") |> 
  mutate(AL_alg = case_when(AL_alg == "rf" ~ "RF",
                            AL_alg == "multinom" ~ "LR"),
         cohort = case_when(cohort == "scRNALung" ~ "scRNASeq\nLung cancer cell lines",
                            cohort == "liverAtlas" ~ "scRNASeq\nLiver",
                            cohort == "tabulaVasc" ~ "scRNASeq\nVasculature")) |> 
  ggplot(aes(x = selection_procedure, y = .estimate, fill = AL_alg)) +
  geom_boxplot() +
  scale_fill_manual(values = al_colours()) +
  labs(x = "Selection procedure", fill = "Active learning\nalgorithm", 
       title = "Performance of selection methods comparing random forest and logistic regression active learning strategies",
       subtitle = "Shown are results obtained when ranking the initial set of cells") +
  facet_grid(.metric ~ cohort + cell_num) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust =1))

pdf(snakemake@output$rf_lr_supp_1, height = 13, width = 15)
  rf_vs_rf_random_1 / rf_vs_rf_ranking_1 + plot_layout(guides = "collect") + 
    plot_annotation(tag_levels = "A")
dev.off()

pdf(snakemake@output$rf_lr_supp_2, height = 13, width = 15)
  rf_vs_rf_random_2 / rf_vs_rf_ranking_2 + plot_layout(guides = "collect") + 
    plot_annotation(tag_levels = "A")
dev.off()



# random vs ranking
rand_vs_rank_lr_1 <- filter(acc, corrupted == 0) |> 
  filter(rand == 0 | is.na(rand)) |> 
  filter(selection_procedure != "random" & !grepl("-AR", selection_procedure)) |> 
  filter(AL_alg == "multinom") |> 
  filter(!grepl("-AR", selection_procedure)) |> 
  filter(cohort == "CyTOF" | cohort == "scRNASeq" | cohort == "snRNASeq") |> 
  mutate(initial = case_when(initial == "random" ~ "Random",
                            initial == "ranking" ~ "Ranking"),
         cohort = case_when(cohort == "CyTOF" ~ "CyTOF\nBone marrow",
                            cohort == "scRNASeq" ~ "scRNASeq\nBreast cancer cell lines",
                            cohort == "snRNASeq" ~ "snRNASeq\nPancreas cancer")) |> 
  ggplot(aes(x = selection_procedure, y = .estimate, fill = initial)) +
  geom_boxplot() +
  labs(x = "Selection procedure", fill = "Initial selection", 
       title = "Performance of selection methods comparing random and ranking based initial cell selections",
       subtitle = "Shown are results obtained using logistic regression as an active learning algorithm") +
  facet_grid(.metric ~ cohort + cell_num) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

rand_vs_rank_lr_2 <- filter(acc, corrupted == 0) |> 
  filter(rand == 0 | is.na(rand)) |> 
  filter(selection_procedure != "random" & !grepl("-AR", selection_procedure)) |> 
  filter(AL_alg == "multinom") |> 
  filter(!grepl("-AR", selection_procedure)) |> 
  filter(cohort == "scRNALung" | cohort == "liverAtlas" | cohort == "tabulaVasc") |> 
  mutate(initial = case_when(initial == "random" ~ "Random",
                             initial == "ranking" ~ "Ranking"),
         cohort = case_when(cohort == "scRNALung" ~ "scRNASeq\nLung cancer cell lines",
                            cohort == "liverAtlas" ~ "scRNASeq\nLiver",
                            cohort == "tabulaVasc" ~ "scRNASeq\nVasculature")) |> 
  ggplot(aes(x = selection_procedure, y = .estimate, fill = initial)) +
  geom_boxplot() +
  labs(x = "Selection procedure", fill = "Initial selection", 
       title = "Performance of selection methods comparing random and ranking based initial cell selections",
       subtitle = "Shown are results obtained using logistic regression as an active learning algorithm") +
  facet_grid(.metric ~ cohort + cell_num) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())


rand_vs_rank_rf_1 <- filter(acc, corrupted == 0) |> 
  filter(rand == 0 | is.na(rand)) |> 
  filter(selection_procedure != "random" & !grepl("-AR", selection_procedure)) |> 
  filter(AL_alg == "rf") |> 
  filter(!grepl("-AR", selection_procedure)) |> 
  filter(cohort == "CyTOF" | cohort == "scRNASeq" | cohort == "snRNASeq") |> 
  mutate(initial = case_when(initial == "random" ~ "Random",
                            initial == "ranking" ~ "Ranking"),
         cohort = case_when(cohort == "CyTOF" ~ "CyTOF\nBone marrow",
                            cohort == "scRNASeq" ~ "scRNASeq\nBreast cancer cell lines",
                            cohort == "snRNASeq" ~ "snRNASeq\nPancreas cancer")) |> 
  ggplot(aes(x = selection_procedure, y = .estimate, fill = initial)) +
  geom_boxplot() +
  labs(x = "Selection procedure", fill = "Initial selection", 
       title = "Performance of selection methods comparing random and ranking based initial cell selections",
       subtitle = "Shown are results obtained using random forest as an active learning algorithm") +
  facet_grid(.metric ~ cohort + cell_num) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust =1))


rand_vs_rank_rf_2 <- filter(acc, corrupted == 0) |> 
  filter(rand == 0 | is.na(rand)) |> 
  filter(selection_procedure != "random" & !grepl("-AR", selection_procedure)) |> 
  filter(AL_alg == "rf") |> 
  filter(!grepl("-AR", selection_procedure)) |> 
  filter(cohort == "scRNALung" | cohort == "liverAtlas" | cohort == "tabulaVasc") |> 
  mutate(initial = case_when(initial == "random" ~ "Random",
                             initial == "ranking" ~ "Ranking"),
         cohort = case_when(cohort == "scRNALung" ~ "scRNASeq\nLung cancer cell lines",
                            cohort == "liverAtlas" ~ "scRNASeq\nLiver",
                            cohort == "tabulaVasc" ~ "scRNASeq\nVasculature")) |> 
  ggplot(aes(x = selection_procedure, y = .estimate, fill = initial)) +
  geom_boxplot() +
  labs(x = "Selection procedure", fill = "Initial selection", 
       title = "Performance of selection methods comparing random and ranking based initial cell selections",
       subtitle = "Shown are results obtained using random forest as an active learning algorithm") +
  facet_grid(.metric ~ cohort + cell_num) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust =1))




pdf(snakemake@output$rand_rank_supp_1, height = 13, width = 15)
  rand_vs_rank_lr_1 / rand_vs_rank_rf_1 + plot_layout(guides = "collect") + 
    plot_annotation(tag_levels = "A")
dev.off()


pdf(snakemake@output$rand_rank_supp_2, height = 13, width = 15)
  rand_vs_rank_lr_2 / rand_vs_rank_rf_2 + plot_layout(guides = "collect") + 
    plot_annotation(tag_levels = "A")
dev.off()
