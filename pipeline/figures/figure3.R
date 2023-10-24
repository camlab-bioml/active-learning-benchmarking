suppressPackageStartupMessages({
  library(scater)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(scales)
  library(magick)
  library(tidyr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(colorspace)
  library(ggplot2)
})
source("pipeline/whatsthatcell-helpers.R")

### [ ACCURACIES ] ####
acc <- lapply(snakemake@input$accs, 
              function(x){
  df <- read_tsv(x) |> 
    mutate(cohort = case_when(grepl("CyTOF", basename(x)) ~ "CyTOF",
                              grepl("snRNASeq", basename(x)) ~ "snRNASeq",
                              grepl("scRNASeq", basename(x)) ~ "scRNASeq",
                              grepl("scRNALung", basename(x)) ~ "scRNALung",
                              grepl("tabulaLiver", basename(x)) ~ "tabulaLiver",
                              grepl("tabulaVasc", basename(x)) ~ "tabulaVasc",
                              grepl("liverAtlas-", basename(x)) ~ "liverAtlas"))
  df
}) |> 
  bind_rows() |> 
  mutate(selection_procedure = case_when(selection_procedure == "random" ~ "Random",
                                         selection_procedure == "NoMarkerSeurat_clustering" ~ "AR No Marker",
                                         selection_procedure == "MarkerSeurat_clustering" ~ "AR Marker",
                                         selection_procedure == "highest-entropy-AL" ~ "AL Highest entropy",
                                         selection_procedure == "lowest-maxp-AL" ~ "AL Lowest maxp",
                                         selection_procedure == "0.95-entropy-AL" ~ "AL 0.95-entropy",
                                         selection_procedure == "0.05-maxp-AL" ~ "AL 0.05-maxp",
                                         selection_procedure == "0.25-maxp-AL" ~ "AL 0.25-maxp",
                                         selection_procedure == "0.75-entropy-AL" ~ "AL 0.75-entropy",
                                         TRUE ~ selection_procedure))


sel_meth_cols <- sel_met_cols

col_fun <- colorRamp2(c(-1,0,1), c("skyblue","white", "brown1"))

hm_mat <- acc |> 
  group_by(method, knn, res, cell_num, initial, selection_procedure, AL_alg, .metric, cohort) |>
  summarize(mean.estimate = mean(na.omit(.estimate))) |> 
  ungroup() |> 
  pivot_wider(names_from = .metric, values_from = mean.estimate) |> 
  select(bal_accuracy:sensitivity) |> 
  dplyr::rename("Balanced accuracy" = "bal_accuracy",
                "F1-score" = "f_meas",
                "Kappa" = "kap",
                "MCC" = "mcc",
                "Sensitivity" = "sensitivity") |> 
  as.matrix() |> 
  cor(use = "pairwise.complete.obs")
  
pdf(snakemake@output$hm, height = 6, width = 7)
  Heatmap(hm_mat, 
          col = col_fun,
          name = "Correlation\ncoefficient",
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.1f", hm_mat[i, j]), x, y, gp = gpar(fontsize = 10))
            }
          )
dev.off()


random <- filter(acc, selection_procedure == "Random") |> 
  select(-c(knn, res, rand, corrupted, initial, selection_procedure, AL_alg, .estimator)) |> 
  dplyr::rename("rand_estimate" = ".estimate")

rand_improvement <- acc |> 
  filter(initial == "ranking" | is.na(initial)) |> 
  filter(AL_alg == "rf" | is.na(AL_alg)) |> 
  select(-c(.estimator, rand, corrupted, initial)) |> 
  left_join(random, by = c("method", "cell_num", "seed", ".metric", "cohort")) |> 
  mutate(rand_improvement = (.estimate - rand_estimate) / rand_estimate)


hm <- rand_improvement |> 
  group_by(selection_procedure, .metric, cohort) |> 
  summarize(median_improvement = median(na.omit(rand_improvement))) |> 
  group_by(.metric, cohort) |> 
  arrange(median_improvement) |> 
  mutate(rank = length(unique(rand_improvement$selection_procedure)):1) |> 
  filter(.metric == "bal_accuracy" | .metric == "sensitivity") |> 
  mutate(.metric = case_when(.metric == "bal_accuracy" ~ "Balanced accuracy",
                             .metric == "sensitivity" ~ "Sensitivity"),
         cohort = case_when(cohort == "tabulaVasc" ~ "scRNASeq - Vasculature",
                            cohort == "snRNASeq" ~ "snRNASeq - Pancreas cancer",
                            cohort == "scRNASeq" ~ "scRNASeq - Breast cancer cell lines",
                            cohort == "scRNALung" ~ "scRNASeq - Lung cancer cell lines",
                            cohort == "liverAtlas" ~ "scRNASeq - Liver",
                            cohort == "CyTOF" ~ "CyTOF - Bone marrow")) |> 
  ggplot(aes(x = selection_procedure, y = cohort, fill = as.factor(rank))) +
  geom_tile() +
  scale_fill_viridis_d(direction = -1) +
  labs(x = "Selection procedure", y = "", fill = "Rank") +
  facet_wrap(~.metric, nrow = 1) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


pdf(snakemake@output$overall_fig, height = 4, width = 9)
  hm
dev.off()

