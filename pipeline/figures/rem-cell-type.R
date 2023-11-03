suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})
source("pipeline/whatsthatcell-helpers.R")

### Summary plots
cytof <- read_tsv(snakemake@input$cytof_summary) |> 
  mutate(cohort = "CyTOF")
scrnaseq <- read_tsv(snakemake@input$scrna_summary) |> 
  mutate(cohort = "scRNASeq")
snrnaseq <- read_tsv(snakemake@input$snrna_summary) |> 
  mutate(cohort = "snRNASeq")
scrnalung <- read_tsv(snakemake@input$scrna_lung_summary) |> 
  mutate(cohort = "scRNALung")
liverAtlas <- read_tsv(snakemake@input$liver_summary) |> 
  mutate(cohort = "liverAtlas")
tabulavasc <- read_tsv(snakemake@input$vasc_summary) |> 
  mutate(cohort = "tabulaVasc")

summary <- bind_rows(cytof, scrnaseq, snrnaseq, scrnalung, liverAtlas, tabulavasc) |> 
  mutate(AL_alg = case_when(AL_alg == "rf" ~ "RF",
                            AL_alg == "multinom" ~ "LR"),
         cohort = case_when(cohort == "tabulaVasc" ~ "scRNASeq\nVasculature",
                            cohort == "snRNASeq" ~ "snRNASeq\nPancreas cancer",
                            cohort == "scRNASeq" ~ "scRNASeq\nBreast cancer cell lines",
                            cohort == "scRNALung" ~ "scRNASeq\nLung cancer cell lines",
                            cohort == "liverAtlas" ~ "scRNASeq\nLiver",
                            cohort == "CyTOF" ~ "CyTOF\nBone marrow")) |> 
  ggplot(aes(as.factor(num_missing_cells), y = criterion_val, colour = gt_rem, group = gt_rem)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 4) +
  ylim(0, 1) +
  scale_colour_manual(values = c("#ffa600", "#9763ff")) +
  labs(x = "Number of missing cells", y = "Scaled median entropy") +
  facet_grid(AL_alg ~ cohort) +
  whatsthatcell_theme() +
  theme(legend.position = "bottom",
        legend.title = element_blank())

# Read in data
### scRNASeq
scrna <- lapply(snakemake@input$scrna, read_tsv) |> 
  bind_rows() |> 
  mutate(al = case_when(grepl('rf', params) ~ "rf",
                        grepl('multinom', params) ~ "multinom"),
         params = gsub(".*rem_celltype-", "", params),
         params = gsub("-seed-[0-9]", "", params))


scrna_box <- scrna |> 
  mutate(al = case_when(al == "multinom" ~ "LR",
                        al == "rf" ~ "RF")) |> 
  ggplot(aes(x = gt_cell_type, y = criterion_val, 
             fill = as.character(num_missing_cells))) +
  geom_boxplot(outlier.colour = "lightgrey", outlier.size = 0.4) +
  labs(y = "Scaled entropy", x = "Ground truth cell type",
       fill = "Number of cells of the type removed in dataset") +
  facet_grid(al ~ params) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


### snRNASeq
snrna <- c(
  snakemake@input$snrna,
  snakemake@input$snrna_group
)

snrna <- lapply(snrna, read_tsv) |> 
  bind_rows() |> 
  mutate(al = case_when(grepl('rf', params) ~ "rf",
                        grepl('multinom', params) ~ "multinom"),
         params = gsub(".*rem_celltype-", "", params),
         params = gsub("-seed-[0-9]", "", params))

snrna_box <- snrna |> 
  mutate(params = factor(params, levels = c("Ductal", "Endothelial", "Schwann", 
                                            "Endothelial, Schwann"))) |> 
  ggplot(aes(x = gt_cell_type, y = criterion_val, 
             fill = as.character(num_missing_cells))) +
  geom_boxplot(outlier.colour = "lightgrey", outlier.size = 0.4) +
  labs(y = "Scaled entropy", x = "Ground truth cell type",
       fill = "Number of cells of the type removed in dataset") +
  facet_wrap(~ params) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

### Cytof
cytof <- lapply(snakemake@input$cytof, read_tsv) |> 
  bind_rows() |> 
  mutate(al = case_when(grepl('rf', params) ~ "rf",
                        grepl('multinom', params) ~ "multinom"),
         params = gsub(".*rem_celltype-", "", params),
         params = gsub("-seed-[0-9]", "", params))

cytof_box <- cytof |> 
  ggplot(aes(x = gt_cell_type, y = criterion_val, 
             fill = as.character(num_missing_cells))) +
  geom_boxplot(outlier.colour = "lightgrey", outlier.size = 0.4) +
  labs(y = "Scaled entropy", x = "Ground truth cell type",
       fill = "Number of cells of the type removed in dataset") +
  facet_wrap(~ params, ncol = 1) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


## Combine
entropies <- (scrna_box | snrna_box | cytof_box) +
  plot_layout(guides = "collect", widths = c(2, 2, 1)) &
  theme(legend.position = 'bottom', plot.tag = element_text(size = 18))

pdf(snakemake@output$main, height = 11, width = 14)
  summary / entropies + 
    plot_layout(heights = c(0.9, 1)) +
    plot_annotation(tag_levels = "A")
dev.off()



## Supplemental for RF for snRNASeq and CyTOF
snrna_sup <- c(
  snakemake@input$snrna_supp,
  snakemake@input$snrna_supp_group
)

snrna_sup <- lapply(snrna_sup, read_tsv) |> 
  bind_rows() |> 
  mutate(al = case_when(grepl('rf', params) ~ "rf",
                        grepl('multinom', params) ~ "multinom"),
         params = gsub(".*rem_celltype-", "", params),
         params = gsub("-seed-[0-9]", "", params))


cytof_sup <- lapply(snakemake@input$cytof_supp, read_tsv) |> 
  bind_rows() |> 
  mutate(al = case_when(grepl('rf', params) ~ "rf",
                        grepl('multinom', params) ~ "multinom"),
         params = gsub(".*rem_celltype-", "", params),
         params = gsub("-seed-[0-9]", "", params))


snrna_sup_box <- snrna_sup |> 
  mutate(params = factor(params, levels = c("Ductal", "Endothelial", "Schwann", 
                                            "Endothelial, Schwann"))) |> 
  ggplot(aes(x = gt_cell_type, y = criterion_val, 
             fill = as.character(num_missing_cells))) +
  geom_boxplot(outlier.colour = "lightgrey", outlier.size = 0.4) +
  labs(y = "Scaled entropy", x = "Ground truth cell type",
       fill = "Number of cells of the type removed in dataset") +
  facet_wrap(~ params) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


cytof_sup_box <- cytof_sup |> 
  ggplot(aes(x = gt_cell_type, y = criterion_val, 
             fill = as.character(num_missing_cells))) +
  geom_boxplot(outlier.colour = "lightgrey", outlier.size = 0.4) +
  labs(y = "Scaled entropy", x = "Ground truth cell type",
       fill = "Number of cells of the type removed in dataset") +
  facet_wrap(~ params, ncol = 1) +
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

pdf(snakemake@output$supp, height = 6, width = 10)
  (snrna_sup_box | cytof_sup_box) +
    plot_layout(guides = "collect", widths = c(2, 1)) +
    plot_annotation(tag_levels = "A") & theme(legend.position = 'bottom')
dev.off()
