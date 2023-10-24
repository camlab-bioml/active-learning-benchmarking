suppressPackageStartupMessages({
  library(scater)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(scales)
  library(magick)
  library(ggrepel)
})
source("pipeline/whatsthatcell-helpers.R")


### [ FIGURE 1A - BENCHMARKING OVERVIEW ] #####
# had to remove because magick causes issues with docker
# schematic <- image_read(snakemake@input$schematic) |>
#   image_ggplot()

### [ FIGURE 1B - DATASET COMPOSOTION ] #####
set.seed(42)

CyTOF <- readRDS(snakemake@input$cytof)
CyTOF <- runTSNE(CyTOF)
scRNA <- readRDS(snakemake@input$scrna)
scRNA <- runTSNE(scRNA)
snRNA <- readRDS(snakemake@input$snrna)
snRNA <- runTSNE(snRNA)

scRNALung <- readRDS(snakemake@input$scRNALung)
scRNALung <- runTSNE(scRNALung)
liverAtlas <- readRDS(snakemake@input$liverAtlas)
liverAtlas <- runTSNE(liverAtlas)
tabulaVasc <- readRDS(snakemake@input@tabulaVasc)
tabulaVasc <- runTSNE(tabulaVasc)

plot_dim_red <- function(sce, mod, include_axis = FALSE,
                         a1_start = 0, a1_end = 0, a1_y = 0,
                         a2_start = 0, a2_end = 0, a2_x = 0, l_nrow = 4){
  if(!("CellType" %in% names(colData(sce)))){
    sce$CellType <- sce$cell_type
  }
  p <- plotTSNE(sce, colour_by = 'CellType', point_alpha = 0.3) +
    scale_color_manual(values = cell_type_colours(mod, FALSE)) +
    coord_fixed() +
    theme(legend.position = "none",
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.title = element_blank()) + 
    guides(color = guide_legend(nrow = l_nrow, byrow = TRUE))
  
  if(include_axis){
    p <- p + geom_segment(x = a1_start, xend = a1_end, y = a1_y, yend = a1_y,
                     arrow = arrow(length = unit(0.2, "cm")),
                     lineend = "round", linejoin = "mitre", size = 1) +
      geom_segment(x = a2_x, xend = a2_x, y = a2_start, yend = a2_end,
                   arrow = arrow(length = unit(0.2, "cm")),
                   lineend = "round", linejoin = "mitre", size = 1) +
      theme(axis.title.x = element_text(hjust = 0.035),
            axis.title.y = element_text(hjust = 0.035, vjust = -1.5))
  }else{
    p <- p +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank())
  }
  
  p
}


### CYTOF
cytof_tsne <- plot_dim_red(CyTOF, "CyTOF", TRUE,
             a1_start = -35, a1_end = -25, a1_y = -33,
             a2_start = -33, a2_end = -23, a2_x = -35)
cytof_bar <- CyTOF$cell_type |> 
  table() |> 
  as.data.frame() |> 
  ggplot(aes(x = reorder(Var1, -Freq), y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Freq, x = reorder(Var1, -Freq), y = Freq + 100,
                angle = 45, hjust = 0, vjust = 0),
            position = position_stack(vjust = 1.03)) +
  scale_fill_manual(values = cell_type_colours("CyTOF", FALSE)) +
  ylim(0, 2800) +
  labs(x = "Cell type", y = "Number of cells", fill = "Cell type") +
  whatsthatcell_theme() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

cytof_plot <- (wrap_elements(full = cytof_tsne, ignore_tag = TRUE) & 
                 labs(title = "CyTOF - Bone marrow")) /
  cytof_bar + plot_layout(heights = c(3,1.8))



### scRNASeq
scrna_tsne <- plot_dim_red(scRNA, "scRNASeq", TRUE,
                           a1_start = -40, a1_end = -30, a1_y = -49,
                           a2_start = -49, a2_end = -39, a2_x = -40, 
                           l_nrow = 5)

scrna_bar <- scRNA$CellType |> 
  table() |> 
  as.data.frame() |> 
  ggplot(aes(x = reorder(Var1, -Freq), y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Freq, x = reorder(Var1, -Freq), y = Freq + 50,
                angle = 45, hjust = 0, vjust = 0),
            position = position_stack(vjust = 1.01)) +
  scale_fill_manual(values = cell_type_colours("scRNASeq", FALSE)) +
  ylim(0, 2250) +
  labs(x = "Cell type", y = "Number of cells", fill = "Cell type") +
  whatsthatcell_theme() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


scrna_plot <- ((wrap_elements(full = scrna_tsne, ignore_tag = TRUE) & 
                  labs(title = "scRNASeq - Breast cancer cell lines")) /
  scrna_bar) + 
  plot_layout(heights = c(3,1))


### snRNASeq
snrna_tsne <- plot_dim_red(snRNA, "snRNASeq")

snrna_bar <- snRNA$cell_type |> 
  table() |> 
  as.data.frame() |> 
  ggplot(aes(x = reorder(Var1, -Freq), y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Freq, x = reorder(Var1, -Freq), y = Freq + 200,
                angle = 45, hjust = 0, vjust = 0),
            position = position_stack(vjust = 1.03)) +
  scale_fill_manual(values = cell_type_colours("snRNASeq", FALSE)) +
  ylim(0, 4200) +
  labs(x = "Cell type", y = "Number of cells", fill = "Cell type") +
  whatsthatcell_theme() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


snrna_plot <- ((wrap_elements(full = snrna_tsne, ignore_tag = TRUE) & 
                  labs(title = "snRNASeq - Pancreas cancer")) / 
                 snrna_bar) + 
  plot_layout(heights = c(3,1))

## scRNALung
scrna_lung_tsne <- plot_dim_red(scRNALung, "scRNALung")

scRNALung_bar <- scRNALung$CellType |> 
  table() |> 
  as.data.frame() |> 
  ggplot(aes(x = reorder(Var1, -Freq), y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Freq, x = reorder(Var1, -Freq), y = Freq + 50,
                angle = 45, hjust = 0, vjust = 0),
            position = position_stack(vjust = 1.03)) +
  scale_fill_manual(values = cell_type_colours("scRNALung", FALSE)) +
  ylim(0, 1800) +
  labs(x = "Cell type", y = "Number of cells", fill = "Cell type") +
  whatsthatcell_theme() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

scrna_lung_plot <- ((wrap_elements(full = scrna_lung_tsne, ignore_tag = TRUE) & 
                       labs(title = "scRNASeq - Lung cancer cell lines")) / 
                      scRNALung_bar) + 
  plot_layout(heights = c(3,1))

## Tabula Liver
liverAtlas_tsne <- plot_dim_red(liverAtlas, "liverAtlas")

liverAtlas_bar <- liverAtlas$CellType |> 
  table() |> 
  as.data.frame() |> 
  ggplot(aes(x = reorder(Var1, -Freq), y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Freq, x = reorder(Var1, -Freq), y = Freq + 40,
                angle = 45, hjust = 0, vjust = 0),
            position = position_stack(vjust = 1.03)) +
  scale_fill_manual(values = cell_type_colours("liverAtlas", FALSE)) +
  ylim(0, 4200) +
  labs(x = "Cell type", y = "Number of cells", fill = "Cell type") +
  whatsthatcell_theme() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

liverAtlas_plot <- ((wrap_elements(full = liverAtlas_tsne, ignore_tag = TRUE) & 
                       labs(title = "scRNASeq - Liver")) / 
                       liverAtlas_bar) + 
  plot_layout(heights = c(3,1.8))

## Tabula Vasc
tabulaVasc_tsne <- plot_dim_red(tabulaVasc, "tabulaVasc")

tabulaVasc_bar <- tabulaVasc$CellType |> 
  table() |> 
  as.data.frame() |> 
  ggplot(aes(x = reorder(Var1, -Freq), y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Freq, x = reorder(Var1, -Freq), y = Freq + 150, 
                angle = 45, hjust = 0, vjust = 0),
            position = position_stack(vjust = 1.03)) +
  scale_fill_manual(values = cell_type_colours("tabulaVasc", FALSE)) +
  ylim(0, 7000) +
  labs(x = "Cell type", y = "Number of cells", fill = "Cell type") +
  whatsthatcell_theme() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

tabulaVasc_plot <- ((wrap_elements(full = tabulaVasc_tsne, ignore_tag = TRUE) & 
                       labs(title = "scRNASeq - Vasculature")) / 
                      tabulaVasc_bar) + 
  plot_layout(heights = c(3,1.8))


## Combine figures
row1 <- wrap_elements((scrna_plot | scrna_lung_plot | snrna_plot) + plot_layout(widths = c(1, 1.2, 1)))
row2 <- wrap_elements((cytof_plot | liverAtlas_plot | tabulaVasc_plot) + plot_layout(widths = c(1, 1.2, 1)))

pdf(snakemake@output$fig1, height = 18, width = 17)
  wrap_elements(row1 / plot_spacer() / row2 + plot_layout(heights= c(1, 0.05, 1))) /
    theme(plot.tag = element_text(size = 22))
dev.off()

