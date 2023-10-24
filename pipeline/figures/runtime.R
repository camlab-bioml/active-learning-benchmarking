suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})
source("pipeline/whatsthatcell-helpers.R")

runtimes <- list.files("output/v8/results/runtime/", full.names = TRUE)
runtimes <- lapply(runtimes, read_tsv) |> 
  bind_rows()

min_y <- min(runtimes$time)
max_y <- max(runtimes$time)

rf_runtime_p <- filter(runtimes, AL_alg == "rf") |> 
  ggplot(aes(x = method, y = time, fill = method)) +
  geom_boxplot() +
  ylim(min_y, max_y) +
  labs(x = "Selection procedure", y = "Runtime (s)", title = "Random forest AL") +
  scale_fill_manual(values = sel_met_cols) +
  facet_wrap(~ modality) + 
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none",
        axis.title.x = element_blank())

lr_runtime_p <- filter(runtimes, AL_alg == "multinom") |> 
  ggplot(aes(x = method, y = time, fill = method)) +
  geom_boxplot() +
  ylim(min_y, max_y) +
  labs(x = "Selection procedure", y = "Runtime (s)", title = "Logistic regression AL") +
  scale_fill_manual(values = sel_met_cols) +
  facet_wrap(~ modality) + 
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none",
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

ar_runtime_p <- filter(runtimes, is.na(AL_alg)) |> 
  ggplot(aes(x = method, y = time, fill = method)) +
  geom_boxplot() +
  ylim(min_y, max_y) +
  labs(x = "Selection procedure", y = "Runtime (s)", title = "Adaptive reweighting") +
  scale_fill_manual(values = sel_met_cols) +
  facet_wrap(~ modality) + 
  whatsthatcell_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none",
        axis.title = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

pdf("output/v8/paper-figures/runtime.pdf", height = 6, width = 15)
  (rf_runtime_p | lr_runtime_p | ar_runtime_p) +
    plot_layout(widths = c(1, 1, 0.4))
dev.off()


