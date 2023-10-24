library(tidyverse)
source("pipeline/whatsthatcell-helpers.R")

AL_rem_files <- list.files("output/v1/AL_with_removed_input/", pattern = "removed", full.names = TRUE)
AL_kept_files <- list.files("output/v1/AL_with_removed_input/", pattern = "kept", full.names = TRUE)

AL_removed <- lapply(AL_rem_files, read_tsv) %>% 
  bind_rows()
AL_kept <- lapply(AL_kept_files, read_tsv) %>% 
  bind_rows()

#AL <- read_tsv("data/scRNASeq/Active-Learning/Active-Learning-annotation-GroundTruth.tsv")


# AL_removed %>%
#   group_by(seed, removed_markers, iteration, cell_type) %>%
#   tally() %>%
#   ggplot(aes(x = iteration, y = n, fill = cell_type)) +
#   geom_bar(stat = 'identity') +
#   labs(x = "Logistic regression training iteration", y = "Number of cells",
#        color = "Ground Truth cell type",
#        title = "Number of cells of each type that are selected when markers are removed") +
#   scale_fill_manual(values = cell_type_colours()) +
#   facet_grid(seed~removed_markers) +
#   whatsthatcell_theme()


prop_assigned_to_missing <- AL_removed %>%  # formerly a
  mutate(selected_of_removed_type = cell_type == removed_markers) %>% 
  group_by(seed, removed_markers, iteration, selected_of_removed_type) %>% 
  tally() %>% 
  mutate(total = sum(n)) %>% 
  mutate(prop = n / total) %>% 
  ungroup() %>% 
  group_by(removed_markers, iteration, selected_of_removed_type) %>% 
  summarize(mean_proportion = mean(prop))

prop_assigned_to_missing$selected_of_removed_type <- as.character(prop_assigned_to_missing$selected_of_removed_type)

prop_assigned_to_kept <- AL_kept %>%  # formerly a
  mutate(selected_of_removed_type = cell_type == removed_markers) %>% 
  group_by(seed, removed_markers, iteration, selected_of_removed_type) %>% 
  tally() %>% 
  mutate(total = sum(n)) %>% 
  mutate(prop = n / total) %>% 
  ungroup() %>% 
  group_by(removed_markers, iteration, selected_of_removed_type) %>% 
  summarize(mean_proportion = mean(prop))

prop_assigned_to_kept$selected_of_removed_type <- as.character(prop_assigned_to_kept$selected_of_removed_type)

# Create a function to get per cell type baseline proportions
baseline_assigned_prop <- function(df, marker){
  df <- df %>% 
    mutate(selected_of_removed_type = cell_type == marker) %>% 
    group_by(seed, iteration, selected_of_removed_type) %>% 
    tally() %>% 
    mutate(total = sum(n)) %>% 
    mutate(prop = n /total) %>% 
    ungroup() %>% 
    group_by(iteration, selected_of_removed_type) %>% 
    summarise(mean_proportion = mean(prop)) %>% 
    filter(selected_of_removed_type == TRUE)
  
  df$selected_of_removed_type <- "Baseline"
  df$removed_markers <- marker
  
  df[, c("removed_markers", "iteration", "selected_of_removed_type", "mean_proportion")]
}

prop_assigned_to_kept <- lapply(unique(AL_removed$removed_markers), function(x){
  baseline_assigned_prop(AL_kept, x)
}) %>% bind_rows()

AL %>%
  group_by(cell_type, iteration) %>% 
  tally() %>% 
  ggplot(aes(x = iteration, y = n, fill = cell_type)) +
  geom_bar(stat = 'identity') +
  labs(x = "Logistic regression training iteration", y = "Number of cells",
       color = "Ground Truth cell type",
       title = "Number of cells of each type that are selected when no markers are removed") +
  scale_fill_manual(values = cell_type_colours()) +
  whatsthatcell_theme()
  
bind_rows(prop_assigned_to_missing, prop_assigned_to_kept) %>% 
  ggplot(aes(x = iteration, y = mean_proportion, color = selected_of_removed_type)) +
  geom_line() +
  geom_point() +
  ylim(c(0,1)) +
  facet_wrap(~removed_markers) +
  labs(x = "Logistic regression training iteration", y = "Proportion of cells",
       color = "Is the selected cell\nof the removed type?",
       title = "Proportion of cells selected when markers are removed") +
  whatsthatcell_theme()



