### [ WHATSTHATCELL FUNCTIONS ] #####

### [ ACTIVE LEARNING ] #####
select_initial_cells <- function(df_expression,
                                 marker_dict,
                                 number_cells = 20) {
  
  cell_types <- names(marker_dict)
  
  selected_cells <- lapply(cell_types, function(x) {
    #Negative markers
    neg_markers <- marker_dict[[x]]$negative
    if (!is.null(neg_markers)) {
      neg_marker_expr <- df_expression[, neg_markers]
      neg_avg_expression <-
        apply(neg_marker_expr, 1, calculate_expression_average)
    } else{
      neg_avg_expression <- 0
    }
    #Positive markers
    pos_markers <- marker_dict[[x]]$positive
    if (!is.null(pos_markers)) {
      pos_marker_expr <- df_expression[, pos_markers]
      pos_avg_expression <-
        apply(pos_marker_expr, 1, calculate_expression_average)
    } else{
      pos_avg_expression <- 0
    }
    
    average_expression <- pos_avg_expression - neg_avg_expression
    df <- tibble(marker = average_expression)
    colnames(df) <- names(marker_dict[x])
    df
  }) %>% bind_cols()
  
  selected_cells <- selected_cells |> 
    rowwise() |> 
    mutate(max_mean_expression = max(across()),
           max_cell_type = cell_types[which.max(c_across(everything(selected_cells)))]) |> 
    ungroup() |> 
    mutate(X1 = df_expression$X1) |> 
    select(X1, max_cell_type, max_mean_expression) |> 
    arrange(-max_mean_expression)
    
  cell_types_to_sample <- recycle(sample(cell_types), nrow(selected_cells))
  
  lapply(cell_types_to_sample, function(i){
    cell <- filter(selected_cells, max_cell_type == i) |> 
      slice_head(n = 1)
    
    selected_cells <<- filter(selected_cells, !(X1 %in% cell$X1))
    cell
  }) |> bind_rows() |> 
    slice_head(n = number_cells) |> 
    pull(X1)
}

calculate_expression_average <- function(expressions) {
  return(sum(expressions) / length(expressions))
}

load_scs <- function(scs_expression) {
  df_expression <- as.data.frame(as.matrix(t(assays(scs_expression)$logcounts)))
  setDT(df_expression, keep.rownames = "X1")[]
  df_expression <- as_tibble(df_expression)
  return(df_expression)
}

calculate_entropy <- function(predictions) {
  entropy <- -(lapply(predictions, function(x) {
    if (x == 0) {
      return(0)
    } else {
      return(x * log(x, 2))
    }
  }) %>% unlist() %>% sum())
  return(entropy)
}

# entropy based selection of cells
select_entropy <- function(entropy_df, method, amount, random_selection){
  # Determine how many cells should be selected non-randomly
  select_amount <- round(amount * (1-random_selection), 0)
  
  # Determine how many cells should be selected randomly
  random_amount <- amount - select_amount
  
  # Select the cells with the highest entropy
  if(method == "highest_entropy"){
    cells <- c(entropy_df[order(entropy_df$entropy, 
                                decreasing = T),][1:select_amount, ]$X1)
  }
  
  # Find appropriate quantile threshold
  else if(method == "0.95_quant_entropy"){
    quantile <- quantile(entropy_df$entropy, 0.95)
  }else if(method == "0.75_quant_entropy"){
    quantile <- quantile(entropy_df$entropy, 0.75)
  }
  
  # Filter cells based on quantiles
  if(grepl("_quant_entropy", method)){
    ordered_cells <- select(entropy_df, X1, entropy) %>% 
      arrange(entropy) %>% 
      filter(entropy >= quantile)
    
    if(nrow(ordered_cells) > amount){
      cells <- ordered_cells$X1[1:select_amount]
    }else{
      cells <- ordered_cells$X1
    }
  }
  
  # Randomly select additional cells
  if(random_selection != 0){
    if(nrow(entropy_df) > random_amount){
      random_cells <- sample(entropy_df$X1, random_amount)
    }else{
      random_cells <- entropy_df$X1
    }
    cells <- c(cells, random_cells)
  }
  
  cells
}

# maximum probability based selection of cells
select_maxp <- function(maxp_df, method = "lowest_maxp", amount, random_selection){
  # Determine how many cells should be selected non-randomly
  select_amount <- round(amount * (1-random_selection), 0)
  
  # Determine how many cells should be selected randomly
  random_amount <- amount - select_amount
  
  # Select the cells with the highest maxp
  if(method == "lowest_maxp"){
    cells <- maxp_df[order(maxp_df$maxp, decreasing = F),][1:select_amount, ]$X1
  }
  
  # Find appropriate quantile threshold
  else if(method == "0.05_quant_maxp"){
    quantile <- quantile(maxp_df$maxp, 0.05)
  }else if(method == "0.25_quant_maxp"){
    quantile <- quantile(maxp_df$maxp, 0.25)
  }
  
  # Filter cells based on quantiles
  if(grepl("_quant_maxp", method)){
    ordered_cells <- select(maxp_df, X1, maxp) %>%
      arrange(maxp) %>%
      filter(maxp >= quantile)

    if(nrow(ordered_cells) > amount){
      cells <- ordered_cells$X1[1:select_amount]
    }else{
      cells <- ordered_cells$X1
    }
  }
  
  # Randomly select additional cells
  if(random_selection != 0){
    if(nrow(maxp_df) > random_amount){
      random_cells <- sample(maxp_df$X1, random_amount)
    }else{
      random_cells <- maxp_df$X1
    }
    cells <- c(cells, random_cells)
  }
  
  cells
}

fit_AL_classifier <- function(annotated_cells, AL_method){
  training_set_size <- nrow(annotated_cells)
  if(training_set_size < 50){
    bootstrap_reps <- 1000
  }else if(training_set_size >=50 | training_set_size < 100){
    bootstrap_reps <- 125
  }else{
    bootstrap_reps <- 100
  }
  tctrl <- trainControl(method = "boot", number = bootstrap_reps)
  ModelFit <- train(cell_type ~ ., 
                    data = select(annotated_cells, -X1),
                    method = AL_method,
                    trace = FALSE,
                    trControl = tctrl)

  ModelFit
}

entropy_maxp_cell_selection <- function(selection_criterion, predicted_scores, left_cells, selection_method, 
                                        amount, random_selection, annotated_cells){
  if(selection_criterion == "entropy"){
    entropies <- apply(predicted_scores, 1, calculate_entropy)
    
    left_cells$entropy <- entropies
    
    selected_cells <- select_entropy(entropy_df = left_cells, 
                                     method = selection_method, 
                                     amount = amount, 
                                     random_selection = random_selection)
    
    criterion_table <- data.frame(cell_id = left_cells$X1, 
                                  criterion_val = left_cells$entropy, 
                                  no_cells_annotated = nrow(annotated_cells),
                                  criterion = selection_criterion)
    
  }else if(selection_criterion == "maxp"){
    max_p_idx <- apply(predicted_scores, 1, which.max)
    max_p <- lapply(seq_len(nrow(predicted_scores)), function(x){
      predicted_scores[x, max_p_idx[x]]
    }) %>% unlist()
    
    left_cells$maxp <- max_p
    
    selected_cells <- select_maxp(maxp_df = left_cells, 
                                  method = selection_method, 
                                  amount = amount, 
                                  random_selection = random_selection)
    
    criterion_table <- data.frame(cell_id = left_cells$X1, 
                                  criterion_val = left_cells$maxp, 
                                  no_cells_annotated = nrow(annotated_cells),
                                  criterion = selection_criterion)
  }
  
  list(selected_cells = selected_cells, criterion_table = criterion_table)
}

select_cells_classifier <- function(df_expression, AL_method, selection_method, amount = 10, 
                                    random_selection, selection_criterion = "entropy") {
  annotated_cells <- df_expression %>% 
    filter(!is.na(cell_type)) %>% 
    filter(cell_type != "Skipped", cell_type != "Unclear")
  
  left_cells <- df_expression %>% 
    filter(is.na(cell_type))
  
  ModelFit <- fit_AL_classifier(annotated_cells, AL_method)
  
  predicted_scores <- predict(ModelFit, 
                              select(left_cells, -X1, -cell_type),
                              type = "prob")
  
  sel_cells <- entropy_maxp_cell_selection(selection_criterion, predicted_scores, left_cells, selection_method,
                                           amount, random_selection, annotated_cells)
  
  
  return_list <- list(selected_cells = sel_cells$selected_cells, 
                      criterion_table = sel_cells$criterion_table)
  return(return_list)
}


cell_ranking_wrapper <- function(df, markers, number_cells = 20){
  # Get initial set of cells based on their marker expression ranking
  ranked_cells <- select_initial_cells(df, markers$cell_types, number_cells)
  
  # What index do the selected cells correspond to?
  to_assign_index <- match(ranked_cells, df$X1)
  df$cell_type[to_assign_index] <- df$gt_cell_type[to_assign_index]
  df$iteration[to_assign_index] <- 0
  
  df
}


active_learning_wrapper <- function(df, AL_method, selection_method, iteration, 
                                    entropies = NULL, random_selection, selection_criterion = "entropy", num_cells = 10){
  # AL selected cells
  AL <- select_cells_classifier(df_expression = df, 
                                AL_method = AL_method, 
                                selection_method = selection_method,
                                amount = num_cells, 
                                random_selection = random_selection,
                                selection_criterion = selection_criterion)
  
  list(new_cells = AL$selected_cells, criterion_table = AL$criterion_table)
}

corrupt_labels <- function(df_expression, all_cell_types, prop_corrupt){
  num_corrupt <- round((nrow(df_expression) * as.numeric(prop_corrupt)), 0)
  corrupt_idx <- sample(1:nrow(df_expression), num_corrupt)

  for(i in 1:length(all_cell_types)){
    tc <- all_cell_types[i]

    to_corrupt <- which(df_expression$gt_cell_type == tc)
    to_corrupt <- intersect(to_corrupt, corrupt_idx)
    if(length(to_corrupt) > 0){
      subset_cell_types <- all_cell_types[!(tc == all_cell_types)]
      df_expression$corrupted[to_corrupt] <- sample(subset_cell_types, 
                                                    length(to_corrupt), 
                                                    replace = TRUE)
    }
  }
  
  df_expression$gt_cell_type[!is.na(df_expression$corrupted)] <- 
    df_expression$corrupted[!is.na(df_expression$corrupted)]
  
  df_expression
}


### [ REMOVING CELL TYPES ] #####
# Used to read in a sce and convert to PCA/expression
create_features <- function(sce){
  if(any(grepl("CD16_32", rownames(sce)))){
    rownames(sce)[grep("CD16_32", rownames(sce))] <- "CD16-32"
  }
  
  df_expression <- load_scs(sce)
  df_expression$cell_type <- NA
  df_expression$iteration <- NA
  df_expression$gt_cell_type <- sce$CellType
  
  # Remove genes with 0 expression
  expression <- df_expression[, c(TRUE, 
                                     colSums(df_expression[,2:(ncol(df_expression)-3)]) > 0,
                                     rep(TRUE, 3))]
  
  df_PCA <- select(expression, -c(X1, cell_type, iteration, gt_cell_type)) |> 
    as.matrix() |> 
    prcomp(center = TRUE, scale. = TRUE)
  
  df_PCA <- df_PCA$x |> 
    as.data.frame()
  
  df_PCA <- bind_cols(
    tibble(X1 = df_expression$X1),
    df_PCA[,1:min(20, ncol(df_PCA))], 
    tibble(gt_cell_type = df_expression$gt_cell_type, # This is where the ground truth labels are stored
           iteration = df_expression$iteration)
  )
  
  list(expression = df_expression, PCA = df_PCA)
}

rem_cell_type_AL_wrapper <- function(df, AL_alg, strat, rand, criterion, iter = 3, num_cells){
  set.seed(1)
  entropies <- list()
  for(i in 1:iter){
    AL <- active_learning_wrapper(select(df, -gt_cell_type, -iteration), 
                                  AL_alg,
                                  strat,
                                  i, 
                                  entropies, 
                                  rand,
                                  criterion)
    
    if(criterion == 'entropy'){
      num_cell_types <- unique(na.omit(df$cell_type)) |> 
        length()
      AL$criterion_table$criterion_val <- AL$criterion_table$criterion_val / log(num_cell_types, 2)
    }
    entropies[[length(entropies) + 1]] <- AL$criterion_table
    
    # What index do the selected cells correspond to?
    to_assign_index <- match(AL$new_cells, df$X1)
    
    # Get ground truth labels based on the index
    df$cell_type[to_assign_index] <- df$gt_cell_type[to_assign_index]
    df$iteration[to_assign_index] <- i
  }
  entropies
}

get_training_type_rem <- function(df, initial, markers, cell_type_to_rem, needed_cells = 20){
  selected_cells <- 0
  discarded_cells <- c()
  
  selected_cells_df <- tibble(
    X1 = character(),
    cell_type = character()
  )
  
  # This loop iteratively runs until the determined size of the training set has been found
  # Cells that are of the cell type that is not selected are then removed
  while(selected_cells < needed_cells){
    left_to_rank_num <- needed_cells - selected_cells
    
    if(initial == "ranking"){
      df <- cell_ranking_wrapper(df, markers, left_to_rank_num) |>
        filter(!(cell_type %in% cell_type_to_rem)  | is.na(cell_type))
    }else if(initial == "random"){
      random_cell_idx <- sample(1:nrow(df), left_to_rank_num)
      df$cell_type[random_cell_idx] <- df$gt_cell_type[random_cell_idx]
      
      df <- filter(df, !(cell_type %in% cell_type_to_rem)  | is.na(cell_type))
    }
    
    # Save selected cells minus cells to remove
    selected_cells_df <- bind_rows(
      selected_cells_df,
      filter(df, !is.na(cell_type)) |> 
        select(X1, cell_type)
    )
    
    # Remove already selected cells so that new cells are selected with new ranking
    df <- filter(df, is.na(cell_type))
    
    # How many cells have been selected total?
    selected_cells <- nrow(selected_cells_df)
  }
  selected_cells_df
}

get_training_type_kept <- function(df, cell_type_to_rem, initial, markers,
                                   selected_cells_df){
  # Find the set of cells of the type that were removed above with the highest marker expression
  missing_df_expression <- filter(df, gt_cell_type %in% cell_type_to_rem)

  if(nrow(missing_df_expression) < 3){
    ranked_cell_num <- nrow(missing_df_expression)
  }else{
    ranked_cell_num <- 3
  }
  
  if(initial == "ranking"){
    # Which 3 cells of this type are selected by ranking?
    missing_df_expression <- cell_ranking_wrapper(missing_df_expression, 
                                                  markers, ranked_cell_num) |> 
      filter(!is.na(cell_type)) |> 
      select(X1, cell_type)
  }else if(initial == "random"){
    random_cell_idx <- sample(1:nrow(missing_df_expression), ranked_cell_num)
    missing_df_expression$cell_type[random_cell_idx] <- missing_df_expression$gt_cell_type[random_cell_idx]
    
    missing_df_expression <- filter(missing_df_expression, !is.na(cell_type)) |> 
      select(X1, cell_type)
  }
  
  # Create three datasets with 1, 2 and 3 cells of the removed type
  kept_cell_type <- lapply(1:ranked_cell_num, function(x){
    bind_rows(
      slice_head(selected_cells_df, n = (nrow(selected_cells_df) - x)), # original dataset with cell type removed
      slice_head(missing_df_expression, n = x) # Dataset created with only said cell type
    ) |> mutate(num_missing_cells = x)
  })
  
  kept_cell_type
}


### [ ACCURACIES ] #####
acc_wrap <- function(tt) {
  cell_types <- unique(union(tt$predicted_cell_type, tt$annotated_cell_type))
  
  tt$annotated_cell_type <- factor(tt$annotated_cell_type, levels = cell_types)
  tt$predicted_cell_type <- factor(tt$predicted_cell_type, levels = cell_types)
  
  bind_rows(
    tryCatch({kap(tt, annotated_cell_type, predicted_cell_type)}, error=function(e) NULL),
    tryCatch({bal_accuracy(tt, annotated_cell_type, predicted_cell_type)}, error=function(e) NULL),
    tryCatch({mcc(tt, annotated_cell_type, predicted_cell_type)}, error=function(e) NULL),
    tryCatch({sensitivity(tt, annotated_cell_type, predicted_cell_type)}, error=function(e) NULL),
    tryCatch({f_meas(tt, annotated_cell_type, predicted_cell_type)}, error=function(e) NULL),
  )
}


# ordering of methods for benchmarking
get_ranked_mean_estimate_by_cohort <- function(acc){
  lapply(unique(acc$cohort), function(x){
    mean_estimate <- filter(acc, cohort == x) |> 
      group_by(selection_procedure) |> 
      summarise(mean_estimate = mean(na.omit(.estimate))) |> 
      arrange(-mean_estimate)
    
    mean_estimate$rank <- 1:nrow(mean_estimate)
    select(mean_estimate, -mean_estimate) |> 
      mutate(cohort = x)
  }) |> 
    bind_rows()
}

get_cross_cohort_mean_ranked_estimate <- function(acc){
  acc |> 
    group_by(selection_procedure) |> 
    summarise(mean_rank = mean(rank)) |> 
    arrange(-mean_rank)
}


### Functions
al_selection <- function(acc, metric, active_learner, initial_sel, title, n_cells = NULL){
  sel_acc <- filter(acc, corrupted == 0) |> 
    filter(AL_alg == active_learner | is.na(AL_alg)) |> 
    filter(rand == 0 | is.na(rand)) |> 
    filter(.metric == metric) |> 
    filter(initial == initial_sel | is.na(initial))
  
  cross_cohort_ranks <- get_ranked_mean_estimate_by_cohort(sel_acc)
  method_order <- get_cross_cohort_mean_ranked_estimate(cross_cohort_ranks) |> 
    pull(selection_procedure)
  
  # p <- sel_acc |> 
  #   group_by(selection_procedure, cell_num, cohort) |> 
  #   mutate(selection_procedure = factor(selection_procedure, levels = method_order),
  #          mean_estimate = mean(.estimate),
  #          lower_quant = quantile(na.omit(.estimate), 0.25),
  #          upper_quant = quantile(na.omit(.estimate), 0.75)) |> 
  #   ggplot() +
  #   #aes(x = as.character(cell_num), y = .estimate, fill = selection_procedure) +
  #   #geom_boxplot() +
  #   geom_segment(aes(x = as.character(cell_num), xend = as.character(cell_num),
  #                    y = lower_quant, yend = upper_quant)) +
  #   geom_point(aes(x = as.character(cell_num), y = .estimate, colour = selection_procedure)) +
  #   scale_fill_manual(values = sel_meth_cols) +
  #   whatsthatcell_theme() +
  #   facet_wrap(~cohort, scales = "free_y", nrow = 1) +
  #   labs(title = metric) +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  p <- sel_acc |> 
    group_by(selection_procedure, cell_num, cohort) |> 
    mutate(selection_procedure = factor(selection_procedure, levels = method_order),
           median_estimate = median(na.omit(.estimate)),
           lower_quant = quantile(na.omit(.estimate), 0.25),
           upper_quant = quantile(na.omit(.estimate), 0.75)) |> 
    ungroup() |> 
    ggplot() +
    aes() +
    geom_segment(aes(x = as.character(cell_num), group = selection_procedure,
                     y = lower_quant, yend = upper_quant),
                 position = position_dodge(width = 0.9)) +
    geom_point(aes(x = as.character(cell_num), y = median_estimate, colour = selection_procedure),
               position = position_dodge(width = 0.9)) +
    scale_color_manual(values = sel_meth_cols) +
    whatsthatcell_theme() +
    facet_wrap(~cohort, scales = "free_y", nrow = 1) +
    labs(title = title, y = ".estimate", colour = "Selection method") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  if(!is.null(n_cells)){
    p <- p + labs(x = n_cells)
  }else{
    p <- p +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank())
  }
  
  p
}

full_acc_plot_wrapper <- function(acc, AL_alg, initial_sel, title){
  bal_acc <- al_selection(acc, "bal_accuracy", AL_alg, initial_sel, "Balanced accuracy")
  f1 <- al_selection(acc, "f_meas", AL_alg, initial_sel, "F1-score")
  kap <- al_selection(acc, "kap", AL_alg, initial_sel, "Kappa")
  mcc <- al_selection(acc, "mcc", AL_alg, initial_sel, "Matthews correlation coefficient") 
  sens <- al_selection(acc, "sensitivity", AL_alg, initial_sel, "Sensitivity", n_cells = "Number of cells")
  
  ((bal_acc + theme(legend.position = "none")) / 
  (f1 + theme(legend.position = "none")) / 
  (kap + theme(legend.position = "none")) /
  (mcc + theme(legend.position = "none")) / 
  (sens + theme(legend.position = "bottom")))
}

### [ PLOTTING ] ####
createHeatmap <- function(sce,
                          cell_type_column = "cell_type",
                          assay = "logcounts",
                          thresh = 2,
                          title = NULL,
                          modality) {
  
  lc <- t(as.matrix(assay(sce, assay)))
  lc <- scale(lc)
  keep_cols <- apply(lc, 2, function(x) all(is.na(x) == FALSE)) %>% as.vector
  lc <- lc[, keep_cols]
  
  lc[lc > thresh] <- thresh
  lc[lc < -thresh] <- -thresh
  
  cell_types = colData(sce)[[ cell_type_column ]]
  
  celltype_annot <- HeatmapAnnotation(`Cell type` = cell_types, 
                                      which="column",
                                      col = list(`Cell type` = cell_type_colours(modality)))  
  
  if(is.null(title)){
    title <- title
  }else{
    title <- paste(title, "\nCell")
  }
  
  type_exprs <- Heatmap(t(lc), 
                        name = "Expression",
                        column_title = title,
                        col=viridis(100),
                        top_annotation = celltype_annot,
                        show_column_names = FALSE,
                        column_order = order(cell_types))
  type_exprs
}

create_entropy_boxplot <- function(vals, y_lab, colors){
  box <- vals %>% 
    mutate(no_cells_annotated = factor(no_cells_annotated)) %>% 
    ggplot(aes(x=no_cells_annotated, y=criterion_val, 
               group = no_cells_annotated, fill = no_cells_annotated)) + 
    geom_boxplot() +
    theme_bw() +
    scale_fill_manual(values = colors) +
    labs(y = y_lab) +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())
  
  box
}



cell_type_colours <- function(modality, include_unassigned = TRUE) {
  pal <- c("#8B5B42", "#AF4EA9", "#FFB60A", "#79E200", "#0024DD", "#6CC1FF",
           "#0496FF", "#1DA05B", "#E11E00", "#A78882", "#BD93D8", "#fff53d", 
           '#FD4FBD', '#9A031E', '#E36414', '#6FBAAC', '#0F4C5C', '#4D76B8',
           '#98473E', '#97DFFC')
  
  scRNASeq_colours <- c(
    "CAL51" = pal[2],
    "BT483" = pal[3],
    "MCF7" = pal[12],
    "JIMT1" = pal[4],
    "AU565" = pal[5],
    "HCC1937" = pal[8], 
    "HDQP1" = pal[9],
    "CAL851" = pal[13],
    "CAMA1" = pal[10],
    "MDAMB468" = pal[6]
  )
  snRNASeq_colours <- c(
    "Fibroblast" = pal[1],
    "Tumor" = pal[2],
    "Endothelial" = pal[9],
    "Immune" = pal[3],
    "Ductal" = pal[5],
    "SmoothMuscle" = pal[10],
    "Schwann" = pal[4],
    "Atypical_Ductal" = pal[7],
    "Endocrine" = pal[8],
    "Acinar" = pal[6]
  )

  CyTOF_colours <- c(
    "Eosinophils" = pal[9],
    "pDCs" = pal[20],
    "NKT cells" = pal[15],
    "Macrophages" = pal[17],
    "Classical Monocytes" = pal[4],
    "Intermediate Monocytes" = pal[16],
    "Non-Classical Monocytes" = pal[8],
    "MEP" = pal[5],
    "CMP" = pal[6],
    "GMP" = pal[7],
    "MPP" = pal[18],
    "Basophils" = pal[19],
    "NK cells" = pal[10],
    "CD4 T cells" = pal[12],
    "CD8 T cells" = pal[3],
    "IgDpos IgMpos B cells" = pal[14],
    "IgD- IgMpos B cells" = pal[11],
    "IgM- IgD- B-cells" = pal[13],
    "B-cell Frac A-C (pro-B cells)" = pal[2]
  )
  
  scRNALung_colours <- c(
    "A549" = pal[1],
    "H1975" = pal[2],
    "H2228" = pal[7],
    "H838" = pal[8],
    "HCC827" = pal[12]
  )
  
  liverAtlas_colours <- c(
    "Endothelial cells" = pal[9],
    #"endothelial cell of hepatic sinusoid" = pal[13],
    "Basophils" = pal[14],
    "Fibroblasts" = pal[19],
    "Hepatocytes" = pal[3],
    "Cholangiocytes" = pal[18],
    "Dendritic cells" = pal[4],
    "Macrophages" = pal[8],
    "Mono+mono derived cells" = pal[16],
    "Neutrophils" = pal[13],
    "Resident NK" = pal[10],
    "Plasma cells" = pal[2],
    "T cells" = pal[12],
    "B cells" = pal[11],
    "pDCs" = pal[20]
  )
  
  tabulaVasc_colours <- c(
    "artery endothelial cell" = pal[9],
    "b cell" = pal[2],
    "endothelial cell" = pal[13],
    "epithelial cell" = pal[18],
    "erythrocyte" = pal[14],
    "fibroblast" = pal[19],
    "lymphatic endothelial cell" = pal[15],
    "macrophage" = pal[8],
    "mast cell" = pal[11],
    "nk cell" = pal[10],
    "pericyte cell" = pal[3],
    "plasma cell" = pal[20],
    "smooth muscle cell" = pal[5],
    "t cell" = pal[12]
  )

  scRNALung_colours <- c(
    "A549" = pal[1],
    "H1975" = pal[2],
    "H2228" = pal[7],
    "H838" = pal[8],
    "HCC827" = pal[12]
  )

  if(include_unassigned){
    scRNASeq_colours <- c(scRNASeq_colours, c('unassigned' = "grey60"))
    snRNASeq_colours <- c(snRNASeq_colours, c('unassigned' = "grey60"))
    CyTOF_colours <- c(CyTOF_colours, c('unassigned' = "grey60"))
    scRNALung_colours <- c(scRNALung_colours, c('unassigned' = "grey60"))
    liverAtlas_colours <- c(liverAtlas_colours, c('unassigned' = "grey60"))
    tabulaVasc_colours <- c(tabulaVasc_colours, c('unassigned' = "grey60"))
  }
  
  if(modality == "scRNASeq"){
    scRNASeq_colours
  }else if(modality == "snRNASeq"){
    snRNASeq_colours
  }else if(modality == "CyTOF"){
    CyTOF_colours
  }else if(modality == "scRNALung"){
    scRNALung_colours
  }else if(modality == "liverAtlas"){
    liverAtlas_colours
  }else if(modality == "tabulaVasc"){
    tabulaVasc_colours
  }
}

al_colours <- function(){
  c("LR" = "#DA94D4", "RF" = "#7EA3CC")
}

sel_met_cols <- c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F", "#00B9E3", "#619CFF", "#DB72FB", "#FF61C3") #hue_pal()(9)
names(sel_met_cols) <- c("AR Marker", "AR No Marker", "Random", "AL 0.05-maxp", "AL 0.25-maxp",
                         "AL 0.75-entropy", "AL 0.95-entropy", "AL Highest entropy", "AL Lowest maxp")

whatsthatcell_theme <- function(){
  theme_bw() +
    theme(strip.background=element_rect(fill="white"))
}


plot_top_markers <- function(sce, clusters, n=20) {
  fm <- findMarkers(sce, clusters)
  
  top_markers_per_cluster <- lapply(fm, function(f) {
    rownames(head(f, n = n))
  })
  
  top_markers <- unique(unlist(top_markers_per_cluster))
  
  expr_mat <- sce[top_markers,] %>% 
    logcounts() %>% 
    t() %>% 
    scale() 
  
  thresh <- 2
  expr_mat[expr_mat > thresh] <- thresh
  expr_mat[expr_mat < -thresh] <- -thresh
  
  dfc <- expr_mat %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    mutate(cluster = clusters) %>% 
    gather(gene, expression, -cluster)
  
  dfc2 <- group_by(dfc, cluster, gene) %>% 
    summarise(mean_expression = mean(expression), 
              pct_expressing = 100 * mean(expression > 0)) %>% 
    ungroup()
  
  
  dfc2$gene <- sub("(.*?)[^-]-", "", dfc2$gene)
  
  ggplot(dfc2, aes(y = cluster, x = gene)) +
    geom_point(aes(colour = mean_expression, size = pct_expressing)) +
    scale_colour_viridis_c(name = "Expression") +
    scale_size(name = "%  cells expressing") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "top") +
    labs(x = "Gene", y = "Cluster")
}
