suppressPackageStartupMessages({
    library(tidyverse)
})

count_cells <- function(df){
    df |>
        group_by(predicted_cell_type, prediction_params, selection_procedure) |>
        tally() |>
        ungroup() |>
        mutate(total_cells = sum(n))
}

total_counts <- lapply(snakemake@input$high_prob_cells, function(x){
    tsv <- read_tsv(x)
    counts <- count_cells(tsv)
}) |> bind_rows()

write_tsv(total_counts, snakemake@output$pred_lab_cell_nums)
