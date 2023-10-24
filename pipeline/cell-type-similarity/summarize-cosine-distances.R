library(tidyverse)

dist <- lapply(snakemake@input$tsvs, read_tsv) |>
    bind_rows()

dist |> 
    filter(cohort == "scRNASeq") |>
    group_by(cell_type1, cell_type2) |>
    summarise(mean_cosine_similarity = mean(cosine_similarity)) |>
    write_tsv(snakemake@output$avg_dist_sc)

dist |> 
    filter(cohort == "snRNASeq") |>
    group_by(cell_type1, cell_type2) |>
    summarise(mean_cosine_similarity = mean(cosine_similarity)) |>
    write_tsv(snakemake@output$avg_dist_sn)

dist |> 
    filter(cohort == "CyTOF") |>
    group_by(cell_type1, cell_type2) |>
    summarise(mean_cosine_similarity = mean(cosine_similarity)) |>
    write_tsv(snakemake@output$avg_dist_cy)

dist |> 
    filter(cohort == "scRNALung") |>
    group_by(cell_type1, cell_type2) |>
    summarise(mean_cosine_similarity = mean(cosine_similarity)) |>
    write_tsv(snakemake@output$avg_dist_scLung)
