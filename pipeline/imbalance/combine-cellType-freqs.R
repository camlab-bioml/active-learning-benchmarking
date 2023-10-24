suppressPackageStartupMessages({
    library(tidyverse)
})

lapply(snakemake@input$tsvs, read_tsv) |>
    bind_rows() |> 
    separate(params, c('modality', 'bal', 'similarity', 's', 'seed')) |>
    select(-s) |>
    write_tsv(snakemake@output$tsv)
