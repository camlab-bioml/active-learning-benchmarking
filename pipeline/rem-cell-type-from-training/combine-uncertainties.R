suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

lapply(snakemake@input$tsv, read_tsv) |>
    bind_rows() |>
    write_tsv(snakemake@output$tsv)

