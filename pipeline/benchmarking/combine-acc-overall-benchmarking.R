suppressPackageStartupMessages({
    library(tidyverse)
})

accs <- lapply(snakemake@input$tsvs, read_tsv) |>
    bind_rows()

write_tsv(accs, snakemake@output$acc)