suppressPackageStartupMessages({
    library(tidyverse)
})

acc <- lapply(snakemake@input$accs, read_tsv) |>
    bind_rows()

write_tsv(acc, snakemake@output$acc)