library(readr)
library(dplyr)
library(data.table)

input_file  <- snakemake@input[["fs"]]
output_file <- snakemake@output[["pca"]]

flux_table <- suppressMessages(
    read_tsv(file = input_file)
) |>
    dplyr::select(-c(1))

# transpose because:
# rows = features
# columns = samples

mat <- t(flux_table)

pca <- prcomp(
    mat,
    center = TRUE,
    scale. = TRUE
)

pca_df <- as.data.frame(pca$x)

dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
)

fwrite(
    pca_df,
    file = output_file,
    sep = "\t"
)