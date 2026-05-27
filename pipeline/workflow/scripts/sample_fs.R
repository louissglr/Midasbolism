library(readr)
library(dplyr)
library(data.table)

input_file  <- snakemake@input[["fs"]]
output_file <- snakemake@output[["sampled"]]

n_samples <- snakemake@params[["n"]]
seed <- snakemake@params[["seed"]]

set.seed(seed)

flux_table <- suppressMessages(
    read_tsv(input_file)
)

# Première colonne = IDs des solutions
solution_ids <- flux_table[[1]]

# Sampling reproductible des IDs
sampled_ids <- sample(
    solution_ids,
    size = min(n_samples, length(solution_ids)),
    replace = FALSE
)

sampled_df <- data.frame(
    solution_id = sampled_ids
)

dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
)

fwrite(
    sampled_df,
    file = output_file,
    sep = "\t"
)