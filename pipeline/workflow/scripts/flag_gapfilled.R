#!/usr/bin/env Rscript

library(tidyverse)

#Input / Output
carveme_file <- snakemake@input[["carveme"]]
srp_file     <- snakemake@input[["srp"]]

output_file  <- snakemake@output[["tsv"]]

# Load data

model_reactions <- read_tsv(
  carveme_file,
  col_names = FALSE,
  show_col_types = FALSE
) %>%
  pull(1) %>%
  unique()

srp_reactions <- read_tsv(
  srp_file,
  col_names = FALSE,
  show_col_types = FALSE
) %>%
  pull(1) %>%
  unique()

# Build presence table

df <- tibble(
  reaction = model_reactions,
  gp_status = if_else(
    !(model_reactions %in% srp_reactions),
    "yes",
    "no"
  )
)

#Save TSV
write_tsv(df, output_file)
