#!/usr/bin/env Rscript

# Description:
#   This script compares the reactions present after the carving step
#   with the reactions present after the gap-filling step (SRP). It identifies
#   reactions retained after carving but absent after gap filling, and generates
#   a presence table indicating whether each reaction is specific to the
#   post-carving reaction set.
#
# Inputs:
#   - snakemake@input[["carveme"]]:
#       TSV file containing the reactions present after the carving step, where:
#         * Each row corresponds to one reaction ID.
#         * Only the first column is used.
#
#   - snakemake@input[["srp"]]:
#       TSV file containing the reactions present after the gap-filling step (SRP), where:
#         * Each row corresponds to one reaction ID.
#
# Outputs:
#   - snakemake@output[["tsv"]]:
#       TSV file containing the reaction presence table, where:
#         * The "reaction" column contains reaction ID from the
#           post-carving reaction set.
#         * The "gp_status" column indicates whether the reaction is absent
#           after the gap-filling step (SRP):
#             - "yes": reaction is present after carving but absent after gap filling.
#             - "no": reaction is present after both carving and gap filling.

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
