#!/usr/bin/env Rscript

# Description:
#   This script compares the reaction content of a metabolic model generated
#   with CarveMe against a reference reaction set (SRP). It identifies
#   reactions present in the metabolic model but absent from the reference
#   set and generates a presence table indicating whether each reaction is
#   specific to the model.
#
# Inputs:
#   - snakemake@input[["carveme"]]:
#       TSV file containing the list of reactions from CarveMe, where:
#         * Each row corresponds to one reaction identifier.
#         * Only the first column is used.
#
#   - snakemake@input[["srp"]]:
#       TSV file containing the reference reaction set, where:
#         * Each row corresponds to one reaction identifier.
#         * Only the first column is used.
#
# Outputs:
#   - snakemake@output[["tsv"]]:
#       TSV file containing the reaction presence table, where:
#         * The "reaction" column contains reaction identifiers from the
#           metabolic model.
#         * The "gp_status" column indicates whether the reaction is absent
#           from the SRP reference set:
#             - "yes": reaction is present in the model but absent from SRP.
#             - "no": reaction is also present in the SRP set.

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
