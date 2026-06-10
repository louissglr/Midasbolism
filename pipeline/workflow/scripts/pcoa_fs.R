library(readr)
library(dplyr)
library(data.table)
library(vegan)
library(ape)

input_file   <- snakemake@input[["fs"]]
sampled_file <- snakemake@input[["sampled"]]
output_file  <- snakemake@output[["pcoa"]]

# Load full flux table
flux_table <- suppressMessages(
    read_tsv(input_file)
)

# Load sampled IDs
sampled_ids <- fread(sampled_file)[[1]]

# Identify ID column
id_col <- names(flux_table)[1]

# Filter by sampled IDs 
flux_table <- flux_table[
    flux_table[[id_col]] %in% sampled_ids,
]

# Remove ID column
flux_values <- flux_table[, -1]

# Convert to matrix (solutions x features)
mat <- as.matrix(flux_values)

# Bray-Curtis distance
dist_mat <- vegdist(mat, method = "bray")

# PCoA
pcoa_res <- ape::pcoa(dist_mat)

# Keep first 2 axes
pcoa_df <- as.data.frame(pcoa_res$vectors[, 1:2])
colnames(pcoa_df) <- c("PCoA1", "PCoA2")


pcoa_df$solution_id <- flux_table[[id_col]]

dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
)

fwrite(
    pcoa_df,
    file = output_file,
    sep = "\t"
)