# Description:
#   This script generates a two-dimensional scatter plot from UMAP coordinates
#   obtained from PCA-reduced flux sampling data. 
#
# Inputs:
#   - snakemake@input[["umap"]]:
#       TSV file containing UMAP coordinates, where:
#         * Rows correspond to flux samples (sampling solutions) from a single
#           metabolic model.
#         * Columns contain the two UMAP dimensions:
#             - "UMAP1"
#             - "UMAP2"
#
#   - snakemake@wildcards[["models"]]:
#       ID of the single metabolic model associated with the UMAP
#       coordinates.
#
# Outputs:
#   - snakemake@output[["png"]]:
#       PNG image containing the UMAP scatter plot, where:
#         * The X-axis represents the first UMAP dimension (UMAP1).
#         * The Y-axis represents the second UMAP dimension (UMAP2).
#         * Each point corresponds to one flux sample.

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

input_file  <- snakemake@input[["umap"]]
output_file <- snakemake@output[["png"]]

model_name <- snakemake@wildcards[["models"]]

umap_df <- fread(input_file)

p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(size = 1, alpha = 0.8) +
    theme_minimal() +
    ggtitle(paste0("UMAP on the first 30 PCs - model: ", model_name))

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

ggsave(output_file, plot = p, width = 6, height = 5, dpi = 300)