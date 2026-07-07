# Description:
#   This script generates a scatter plot from PCA coordinates
#   obtained using the pca_fs.R script. The plot displays the
#   first two principal components (PC1 and PC2) computed from flux sampling
#   solutions of a single metabolic model. The percentage of variance explained
#   by each principal component is displayed on the axes.
#
# Inputs:
#   - snakemake@input[["pca"]]:
#       TSV file containing PCA coordinates, where:
#         * Each row corresponds to one flux sampling solution.
#         * Columns correspond to principal component scores (PCs).
#         * The PCA coordinates are generated from a single metabolic model.
#
#   - snakemake@input[["var"]]:
#       TSV file containing PCA variance information, where:
#         * The "variance_explained" column contains the proportion of variance
#           explained by each principal component.
#         * Each row corresponds to one principal component.
#
#   - snakemake@wildcards[["models"]]:
#       ID of the metabolic model used to generate the flux sampling
#
# Outputs:
#   - snakemake@output[["png"]]:
#       PNG image containing the PCA scatter plot, where:
#         * The X-axis represents PC1.
#         * The Y-axis represents PC2.
#         * Each point corresponds to one flux sampling solution.

suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(data.table)
})



input_file <- snakemake@input[["pca"]]
var_file <- snakemake@input[["var"]]
output_file <- snakemake@output[["png"]]

model_name <- snakemake@wildcards[["models"]]

pca_df <- fread(input_file)
var_df <- fread(var_file)

if (!all(c("PC1", "PC2") %in% colnames(pca_df))) {
  stop("PC1 et PC2 introuvables dans le fichier PCA")
}

pc1_var <- var_df$variance_explained[1] * 100
pc2_var <- var_df$variance_explained[2] * 100

p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.7, size = 1) +
  theme_minimal() +
  labs(
    title = paste0("PCA (PC1 vs PC2) - ", model_name),
    x = paste0("PC1 (", round(pc1_var, 1), "%)"),
    y = paste0("PC2 (", round(pc2_var, 1), "%)")
  )

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = output_file,
  plot = p,
  width = 6,
  height = 5,
  dpi = 300
)