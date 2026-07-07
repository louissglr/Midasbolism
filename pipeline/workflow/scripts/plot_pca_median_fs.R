# Description:
#   This script generates a scatter plot from PCA coordinates
#   obtained using the pca_median_fs.R script. The plot displays the
#   first two principal components (PC1 and PC2) and labels each point with
#   its corresponding metabolic model. The percentage of variance explained by
#   each principal component is displayed on the axes.
#
# Inputs:
#   - snakemake@input[["pca"]]:
#       TSV file containing PCA coordinates where:
#         * Each row corresponds to a metabolic model.
#         * Columns correspond to principal component scores (PCs).
#         * The first column "model" contains the model ID.
#
#   - snakemake@input[["var"]]:
#       TSV file containing PCA variance information, where:
#         * Contains the variance explained by each PC (col1 = variance_explained).
#         * Contains the cumulative variance explained across PCs (col2 = cumulative_variance).
#
# Outputs:
#   - snakemake@output[["png"]]:
#       PNG image containing the PCA scatter plot, where:
#         * The X-axis represents PC1.
#         * The Y-axis represents PC2.
#         * Each point corresponds to one metabolic model.

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(ggrepel)
})

input_file <- snakemake@input[["pca"]]
var_file <- snakemake@input[["var"]]
output_file <- snakemake@output[["png"]]

# Read PCA coordinates and variance information
pca_df <- fread(input_file)
var_df <- fread(var_file)

# Check that PC1 and PC2 are available
if (!all(c("PC1", "PC2") %in% names(pca_df))) {
  stop("PC1 and PC2 not found in PCA file")
}

# Check that the model identifier column is available
if (!("model" %in% names(pca_df))) {
  stop("Column 'model' not found in PCA file")
}

# Extract the percentage of variance explained by PC1 and PC2
pc1_var <- 100 * var_df$variance_explained[1]
pc2_var <- 100 * var_df$variance_explained[2]

# Generate PCA scatter plot
p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(
    aes(label = model),
    size = 3,
    max.overlaps = Inf
  ) +
  theme_classic() +
  labs(
    title = "PCA of median flux",
    x = sprintf("PC1 (%.1f%%)", pc1_var),
    y = sprintf("PC2 (%.1f%%)", pc2_var)
  )

# Create output directory if needed
dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

# Save PCA plot
ggsave(
  filename = output_file,
  plot = p,
  width = 6,
  height = 5,
  dpi = 300
)