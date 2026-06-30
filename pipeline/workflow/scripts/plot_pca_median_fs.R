library(ggplot2)
library(data.table)
library(ggrepel)

input_file <- snakemake@input[["pca"]]
var_file <- snakemake@input[["var"]]
output_file <- snakemake@output[["png"]]

pca_df <- fread(input_file)
var_df <- fread(var_file)

if (!all(c("PC1", "PC2") %in% names(pca_df))) {
  stop("PC1 et PC2 introuvables dans le fichier PCA")
}

if (!("model" %in% names(pca_df))) {
  stop("Colonne 'model' introuvable dans le fichier PCA")
}

pc1_var <- 100 * var_df$variance_explained[1]
pc2_var <- 100 * var_df$variance_explained[2]

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

dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  filename = output_file,
  plot = p,
  width = 6,
  height = 5,
  dpi = 300
)