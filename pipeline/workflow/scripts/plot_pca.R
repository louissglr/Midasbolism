library(readr)
library(ggplot2)
library(data.table)

input_file  <- snakemake@input[["pca"]]
output_file <- snakemake@output[["png"]]

pca_df <- fread(input_file)

if (!all(c("PC1", "PC2") %in% colnames(pca_df))) {
    stop("PC1 et PC2 introuvables dans le fichier PCA")
}


p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point(alpha = 0.7, size = 1) +
    theme_minimal() +
    labs(
        title = "PCA (PC1 vs PC2)",
        x = "PC1",
        y = "PC2"
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
