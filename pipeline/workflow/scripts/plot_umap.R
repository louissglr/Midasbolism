library(data.table)
library(ggplot2)

input_file  <- snakemake@input[["umap"]]
output_file <- snakemake@output[["png"]]

umap_df <- fread(input_file)

p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(size = 1, alpha = 0.8) +
    theme_minimal()

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

ggsave(output_file, plot = p, width = 6, height = 5, dpi = 300)
