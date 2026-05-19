library(data.table)
library(uwot)

input_file  <- snakemake@input[["pca"]]
output_file <- snakemake@output[["umap"]]

pca_df <- fread(input_file)

n_pc <- min(100, ncol(pca_df))
mat <- as.matrix(pca_df[, 1:n_pc, with = FALSE])

set.seed(123)
umap_res <- umap(
    mat,
    n_neighbors = 30,
    min_dist = 0.1,
    metric = "euclidean"
)

umap_df <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2")

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

fwrite(umap_df, output_file, sep = "\t")
