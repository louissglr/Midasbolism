library(data.table)
library(uwot)

input_file  <- snakemake@input[["pca"]]
output_file <- snakemake@output[["umap"]]

model_name <- snakemake@wildcards[["models"]]

pca_df <- fread(input_file)

n_pc <- min(30, ncol(pca_df))
mat <- as.matrix(pca_df[, 1:n_pc, with = FALSE])

set.seed(123)

knn_umap = 30
umap_res <- uwot::umap(
        X = mat,
        n_neighbors = knn_umap # 30 default for Seurat::RunUMAP() (used in muscadet::cluster_seurat())
    )

umap_df <- as.data.frame(umap_res)
colnames(umap_df) <- c("UMAP1", "UMAP2")


dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

fwrite(umap_df, output_file, sep = "\t")