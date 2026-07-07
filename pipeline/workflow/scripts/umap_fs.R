# Description:
#   This script performs a Uniform Manifold Approximation and Projection (UMAP)
#   dimensionality reduction on PCA coordinates generated from flux sampling
#   solutions of a single metabolic model. The first principal components (=30) are
#   used as input for UMAP, and the resulting two-dimensional embedding is
#   saved as output.
#
# Inputs:
#   - snakemake@input[["pca"]]:
#       TSV file containing PCA coordinates, where:
#         * Rows correspond to flux samples (sampling solutions) from a single
#           metabolic model.
#         * Columns correspond to principal component scores.
#
#   - snakemake@wildcards[["models"]]:
#       ID of the single metabolic model associated with the PCA
#       coordinates.
#
# Outputs:
#   - snakemake@output[["umap"]]:
#       TSV file containing UMAP coordinates, where:
#         * Rows correspond to flux samples.
#         * Columns contain the two-dimensional UMAP embedding:
#             - "UMAP1"
#             - "UMAP2"

suppressPackageStartupMessages({
  library(data.table)
library(uwot)
})

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