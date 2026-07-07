# Description:
#   This script performs a Principal Component Analysis (PCA) on a median flux
#   matrix generated from multiple metabolic models obtained via the HUMESS
#   pipeline. The script also computes the proportion of variance explained by 
#   each principal component and the cumulative explained variance.
#
# Inputs:
#   - snakemake@input[["fs"]]:
#       TSV file containing:
#         * One row per model.
#         * One column per reaction retained after variance filtering.
#         * The first column contains the model ID.
#
# Outputs:
#   - snakemake@output[["pca"]]:
#       TSV file containing PCA coordinates where:
#         * Each row corresponds to a metabolic model.
#         * Columns correspond to principal component scores (PCs).
#         * The first column contains the model ID.
#
#   - snakemake@output[["var"]]:
#       TSV file containing PCA variance information where:
#         * Each row corresponds to one principal component.
#         * Contains the variance explained by each PC (col1 = variance_explained).
#         * Contains the cumulative variance explained across PCs (col2 = cumulative_variance).

# Logging
if (exists("snakemake")) {
  log_file <- snakemake@log[[1]]
  log_dir <- dirname(log_file)

  if (!dir.exists(log_dir))
    dir.create(log_dir, recursive = TRUE)

  mainlog <- file(log_file, open = "wt")

  sink(mainlog, type = "output")
  sink(mainlog, type = "message")

  on.exit(sink(type = "output"))
  on.exit(sink(type = "message"), add = TRUE)
  on.exit(close(mainlog), add = TRUE)
}

suppressPackageStartupMessages({
  library(data.table)
})

# Inputs
input_file <- snakemake@input[["fs"]]

# Outputs
pca_file <- snakemake@output[["pca"]]
var_file <- snakemake@output[["var"]]

cat("[CHECKPOINT] Input file:", input_file, "\n")

# Lecture
flux_table <- fread(input_file)

cat(
  "[CHECKPOINT] Input matrix dim:",
  nrow(flux_table), "x", ncol(flux_table),
  "\n"
)

# Récupération des modèles
models <- flux_table$model

# Matrice numérique
mat <- as.matrix(flux_table[, !"model"])

storage.mode(mat) <- "double"

cat(
  "[CHECKPOINT] PCA matrix dim:",
  nrow(mat), "x", ncol(mat),
  "\n"
)

# Nombre maximal de composantes possibles
rank_max <- min(nrow(mat) - 1, ncol(mat))

cat(
  "[CHECKPOINT] PCA rank:",
  rank_max,
  "\n"
)

# PCA
pca <- prcomp(
  x = mat,
  center = TRUE,
  scale. = FALSE,
)

# Coordonnées
pca_df <- data.table(
  model = models,
  pca$x
)

cat(
  "[CHECKPOINT] PCA scores dim:",
  nrow(pca_df), "x", ncol(pca_df),
  "\n"
)

# Écriture PCA
dir.create(dirname(pca_file),
           recursive = TRUE,
           showWarnings = FALSE)

fwrite(
  pca_df,
  pca_file,
  sep = "\t"
)

cat(
  "[CHECKPOINT] PCA written:",
  pca_file,
  "\n"
)

# Variance expliquée
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)

var_df <- data.table(
  PC = paste0("PC", seq_along(var_exp)),
  variance_explained = var_exp,
  cumulative_variance = cumsum(var_exp)
)

fwrite(
  var_df,
  var_file,
  sep = "\t"
)

cat(
  "[CHECKPOINT] Variance written:",
  var_file,
  "\n"
)

cat("[CHECKPOINT] Finished successfully\n")