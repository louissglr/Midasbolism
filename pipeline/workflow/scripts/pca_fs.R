# Logging
if (exists("snakemake")) {
  log_file <- snakemake@log[[1]]
  log_dir <- dirname(log_file)

  if (!dir.exists(log_dir))
    dir.create(log_dir, recursive = TRUE)

  mainlog <- file(log_file, open = "wt")

  sink(mainlog, append = FALSE, type = "output")
  sink(mainlog, append = FALSE, type = "message")

  on.exit(sink(type = "output"))
  on.exit(sink(type = "message"), add = TRUE)
  on.exit(close(mainlog), add = TRUE)
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(data.table)
})

# Inputs
input_file <- snakemake@input[["fs"]]

# Wildcard Snakemake
model_name <- snakemake@wildcards[["models"]]

# Params
n <- snakemake@params[["n"]]
seed <- snakemake@params[["seed"]]
sampling <- snakemake@params[["sampling"]]

# Output
output_file <- snakemake@output[["pca"]]

set.seed(seed)

cat("[CHECKPOINT] model:", model_name, "\n")
cat("[CHECKPOINT] input_file:", input_file, "\n")
cat("[CHECKPOINT] n:", n, "seed:", seed, "\n")
cat("[CHECKPOINT] sampling:", sampling, "\n")

# Read full matrix
flux_table <- suppressMessages(
  read_tsv(file = input_file)
)

cat("[CHECKPOINT] flux_table dim:", dim(flux_table), "\n")

# --- CONDITIONAL SAMPLING ---
n_total <- nrow(flux_table)

if (isTRUE(sampling)) {

  n_keep <- min(n, n_total)

  sample_idx <- sample(
    seq_len(n_total),
    n_keep,
    replace = FALSE
  )

  flux_table_sub <- flux_table[sample_idx, ]

  cat("[CHECKPOINT] sampling ENABLED\n")

} else {

  flux_table_sub <- flux_table
  n_keep <- n_total

  cat("[CHECKPOINT] sampling DISABLED\n")
}

cat("[CHECKPOINT] flux_table_sub dim:", dim(flux_table_sub), "\n")
cat("[CHECKPOINT] Total rows:", n_total, "\n")
cat("[CHECKPOINT] Used rows:", n_keep, "\n")

# Numeric matrix for PCA
mat <- flux_table_sub %>%
  select(-1) %>%   # première colonne = identifiant
  as.matrix()

storage.mode(mat) <- "double"

cat("[CHECKPOINT] matrix dim:", dim(mat), "\n")

# PCA
pca <- stats::prcomp(
  x = mat,
  rank. = 200, #maximal number of principal components to be used
  center = TRUE, # a logical value indicating whether the variables should be shifted to be zero centered
  scale. = FALSE #a logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place.
)

pca_df <- as.data.frame(pca$x)

cat("[CHECKPOINT] PCA output dim:", dim(pca_df), "\n")
cat("[CHECKPOINT] prcomp requested rank: 200\n")

# Write PCA output
dir.create(
  dirname(output_file),
  recursive = TRUE,
  showWarnings = FALSE
)

data.table::fwrite(
  pca_df,
  file = output_file,
  sep = "\t"
)

cat("[CHECKPOINT] PCA written to:", output_file, "\n")

# Variance expliquée
var_exp <- (pca$sdev^2) / sum(pca$sdev^2)

var_df <- data.frame(
  PC = paste0("PC", seq_along(var_exp)),
  variance_explained = var_exp
)

var_file <- sub(
  "_pca.tsv",
  "_pca_var.tsv",
  output_file
)

data.table::fwrite(
  var_df,
  file = var_file,
  sep = "\t"
)

cat("[CHECKPOINT] variance explained written to:", var_file, "\n")

# Exemple d'utilisation du nom du modèle
cat("[CHECKPOINT] PCA completed for model:", model_name, "\n")