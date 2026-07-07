# Description:
#   This script reads multiple flux sampling matrices (one per metabolic
#   model obtained via the HUMESS pipeline), computes the median flux for
#   each reaction within each model, combines all models into a single
#   matrix, replaces missing reactions with zero flux, removes reactions
#   with zero variance across models, and writes the resulting matrix to
#   a TSV file.
#
# Inputs:
#   - snakemake@input:
#       Multiple TSV files named "fs_<model>.tsv".
#       Each file contains flux sampling results where:
#         * Rows correspond to sampled solutions (~150k).
#         * Columns correspond to reactions.
#         * The first column contains solution ID
#
# Outputs:
#   - snakemake@output[["tsv"]]:
#       TSV file containing:
#         * One row per model.
#         * One column per reaction retained after variance filtering.
#         * The first column contains the model ID.


# Logging
if (exists("snakemake")) {
  log_file <- snakemake@log[[1]]
  log_dir <- dirname(log_file)
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

  mainlog <- file(log_file, open = "wt")
  sink(mainlog, append = FALSE, type = "output")
  sink(mainlog, append = FALSE, type = "message")

  on.exit(sink(type = "output"))
  on.exit(sink(type = "message"), add = TRUE)
  on.exit(close(mainlog), add = TRUE)
}

suppressPackageStartupMessages({
  library(data.table)
})

cat("[CHECKPOINT] Starting median flux aggregation\n")

all_medians <- list()

for (fs_file in snakemake@input) {

  model <- sub("^fs_(.*)\\.tsv$", "\\1", basename(fs_file))

  cat("[CHECKPOINT] Processing model:", model, "\n")

  fs <- fread(fs_file)

  cat(
    "[CHECKPOINT] Raw sampling matrix dim:",
    nrow(fs), "x", ncol(fs),
    "\n"
  )

  # Remove solution ID column
  fs <- fs[, -1]

  cat(
    "[CHECKPOINT] Flux matrix dim:",
    nrow(fs), "x", ncol(fs),
    "\n"
  )

  # Compute the median flux for each reaction
  med <- fs[, lapply(.SD, median, na.rm = TRUE)]

  med[, model := model]
  setcolorder(med, c("model", setdiff(names(med), "model")))

  all_medians[[model]] <- med
}

# Merge all models
res <- rbindlist(all_medians, fill = TRUE)

cat(
  "[CHECKPOINT] Combined matrix dim (models x reactions):",
  nrow(res), "x", ncol(res),
  "\n"
)


# Count missing values before replacement
n_na <- sum(is.na(res))

cat(
  "[CHECKPOINT] Missing values before replacement:",
  n_na,
  "\n"
)

# Missing reactions are assigned a flux of zero
res[is.na(res)] <- 0

cat(
  "[CHECKPOINT] Missing values replaced by 0:",
  n_na,
  "\n"
)

# Extract the numeric matrix without the model column
flux_mat <- as.matrix(res[, !"model"])

storage.mode(flux_mat) <- "double"

cat(
  "[CHECKPOINT] Numeric matrix dim:",
  nrow(flux_mat), "x", ncol(flux_mat),
  "\n"
)

# Compute reaction variance across models
var_reac <- apply(flux_mat, 2, var)

n_total_reactions <- length(var_reac)
n_zero_var <- sum(var_reac == 0)
n_kept <- sum(var_reac > 0)

cat(
  "[CHECKPOINT] Total reactions:",
  n_total_reactions,
  "\n"
)

cat(
  "[CHECKPOINT] Reactions with variance = 0:",
  n_zero_var,
  "\n"
)

cat(
  "[CHECKPOINT] Reactions kept:",
  n_kept,
  "\n"
)

# List removed reactions
removed_reactions <- names(var_reac[var_reac == 0])

cat(
  "[CHECKPOINT] First 20 removed reactions:",
  paste(head(removed_reactions, 20), collapse = ", "),
  "\n"
)

# Filter reactions with non-zero variance
flux_mat_filt <- flux_mat[, var_reac > 0, drop = FALSE]

cat(
  "[CHECKPOINT] Filtered matrix dim:",
  nrow(flux_mat_filt), "x", ncol(flux_mat_filt),
  "\n"
)

# Rebuild the data.table
res_filt <- data.table(
  model = res$model,
  flux_mat_filt
)

cat(
  "[CHECKPOINT] Output matrix dim:",
  nrow(res_filt), "x", ncol(res_filt),
  "\n"
)

out_file <- snakemake@output[["tsv"]]

dir.create(dirname(out_file),
           recursive = TRUE,
           showWarnings = FALSE)

fwrite(res_filt, out_file, sep = "\t")

cat(
  "[CHECKPOINT] Output written to:",
  out_file,
  "\n"
)

cat("[CHECKPOINT] Finished successfully\n")