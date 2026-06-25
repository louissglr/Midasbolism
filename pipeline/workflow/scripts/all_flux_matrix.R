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
  library(data.table)
  library(irlba)
  library(ggplot2)
})

cat("[CHECKPOINT] Starting PCA on flux sampling\n")

all_samples <- list()

for (fs_file in snakemake@input) {

  model <- sub(
    "^fs_(.*)\\.tsv$",
    "\\1",
    basename(fs_file)
  )

  cat("[CHECKPOINT] Processing model:", model, "\n")

  fs <- fread(fs_file)

  cat(
    "[CHECKPOINT] Raw matrix dim:",
    nrow(fs), "x", ncol(fs),
    "\n"
  )

  # Remove solution ID column
  fs <- fs[, -1]

  # Add group label
  fs[, model := model]

  setcolorder(
    fs,
    c("model", setdiff(names(fs), "model"))
  )

  all_samples[[model]] <- fs
}

cat("[CHECKPOINT] Merging sampling matrices\n")

res <- rbindlist(
  all_samples,
  fill = TRUE
)

rm(all_samples)
gc()

cat(
  "[CHECKPOINT] Combined matrix dim:",
  nrow(res), "x", ncol(res),
  "\n"
)

n_na <- sum(is.na(res))

cat(
  "[CHECKPOINT] Missing values before replacement:",
  n_na,
  "\n"
)

res[is.na(res)] <- 0

model_labels <- res$model

flux_mat <- as.matrix(
  res[, !"model"]
)

rm(res)
gc()

storage.mode(flux_mat) <- "double"

cat(
  "[CHECKPOINT] Numeric matrix dim:",
  nrow(flux_mat), "x", ncol(flux_mat),
  "\n"
)

cat("[CHECKPOINT] Removing constant reactions\n")

var_reac <- apply(
  flux_mat,
  2,
  var
)

cat(
  "[CHECKPOINT] Reactions kept:",
  sum(var_reac > 0),
  "/",
  length(var_reac),
  "\n"
)

flux_mat <- flux_mat[
  ,
  var_reac > 0,
  drop = FALSE
]

rm(var_reac)
gc()

cat(
  "[CHECKPOINT] Filtered matrix dim:",
  nrow(flux_mat), "x", ncol(flux_mat),
  "\n"
)

cat("[CHECKPOINT] Running PCA\n")

pca <- prcomp_irlba(
  flux_mat,
  n = 2,
  center = TRUE,
  scale. = TRUE
)

# pca <- prcomp(
#   flux_mat,
#   center = TRUE,
#   scale. = TRUE
# )

cat("[CHECKPOINT] PCA completed\n")

scores_plot <- data.table(
  model = model_labels,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

rm(
  flux_mat,
  model_labels,
  pca
)
gc()

cat("[CHECKPOINT] Sampling points for plotting\n")

scores_plot <- scores_plot[
  ,
  .SD[
    sample(
      .N,
      min(.N, 3000)
    )
  ],
  by = model
]

cat(
  "[CHECKPOINT] Points plotted:",
  nrow(scores_plot),
  "\n"
)

p <- ggplot(
  scores_plot,
  aes(
    x = PC1,
    y = PC2,
    color = model
  )
) +
  geom_point(
    alpha = 0.5,
    size = 0.8
  ) +
  theme_bw() +
  labs(
    title = "PCA of flux sampling solution spaces",
    x = "PC1",
    y = "PC2",
    color = "Group"
  )

ggsave(
  filename = snakemake@output[["png"]],
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)

cat("[CHECKPOINT] Figure saved\n")
cat("[CHECKPOINT] Finished successfully\n")