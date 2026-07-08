#! /usr/bin/env Rscript

# Description:
#   This script (adapted to Humess - Louis Paré) computes the pairwise
#   correlation between reactions from a flux sampling matrix. It handles
#   missing values resulting from non-variable reactions according to the
#   selected strategy, calculates the cumulative absolute correlation score
#   for each reaction, integrates gap-filling status information, rescales
#   the scores to a 0–100 range, and generates a stacked histogram showing
#   the distribution of cumulative correlations for gap-filled and
#   non-gap-filled reactions. The histogram is saved as a PNG file.

#### libraries
suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
library(scales, quietly = TRUE)
library(HiClimR, quietly = TRUE)

#### functions

is_all_na <- function(list_to_check){
  if(length(which(!unlist(lapply(list_to_check, is.na)))) > 0){
    return(FALSE)
  } else {
    return(TRUE)
  }
}

handle_NA <- function(cor_matrix, strategy = "zero"){

  if (strategy == "zero"){
    cor_matrix[is.na(cor_matrix)] <- 0
    return(cor_matrix)
  }

  if (strategy == "delete"){

    to_del = c()

    for (i in 1:ncol(cor_matrix)) {
      if(is_all_na(cor_matrix[i,])){
        to_del = c(to_del, colnames(cor_matrix)[i])
      }
    }

    if(length(to_del) == 0){
      return(cor_matrix)
    } else {
      return(
        cor_matrix[
          -which(rownames(cor_matrix) %in% to_del),
          -which(colnames(cor_matrix) %in% to_del)
        ]
      )
    }
  }

  if (strategy == "keep"){
    return(cor_matrix)
  }

  stop(
    paste0(
      strategy,
      " is not a valid strategy, please choose between 'zero', 'delete' or 'keep'."
    )
  )
}

plot_and_save_histogram <- function(
  res_cor_sum,
  gp_status = NULL,
  model_name = NULL
){

  df_tibble = tibble(
    reaction = names(res_cor_sum),
    cor_sum = res_cor_sum,
    cor_scaled = scales::rescale(
      x = res_cor_sum,
      from = c(min(res_cor_sum), max(res_cor_sum)),
      to = c(0, 100)
    )
  )

  if (!is.null(gp_status)) {
    df_tibble <- df_tibble %>%
      left_join(gp_status, by = "reaction") %>%
      mutate(status = replace_na(status, "yes"))
  } else {
    df_tibble <- df_tibble %>%
      mutate(status = "no")
  }

  # Ordre d'empilement : no en bas, yes au-dessus
  df_tibble <- df_tibble %>%
    mutate(
      status = factor(
        status,
        levels = c("no", "yes")
      )
    )

  cumulative_cor_hist = ggplot(
    df_tibble,
    aes(x = cor_scaled, fill = status)
  ) +
    geom_histogram(
      binwidth = 1,
      color = "darkgreen",
      alpha = 0.9,
      boundary = 0,
      position = position_stack(reverse = TRUE)
    ) +
    scale_fill_manual(
      values = c(
        "no" = "#e07b6a",
        "yes" = "#69b3a2"
      ),
      breaks = c("no", "yes"),
      name = "GapFilling"
    ) +
    labs(
      title = paste(
        "RCC -",
        model_name
      ),
      x = "Scaled cumulative correlation",
      y = "Reaction count"
    ) +
    xlim(c(0, 100)) +
    theme_bw() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      )
    )

  ggsave(
    filename = plot_output_path,
    plot = cumulative_cor_hist,
    device = "png",
    width = 1800,
    height = 1200,
    units = "px"
  )
}

#### parameters

flux_sample_matrix_path = snakemake@input[["fs"]]
plot_output_path        = snakemake@output[["png"]]
gp_status_path          = snakemake@input[["gp_status"]]

# Récupération du nom du modèle depuis le wildcard Snakemake
model_name = snakemake@wildcards[["models"]]

#### main

# Open flux sampling matrix
fs_df = suppressMessages(
  read_tsv(file = flux_sample_matrix_path)
) %>%
  select(-c(1))

# Compute correlation matrix
res_cor = fastCor(
  fs_df,
  nSplit = 1,
  upperTri = FALSE,
  optBLAS = TRUE
)

# Remove NA-only rows/columns
res_cor = handle_NA(res_cor, "delete")

# Calculate cumulative correlation score
res_cor_sum = apply(abs(res_cor), 1, sum)

# Load gapfilling status
gp_status = read_tsv(
  gp_status_path,
  show_col_types = FALSE
) %>%
  rename(status = gp_status)

# Plot and save histogram
plot_and_save_histogram(
  res_cor_sum = res_cor_sum,
  gp_status = gp_status,
  model_name = model_name
)