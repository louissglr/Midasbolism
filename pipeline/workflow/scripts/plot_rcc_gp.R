#! /usr/bin/env Rscript

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

    to_del_it = which(rownames(cor_matrix) %in% to_del)

    if(length(to_del_it) == 0){
      return(cor_matrix)
    } else {
      return(
        cor_matrix[
          -which(rownames(cor_matrix) %in% to_del),
          -which(rownames(cor_matrix) %in% to_del)
        ]
      )
    }
  }

  if (strategy == "keep"){
    return(cor_matrix)
  }

  stop(paste0(strategy, " is not a valid strategy, please choose between 'zero', 'delete' or 'keep'."))
}

plot_and_save_histogram <- function(res_cor_sum, gp_status = NULL){

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
    df_tibble <- df_tibble %>% mutate(status = "no")
  }

  cumulative_cor_hist = ggplot(df_tibble, aes(x = cor_scaled, fill = status)) +
    geom_histogram(
      binwidth = 1,
      color = "darkgreen",
      alpha = 0.9,
      boundary = 0
    ) +
    scale_fill_manual(
      values = c("yes" = "#69b3a2", "no" = "#e07b6a"),
      name = "GapFilling"
    ) +
    xlim(c(0, 100)) +
    xlab("Scaled cumulative correlation") +
    ylab("Reaction count") +
    theme_bw()

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
gp_status_path         = snakemake@input[["gp_status"]]

#### main

# Open flux sampling matrix
fs_df = suppressMessages(read_tsv(file = flux_sample_matrix_path)) %>%
  select(-c(1))

# Compute correlation
res_cor = fastCor(fs_df, nSplit = 1, upperTri = FALSE, optBLAS = TRUE)

# Suppress NA
res_cor = handle_NA(res_cor, "delete")

# Calculate cumulative sum of correlation, row wise
res_cor_sum = apply(abs(res_cor), 1, sum)

# Load SRP status
gp_status = read_tsv(gp_status_path, show_col_types = FALSE) %>%
  rename(status = gp_status)

# Plot and save histogram profile
plot_and_save_histogram(res_cor_sum, gp_status)