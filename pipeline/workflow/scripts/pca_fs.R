library(readr)
library(dplyr)
library(data.table)

# Input  :
# rows = sampling solutions
# columns = flux reactions
input_file  <- snakemake@input[["fs"]]
sampled_file <- snakemake@input[["sampled"]]

# Output (pca$x) :
# rows = flux reactions
# columns = PCs
output_file <- snakemake@output[["pca"]]

flux_table <- suppressMessages(
    read_tsv(file = input_file)
)

# Load sampled IDs
sampled_ids <- fread(sampled_file)[[1]]

# Identify ID column 
id_col <- names(flux_table)[1]

# Subsampling
flux_table <- flux_table %>%
    filter(.data[[id_col]] %in% sampled_ids)

# Keep IDs for later
ids <- flux_table[[id_col]]

flux_table <- flux_table |>
    dplyr::select(-c(1))

# transpose because:
# rows = features 
# columns = samples

#mat <- t(flux_table)

mat <- as.matrix(flux_table)

pca <- prcomp(
    mat,
    center = TRUE,
    scale. = TRUE
)

pca_df <- as.data.frame(pca$x) #loadings

# add IDs back (traceability like in PCoA)
pca_df[[id_col]] <- ids
pca_df <- pca_df %>% dplyr::select(all_of(id_col), everything())

dir.create(
    dirname(output_file),
    recursive = TRUE,
    showWarnings = FALSE
)

fwrite(
    pca_df,
    file = output_file,
    sep = "\t"
)