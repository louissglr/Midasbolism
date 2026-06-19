library(data.table)

all_medians <- list()

for (fs_file in snakemake@input) {

    # Récupération du nom du modèle depuis le nom du fichier
    model <- sub("^fs_(.*)\\.tsv$", "\\1", basename(fs_file))

    # Lecture du fichier
    fs <- fread(fs_file)

    # Suppression de la première colonne (IDs des solutions)
    fs <- fs[, -1]

    # Médiane de chaque réaction
    med <- fs[, lapply(.SD, median, na.rm = TRUE)]

    # Ajout du modèle
    med[, model := model]

    # Mise en première position de la colonne model
    setcolorder(med, c("model", setdiff(names(med), "model")))

    all_medians[[model]] <- med
}

# Fusion de tous les modèles
res <- rbindlist(all_medians, fill = TRUE)

# Écriture du résultat
fwrite(
    res,
    snakemake@output[["tsv"]],
    sep = "\t"
)