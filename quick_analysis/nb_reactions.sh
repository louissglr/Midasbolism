#!/usr/bin/env bash

# BEGIN MANDATORY OPTIONS
#SBATCH --time=0-00:05:00
#SBATCH --qos=debug
# END MANDATORY OPTIONS

# BEGIN INFORMATIONAL OPTIONS
#SBATCH --job-name=nb_reactions
#SBATCH --comment="nb_reactions"
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
# END INFORMATIONAL OPTIONS

# BEGIN RESOURCES OPTIONS
#SBATCH --partition=standard
#SBATCH --ntasks=1
# END RESOURCES OPTIONS

dataDir="/LAB-DATA/GLiCID/projects/CRCI2NA_DATA/ICAGEN/MIDAS/Analyses/RNA-seq/METABOLISM/HUMESS/models/edgeR/MIDAS_exp.cluster/output"

models=(CD1 CD2 HP1 HP2 HP3 LB MF MN MS)

echo -e "model\tsrp\tcarveme\tfs_reactions\tgapfilled\tintersection_carveme_fs"

for model in "${models[@]}"; do

    carveme="${dataDir}/models/${model}/stats/carveme.reactions.list"
    srp="${dataDir}/models/${model}/stats/srp.reactions.list"
    fs="${dataDir}/models/${model}/fs/fs_${model}.tsv"

    n_carveme=$(sort -u "$carveme" | wc -l)
    n_srp=$(sort -u "$srp" | wc -l)

    # Réactions ajoutées par gap-filling
    gapfilled=$(comm -23 <(sort -u "$carveme") <(sort -u "$srp") | wc -l) #nb réactions présentes dans carveMe mais absente de SRP

    # Réactions de FS (noms des colonnes sauf la première)
    fs_tmp=$(mktemp)
    head -n1 "$fs" | tr '\t' '\n' | tail -n +2 | sort -u > "$fs_tmp"

    fs_reactions=$(wc -l < "$fs_tmp")

    # Intersection CarveMe ∩ FS
    intersection_carveme_fs=$(comm -12 \
        <(sort -u "$carveme") \
        "$fs_tmp" | wc -l)

    echo -e "${model}\t${n_srp}\t${n_carveme}\t${fs_reactions}\t${gapfilled}\t${intersection_carveme_fs}"

    rm -f "$fs_tmp"

done