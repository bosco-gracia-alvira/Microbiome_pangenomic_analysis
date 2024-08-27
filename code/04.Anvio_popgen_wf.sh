#!/bin/bash
# This script runs Anvio population genomics analysis for the desired species
# Bosco Gracia Alvira, 2024

### VARIABLES

WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data"

# First argument of the script is the species to analyse
REPLY=$1

if [[ -z "$REPLY" ]]
then
    echo "You need to provide the species you want to analyse as first argument"
    echo
    echo "These are the species that are available for the pangenomic analysis"
    echo "Copy the name of the species you want to analyse and paste it as first argument of the script"
    echo
    cut -f2 "$WORKDIR"/taxonomy.tsv | awk -F 's__' '{print $2}' | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
    echo
    exit
elif [[ "$REPLY" == "h" ]]
then
    echo "These are the species that are available for the pangenomic analysis"
    echo "Copy the name of the species you want to analyse and paste it as first argument of the script"
    echo
    cut -f2 "$WORKDIR"/taxonomy.tsv | awk -F 's__' '{print $2}' | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
    echo
    exit
fi

# Set the paths
SPECIES=$(echo "$REPLY" | sed 's/_/ /')
SAMPLES=$(awk -v s="$SPECIES" -F "\t" '$2 ~ s {print $1}' "$WORKDIR"/taxonomy.tsv | grep -v "user" | sed 's/-/_/g')

### COMMANDS
IFS="
"

if [[ ! -f "$WORKDIR"/samples-txt ]]
then
        echo -e "The species $REPLY is not availabe :("
        echo -e "I don't find the samples-txt file"
        echo -e "You need to run 01.Gen_tables.sh"
        exit
fi

if [[ ! -f "$WORKDIR"/"$REPLY"/Anvio_popgen ]]
then
    mkdir -p "$WORKDIR/$REPLY/Anvio_popgen"
fi

eval "$(conda shell.bash hook)"
conda activate anvio-7.1

# Copy the assembly into the Anvio directory and reformat it
anvi-script-reformat-fasta \
        -l 250 \
        -o "$WORKDIR/$REPLY/Anvio_popgen/${REPLY}_pangenome.fasta" \
        --simplify-names \
        "$WORKDIR/$REPLY/SuperPang/assembly.fasta"

# Create the reference-txt file in the Anvio directory
echo -e "name\tpath" > "$WORKDIR"/$REPLY/Anvio_popgen/reference-txt
echo -e "${REPLY}\t./${REPLY}_pangenome.fasta" >> "$WORKDIR"/$REPLY/Anvio_popgen/reference-txt

# Create a new samples-txt file in the Anvio directory
head -n 1 "$WORKDIR"/samples-txt > "$WORKDIR/$REPLY/Anvio_popgen/samples-txt"

# Populate it with the samples of interest
for i in $SAMPLES
do
        grep -w "$i" "$WORKDIR/samples-txt" >> "$WORKDIR/$REPLY/Anvio_popgen/samples-txt"
done

cd "$WORKDIR/$REPLY/Anvio_popgen" || exit

anvi-run-workflow -w metagenomics \
        -c ../../../code/04.config-popgen.json

# This is highly experimental
# anvi-script-add-default-collection -p "06_MERGED/${REPLY}/PROFILE.db"

# anvi-gen-variability-profile -c "02_CONTIGS/${REPLY}-contigs.db" \
#                              -p "06_MERGED/${REPLY}/PROFILE.db" \
#                              -C DEFAULT \
#                              -b EVERYTHING \
#                              --engine CDN \
#                              --include-site-pnps \
#                              --kiefl-mode \
#                              -o "${REPLY}-SCVs.txt"

# anvi-gen-variability-profile -c "02_CONTIGS/${REPLY}-contigs.db" \
#                              -p "06_MERGED/${REPLY}/PROFILE.db" \
#                              -C DEFAULT \
#                              -b EVERYTHING \
#                              --engine AA \
#                              --kiefl-mode \
#                              -o $REPLY-SAAVs.txt

# awk '{ if ((NR == 1) || ($73!=$74)) { print } }' "${REPLY}-SCVs.txt" > consensus-SCVs.tsv
# awk '{ if ((NR == 1) || ($30!=$31)) { print } }' "${REPLY}-SAAVs.txt" > consensus-SAAVs.tsv
