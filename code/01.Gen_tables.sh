#!/bin/bash
# This script creates the file fasta-txt, that is the input for Anvio pangenomic analysis.
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data"
ASSEMBLY="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly"
GENOMES="$WORKDIR"/Isolates

### COMMANDS
IFS="
"

# Create the taxonomy, checkm and metadata files
 cat "$ASSEMBLY"/Pool_503/07.GTDB-Tk/gtdbtk.bac120.summary.tsv > "$WORKDIR"/taxonomy.tsv
 tail -n +2 "$ASSEMBLY"/Pool_591/07.GTDB-Tk/gtdbtk.bac120.summary.tsv >> "$WORKDIR"/taxonomy.tsv
 tail -n +2 "$ASSEMBLY"/Pool_643/07.GTDB-Tk/gtdbtk.bac120.summary.tsv >> "$WORKDIR"/taxonomy.tsv

 cat "$ASSEMBLY"/Pool_503/04.CheckM2/quality_report.tsv > "$WORKDIR"/checkm.tsv
 tail -n +2 "$ASSEMBLY"/Pool_591/04.CheckM2/quality_report.tsv >> "$WORKDIR"/checkm.tsv
 tail -n +2 "$ASSEMBLY"/Pool_643/04.CheckM2/quality_report.tsv >> "$WORKDIR"/checkm.tsv

 cat "$ASSEMBLY"/Pool_503/metadata.tsv > "$WORKDIR"/metadata.tsv
 tail -n +2 "$ASSEMBLY"/Pool_591/metadata.tsv >> "$WORKDIR"/metadata.tsv
 tail -n +2 "$ASSEMBLY"/Pool_643/metadata.tsv >> "$WORKDIR"/metadata.tsv

# The R script has to be run using the native R. It does not work with Anvio's R
conda deactivate

# This script merges the taxonomy and metadata files
Rscript "$WORKDIR"/../code/01.Reformat_metadata.R

# First argument of the script is the species to analyse
REPLY=$1

if [[ -z "$REPLY" ]]
then
    echo "You need to provide the species you want to analyse as first argument"
    echo
    echo "These are the species that are available for the pangenomic analysis"
    echo "Copy the name of the species you want to analyse and paste it as first argument of the script"
    echo
    cut -f2 "$WORKDIR"/taxonomy.tsv | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
    echo
    exit
elif [[ "$REPLY" == "h" ]]
then
    echo "These are the species that are available for the pangenomic analysis"
    echo "Copy the name of the species you want to analyse and paste it as first argument of the script"
    echo
    cut -f2 "$WORKDIR"/taxonomy.tsv | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
    echo
    exit
fi

# Set the variable SPECIES from the input
SPECIES=$(echo "$REPLY" | sed 's/_/ /')

# Create the Anvio directory
if [[ ! -f "$WORKDIR"/$REPLY/Anvio ]]
then
    mkdir -p "$WORKDIR"/$REPLY/Anvio
fi

# Create the genomes directory and copy all the available genomes to it
if [[ ! -f "$GENOMES" ]]
then
    mkdir -p "$GENOMES"
fi
cp "$ASSEMBLY"/Pool_??3/07.GTDB-Tk/Genomes/*.fasta "$GENOMES"

# Pick the samples that belong to the species of interest
SAMPLES=$(awk -v s="$SPECIES" -F "\t" '$2 ~ s {print $1}' "$WORKDIR"/taxonomy.tsv | grep -v "user" | sed 's/-/_/g')

# Create a new fasta-txt file in the Anvio directory
head -n 1 "$WORKDIR"/fasta-txt > "$WORKDIR"/$REPLY/Anvio/fasta-txt

# Populate it with the samples of interest
for i in $SAMPLES
do
        grep -w "$i" "$WORKDIR"/fasta-txt >> "$WORKDIR"/$REPLY/Anvio/fasta-txt
done

# Create a new samples-txt file in the Anvio directory
head -n 1 "$WORKDIR"/samples-txt > "$WORKDIR"/$REPLY/Anvio/samples-txt

# Populate it with the samples of interest
for i in $SAMPLES
do
        grep -w "$i" "$WORKDIR"/samples-txt >> "$WORKDIR"/$REPLY/Anvio/samples-txt
done

# Create a new misc file in the Anvio directory
head -n 1 "$WORKDIR"/Anvio_misc.tsv > "$WORKDIR"/$REPLY/Anvio/Anvio_misc.tsv

# Populate it with the samples of interest
for i in $SAMPLES
do
        grep -w "$i" "$WORKDIR"/Anvio_misc.tsv >> "$WORKDIR"/$REPLY/Anvio/Anvio_misc.tsv
done
