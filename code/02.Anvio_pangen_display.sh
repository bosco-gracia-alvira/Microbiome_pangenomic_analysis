#!/bin/bash
# This script runs Anvio pangenomic analysis
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data"

### COMMANDS
IFS="
"

# Ask the user which species to analyse
echo 'Which species do you want to analyse with Anvio? Type in the terminal the species with the format: "Genus_species"'
echo
cut -f2 "$WORKDIR"/taxonomy.tsv | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
read

# If the fasta-txt is not available we cannot do the analysis!
if [[ ! -f "$WORKDIR"/$REPLY/Anvio/03_PAN/MYPAN-PAN.db ]]
then    
        echo
        echo
        echo -e "The species $REPLY is not availabe :("
        echo -e "Maybe you have forgotten to run 02.Anvio_pangen_wf.sh"
        echo
        exit
fi

eval "$(conda shell.bash hook)"
conda activate anvio-7.1

anvi-display-pan \
        -g "$WORKDIR"/"$REPLY"/Anvio/03_PAN/MYPAN-GENOMES.db \
        -p "$WORKDIR"/"$REPLY"/Anvio/03_PAN/MYPAN-PAN.db
