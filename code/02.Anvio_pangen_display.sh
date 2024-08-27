#!/bin/bash
# This script runs Anvio pangenomic analysis
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data"

### COMMANDS
IFS="
"

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

# If the fasta-txt is not available we cannot do the analysis!
if [[ ! -f "$WORKDIR"/$REPLY/Anvio_pangen/03_PAN/MYPAN-PAN.db ]]
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
        -g "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/MYPAN-GENOMES.db \
        -p "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/MYPAN-PAN.db
