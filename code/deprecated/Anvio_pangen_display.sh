#!/bin/bash

echo 'Which Anvio profile do you want to display? Type in the terminal the species with the format: "Genus_species"'
read

cd ~/PhD/Microbiome_pangenomic_analysis/data/$REPLY/Anvio_pangen

CONTIG=$(cat fasta-txt | awk '!/name/ {print $1}')

eval "$(conda shell.bash hook)"
conda activate anvio-7.1

anvi-display-pan -g 03_PAN/MYPAN-GENOMES.db \
                 -p 03_PAN/MYPAN-PAN.db \
                 --title $REPLY