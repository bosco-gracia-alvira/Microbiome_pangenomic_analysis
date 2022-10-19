#!/bin/bash

cd ~/PhD/data/Microbiome_pangenomic_analysis/data

echo 'Which Anvio profile do you want to display? Type in the terminal the species with the format: "Genus_species"'
read

CONTIG=$(cat $REPLY/fasta-txt | awk '!/name/ {print $1}')

eval "$(conda shell.bash hook)"
conda activate anvio-7.1

anvi-interactive -c $REPLY/Anvio/03_CONTIGS/$CONTIG-contigs.db \
                 -p $REPLY/Anvio/06_MERGED/$CONTIG/PROFILE.db