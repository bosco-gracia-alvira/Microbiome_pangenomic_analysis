#!/bin/bash

cd ~/PhD/data/Microbiome_pangenomic_analysis/data

echo 'Which Anvio profile do you want to annotate? Type in the terminal the species with the format: "Genus_species"'
read

CONTIG=cat $REPLY/fasta-txt | awk '!/name/ {print $1}'

eval "$(conda shell.bash hook)"
conda activate anvio-7.1

anvi-run-pfams -c $REPLY/03_CONTIGS/$CONTIG-contigs.db -T 4
anvi-run-ncbi-cogs -c $REPLY/03_CONTIGS/$CONTIG-contigs.db -T 4
anvi-run-kegg-kofams -c $REPLY/03_CONTIGS/$CONTIG-contigs.db -T 4