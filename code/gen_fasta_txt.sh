#!/bin/bash

cd ~/PhD

echo 'Which bacteria do you want to analyse? Type it in the terminal with the format: "Genus_species"'
read

WORKDIR=Microbiome_pangenomic_analysis/data/$REPLY

echo -e "name\tpath\n503_G175_R4_T1_MRS_4\t$WORKDIR/503_G175_R4_T1_MRS_4.fasta" > $WORKDIR/fasta-txt