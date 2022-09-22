#!/bin/bash

echo 'Which species do you want to analyse with Anvio? Type in the terminal the species with the format: "Genus_species"'
read

if [[ ! -d ~/PhD/Microbiome_pangenomic_analysis/data/$REPLY ]]
then
        echo -e "The species $REPLY is not availabe :(";
        echo -e "Maybe you have forgot to run gen_samples_txt.sh"
        exit
fi

if [[ ! -f ~/PhD/Microbiome_pangenomic_analysis/data/$REPLY/fasta-txt ]]
then
        echo -e "The species $REPLY is not availabe :(";
        echo -e "Maybe you have forgot to run gen_fasta_txt.sh"
        exit
fi

cd ~/PhD/Microbiome_pangenomic_analysis/data/$REPLY

eval "$(conda shell.bash hook)"
conda activate anvio-7.1

anvi-run-workflow -w metagenomics \
        -c ../../code/config.json