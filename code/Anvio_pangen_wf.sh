#!/bin/bash

cd ~/PhD

echo 'Which species do you want to analyse with Anvio? Type in the terminal the species with the format: "Genus_species"'
echo
cat Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv | cut -f2 | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
read

if [[ ! -d ~/PhD/Microbiome_pangenomic_analysis/data/$REPLY/Anvio_pangen ]]
then
        echo -e "The species $REPLY is not availabe :(";
        echo -e "Maybe you have forgotten to run gen_pangen_fasta_txt.sh"
        exit
fi

cd ~/PhD/Microbiome_pangenomic_analysis/data/$REPLY/Anvio_pangen

eval "$(conda shell.bash hook)"
conda activate anvio-7.1

anvi-run-workflow -w pangenomics \
                  -c ../../../code/config-pangen.json