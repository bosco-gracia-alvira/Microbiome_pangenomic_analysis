#!/bin/bash

cd ~/PhD

echo 'To which species do you want to generate a misc-table? Type in the terminal the species with the format: "Genus_species"'
echo
cat Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv | cut -f2 | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
read

if [[ ! -d ~/PhD/Microbiome_pangenomic_analysis/data/$REPLY/Anvio_pangen/03_PAN ]]
then    
        echo
        echo
        echo -e "The species $REPLY is not availabe :("
        echo -e "Maybe you have forgotten to run Anvio_pangen_wf.sh"
        echo
        exit
fi


WORKDIR=~/PhD/Microbiome_pangenomic_analysis/data/$REPLY/Anvio_pangen
cd $WORKDIR

if [[ ! -d KEGG_metabolism ]]
then    
        mkdir KEGG_metabolism
fi

cd KEGG_metabolism

eval "$(conda shell.bash hook)"
conda activate anvio-7.1
anvi-estimate-metabolism -e $WORKDIR/my-external-genomes_no_bins.txt \
                         --matrix-format \
                         --module-specific-matrices M00569,M00641,M00555,M00542,M00568,M00745 \
                         -O A_indonesiensis

                         