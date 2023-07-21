#!/bin/bash
# This script is used to estimate the metabolism of A. indonesiensis pangenome using Anvi'o. I am following the instructions of the Anvi'o tutorial:
# https://merenlab.org/tutorials/infant-gut/#chapter-v-metabolism-prediction
# Bosco Gracia Alvira, 2023

### VARIABLES
START="/Users/bgracia/Dropbox (PopGen)/Bosco/PhD_Dropbox"

### COMMANDS
eval "$(conda shell.bash hook)"
conda activate anvio-7.1

cd "$START" || exit

echo 'To which species do you want to generate a misc-table? Type in the terminal the species with the format: "Genus_species"'
echo
cat "$START"/Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv | cut -f2 | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
read

WORKDIR="$START"/Microbiome_pangenomic_analysis/data/$REPLY/Anvio_pangen

if [[ ! -d "$WORKDIR"/03_PAN ]]
then    
        echo
        echo
        echo -e "The species $REPLY is not availabe :("
        echo -e "Maybe you have forgotten to run Anvio_pangen_wf.sh"
        echo
        exit
else
        mkdir -p "$WORKDIR"/05_ENRICHMENT
fi


for i in {COG20_FUNCTION,COG20_CATEGORY,KOfam,KEGG_Module,COG20_PATHWAY,KEGG_Class}
do      anvi-compute-functional-enrichment-across-genomes -e "$WORKDIR"/my-external-genomes_no_bins.txt \
                                                  -o "$WORKDIR"/05_ENRICHMENT/Enrichment"${i}".txt \
                                                  -G "$WORKDIR"/Groups.txt \
                                                  --annotation-source "${i}"
done

for i in $(cut -f6 "$WORKDIR"/05_ENRICHMENT/EnrichmentCOG20_FUNCTION.txt)
do grep "$i" "$WORKDIR"/COG20.tsv
done | sort > "$WORKDIR"/COG20-2.tsv

for i in $(cut -f6 "$WORKDIR"/05_ENRICHMENT/EnrichmentCOG20_FUNCTION.txt)
do echo "${i}";
done | sort

awk -F, 'FNR==NR {lines[$2]; next} $2 in lines' file1 file2

#anvi-compute-functional-enrichment-in-pan -p $WORKDIR/03_PAN/MYPAN-PAN.db \
                                          #-g $WORKDIR/03_PAN/MYPAN-GENOMES.db \
                                          #--category-variable clade \
                                          #--annotation-source $WORKDIR/COG20_FUNCTION \
                                          #-o $WORKDIR/05_ENRICHMENT/functional-enrichment.txt