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

  
for i in {COG20_FUNCTION,COG20_CATEGORY,KOfam,KEGG_Module,COG20_PATHWAY,KEGG_Class};
do      anvi-compute-functional-enrichment-across-genomes -e $WORKDIR/my-external-genomes_no_bins.txt \
                                                  -o $WORKDIR/Enrichment${i}.txt \
                                                  -G $WORKDIR/Groups.txt \
                                                  --annotation-source ${i};
done

for i in $(cut -f6 EnrichmentCOG20_FUNCTION.txt);
do grep "$i" COG20.tsv;
done | sort > COG20-2.tsv

for i in $(cut -f6 EnrichmentCOG20_FUNCTION.txt);
do echo $i;
done | sort

awk -F, 'FNR==NR {lines[$2]; next} $2 in lines' file1 file2

#anvi-compute-functional-enrichment-in-pan -p $WORKDIR/Anvio_pangen/03_PAN/MYPAN-PAN.db \
                                          #-g $WORKDIR/Anvio_pangen/03_PAN/MYPAN-GENOMES.db \
                                          #--category-variable clade \
                                          #--annotation-source COG20_FUNCTION \
                                          #-o $WORKDIR/Anvio_pangen/functional-enrichment.txt