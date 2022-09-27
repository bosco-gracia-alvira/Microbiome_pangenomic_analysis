#!/bin/bash

#echo 'Which species do you want to analyse with Anvio? Type in the terminal the species with the format: "Genus_species"'
#read
REPLY="Leuconostoc_pseudomesenteroides"
SPECIES=$(echo $REPLY | sed 's/_/ /')
cd ~/PhD
mkdir Microbiome_pangenomic_analysis/data/temp/
WORKDIR=~/PhD/Microbiome_pangenomic_analysis/data/$REPLY
TEMP=Microbiome_pangenomic_analysis/data/temp

cat Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv > $TEMP/taxonomy.tsv
SAMPLES=($(awk -v s="$SPECIES" -F "\t" '$2 ~ s {print $1}' $TEMP/taxonomy.tsv | tr "\n" " "))

for i in ${SAMPLES[@]};
do cp Isolates_assembly/Pool_$(echo $i | cut -f1 -d "_")/07.GTDB-Tk/Genomes/$i.fa $TEMP;
done

prokka #Extraction of coding sequences for each isolate

roary #Pan genome and alignment of single copy core genes

snp-sites #Reduction of the alignment only to SNP positions

raxml or mrbayes #Construction of the phylogenetic tree

