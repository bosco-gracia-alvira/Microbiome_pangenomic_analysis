#!/bin/bash

cd ~/PhD

mkdir -p Microbiome_pangenomic_analysis/data/temp/genomes
TEMP=Microbiome_pangenomic_analysis/data/temp
cat Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv > $TEMP/taxonomy.tsv

echo 'Which of the following species do you want to analyse? Type it in the terminal with the format: "Genus_species"'
echo 
cut -f2 $TEMP/taxonomy.tsv | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
echo
read

SPECIES=$(echo $REPLY | sed 's/_/ /')

if [[ ! -d Microbiome_pangenomic_analysis/data/$REPLY ]]
then
    mkdir Microbiome_pangenomic_analysis/data/$REPLY
fi

WORKDIR=Microbiome_pangenomic_analysis/data/$REPLY

SAMPLES=($(awk -v s="$SPECIES" -F "\t" '$2 ~ s {print $1}' $TEMP/taxonomy.tsv))

for i in ${SAMPLES[@]};
do cp Isolates_assembly/Pool_$(echo $i | cut -f1 -d "_")/07.GTDB-Tk/Genomes/$i.fa $TEMP/genomes;
done

eval "$(conda shell.bash hook)"
conda activate pyani_env

average_nucleotide_identity.py -v -i $TEMP/genomes -m ANIb -o $TEMP/ANI -g

if [[ ! -d $WORKDIR/ANI ]]
then
    mkdir $WORKDIR/ANI
fi

mv $TEMP/ANI/* $WORKDIR/ANI/
rm -r $TEMP