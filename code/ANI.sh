#!/bin/bash

cd ~/PhD

echo 'Which bacteria do you want to analyse? Type it in the terminal with the format: "Genus_species"'
read

SPECIES=$(echo $REPLY | sed 's/_/ /')

if [[ ! -d Microbiome_pangenomic_analysis/data/$REPLY ]]
then
    mkdir Microbiome_pangenomic_analysis/data/$REPLY
fi
mkdir -p Microbiome_pangenomic_analysis/data/temp/genomes

WORKDIR=Microbiome_pangenomic_analysis/data/$REPLY
TEMP=Microbiome_pangenomic_analysis/data/temp

cat Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv > $TEMP/taxonomy.tsv

SAMPLES=($(awk -v s="$SPECIES" -F "\t" '$2 ~ s {print $1}' $TEMP/taxonomy.tsv))

for i in ${SAMPLES[@]};
do cp Isolates_assembly/Pool_$(echo $i | cut -f1 -d "_")/07.GTDB-Tk/Genomes/$i.fa $TEMP/genomes;
done

eval "$(conda shell.bash hook)"
conda activate pyani_env

average_nucleotide_identity.py -v -i $TEMP/genomes -m ANIb -o $TEMP/ANI -g
mv $TEMP/ANI/ANIb_percentage_identity.png $WORKDIR

rm -r $TEMP