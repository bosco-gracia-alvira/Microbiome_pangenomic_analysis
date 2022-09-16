#!/bin/bash

cd ~/PhD

echo 'Which bacteria do you want to analyse? Type it in the terminal with the format: Genus_species'
read

SPECIES=$(echo $REPLY | sed 's/_/ /')

mkdir Microbiome_pangenomic_analysis/data/$REPLY
mkdir Microbiome_pangenomic_analysis/data/temp

WORKDIR=Microbiome_pangenomic_analysis/data/$REPLY
TEMP=Microbiome_pangenomic_analysis/data/temp

cat Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv > $TEMP/taxonomy.tsv
cat Isolates_assembly/Pool_???/01.Raw_data/Demultiplexed/reads.names > $TEMP/names.csv

SAMPLES=($(awk -v s="$SPECIES" -F "\t" '$2 ~ s {print $1}' $TEMP/taxonomy.tsv | cut -f1-6 -d "_" | sort | uniq | tr "\n" " "))

echo -e 'sample\tr1\tr2' > $TEMP/header.tmp
touch $TEMP/body.tmp
for i in ${SAMPLES[*]};
do awk -F "," -v i="$i" '$2==i {print $2,"Isolates_assembly/Pool_"substr($2,1,3)"/02.Rm_adapters/fastq/"$1"_1.fq.gz","Isolates_assembly/Pool_"substr($2,1,3)"/02.Rm_adapters/fastq/"$1"_2.fq.gz"}' $TEMP/names.csv >> $TEMP/body.tmp;
done

cat $TEMP/header.tmp $TEMP/body.tmp > $WORKDIR/samples-txt

rm -r $TEMP
