#!/bin/bash

cd ~/PhD
1=Acetobacter_oryzifermentans
SPECIES=$(echo $1 | sed 's/_/ /')

mkdir -p Microbiome_pangenomic_analysis/data/$1
mkdir Microbiome_pangenomic_analysis/data/temp

TEMP=Microbiome_pangenomic_analysis/data/temp/

cat Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv > $TEMP/taxonomy.tsv
cat Isolates_assembly/Pool_???/01.Raw_data/Demultiplexed/reads.names > $TEMP/names.csv

SAMPLES=$(awk -v s=$SPECIES -F '\t' '$2 ~ s {print $1}' $TEMP/taxonomy.tsv | cut -f1-6 -d "_" | sort | uniq)
SAMPLES2=($(echo $SAMPLES | tr "\n" " "))

echo "sample\tr1\tr2" > $TEMP/header.tmp

for i in ${SAMPLES2[*]};
do awk -v i=$i -F ',' '$2==i {print $2,"Isolates_assembly/Pool_"substr($2,1,3)"/02.Rm_adapters/fastq/"$1"_1.fq.gz","Isolates_assembly/Pool_"substr($2,1,3)"/02.Rm_adapters/fastq/"$1"_2.fq.gz"}' $TEMP/names.csv >> $TEMP/body.tmp;
done

cat $TEMP/header.tmp $TEMP/body.tmp > $TEMP/samples-txt

