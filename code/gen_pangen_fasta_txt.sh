#!/bin/bash
cd ~/PhD

mkdir -p Microbiome_pangenomic_analysis/data/temp
TEMP=Microbiome_pangenomic_analysis/data/temp
cat Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv > $TEMP/taxonomy.tsv

echo 'Which of the following species do you want to analyse? Type it in the terminal with the format: "Genus_species"'
echo 
cut -f2 $TEMP/taxonomy.tsv | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
echo
read

SPECIES=$(echo $REPLY | sed 's/_/ /')

if [[ ! -d Microbiome_pangenomic_analysis/data/$REPLY/Anvio_pangen ]]
then
    mkdir -p Microbiome_pangenomic_analysis/data/$REPLY/Anvio_pangen;
fi

WORKDIR=~/PhD/Microbiome_pangenomic_analysis/data/$REPLY/Anvio_pangen

SAMPLES=$(awk -v s="$SPECIES" -F "\t" '$2 ~ s {print $1}' $TEMP/taxonomy.tsv | sort)
awk -v s="$SPECIES" -F "\t" '$2 ~ s {print $1}' $TEMP/taxonomy.tsv | sort > $TEMP/samples.tmp
cut -f2- -d "_" $TEMP/samples.tmp > $TEMP/name.tmp

for i in $(cat $TEMP/samples.tmp);
do echo -e ~/PhD/Isolates_assembly/Pool_${i:0:3}/07.GTDB-Tk/Genomes/$i.fa >> $TEMP/path.tmp;
done

for i in $(cat $TEMP/name.tmp);
do echo -e 02_CONTIGS/${i}-contigs.db >> $TEMP/contigs_db_path.tmp;
done

echo -e 'name\tpath' > $TEMP/header.tmp
paste $TEMP/name.tmp $TEMP/path.tmp > $TEMP/body.tmp
cat $TEMP/header.tmp $TEMP/body.tmp > $WORKDIR/fasta-txt

echo -e 'name\tcontigs_db_path' > $TEMP/header2.tmp
paste $TEMP/name.tmp $TEMP/contigs_db_path.tmp > $TEMP/body2.tmp
cat $TEMP/header2.tmp $TEMP/body2.tmp > $WORKDIR/my-external-genomes.txt

rm -r $TEMP