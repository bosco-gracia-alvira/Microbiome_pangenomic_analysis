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

if [[ ! -d Microbiome_pangenomic_analysis/data/$REPLY/Anvio ]]
then
    mkdir Microbiome_pangenomic_analysis/data/$REPLY/Anvio
fi
mkdir Microbiome_pangenomic_analysis/data/temp

WORKDIR=Microbiome_pangenomic_analysis/data/$REPLY/Anvio

cat Isolates_assembly/Pool_???/01.Raw_data/Demultiplexed/reads.names > $TEMP/names.csv

SAMPLES=($(awk -v s="$SPECIES" -F "\t" '$2 ~ s {print $1}' $TEMP/taxonomy.tsv | cut -f1-6 -d "_" | sort | uniq | tr "\n" " "))

echo -e 'sample\tr1\tr2' > $TEMP/header.tmp
touch $TEMP/body.tmp
for i in ${SAMPLES[@]};
do  awk -F "," -v i="$i" '$2==i {print $2,"../../../../Isolates_assembly/Pool_"substr($2,1,3)"/02.Rm_adapters/fastq_clean/"$1".clean_1.fq.gz","../../../../Isolates_assembly/Pool_"substr($2,1,3)"/02.Rm_adapters/fastq_clean/"$1".clean_2.fq.gz"}' $TEMP/names.csv >> $TEMP/body.tmp;
done
cut -d "_" -f2- $TEMP/body.tmp | sort | awk 'FNR==1{
              print
              next
            }
            {
              A[$1]=$1 in A ? A[$1]","$2:$2
            }
            {
              B[$1]=$1 in B ? B[$1]","$3:$3
            }
         END{
              for(i in A)
                   print i,A[i],B[i]
            }' | tr " " "\t" >> $TEMP/body_merged.tmp # This strange awk command that I don't fully understand merges together rows by sample name and merges together the location of the r1 and r2 with commas.

cat $TEMP/header.tmp $TEMP/body_merged.tmp > $WORKDIR/samples-txt

rm -r $TEMP