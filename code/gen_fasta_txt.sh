#!/bin/bash
cd ~/PhD

eval "$(conda shell.bash hook)"
conda activate ragtag-2.1

echo 'Which bacteria do you want to scaffold? Type in the terminal the species with the format: "Genus_species"'
read

SPECIES=$(echo $REPLY | sed 's/_/ /')

if [[ ! -f Microbiome_pangenomic_analysis/data/$REPLY/Anvio/samples-txt ]]
then
    echo 'You might want to run gen_samples_txt.sh first :/';
    exit
fi

mkdir Microbiome_pangenomic_analysis/data/temp/
WORKDIR=~/PhD/Microbiome_pangenomic_analysis/data/$REPLY/Anvio/
TEMP=Microbiome_pangenomic_analysis/data/temp

cat Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv > $TEMP/taxonomy.tsv
SAMPLES=$(awk -v s="$SPECIES" -F "\t" '$2 ~ s {print $1}' $TEMP/taxonomy.tsv | sort)

ncbi-genome-download bacteria \
        -g "$SPECIES" \
        -s refseq \
        -R "reference,representative" \
        -F fasta \
        -l all \
        -P \
        --flat-output \
        -o $TEMP

gunzip $TEMP/*.gz

echo -e "We have downloaded "$(basename -a $TEMP/*.fna | wc -l)" genomes from the desired species"

echo 'Which isolate do you want to scaffold? I would not pick any "bin", but it is up to you. These are the available assemblies:'
echo $SAMPLES
read ISOLATE

cp Isolates_assembly/Pool_$(echo $ISOLATE | cut -f1 -d "_")/07.GTDB-Tk/Genomes/$ISOLATE.fa $TEMP

ragtag.py scaffold -o $TEMP/Scaffold $TEMP/*.fna $TEMP/$ISOLATE.fa;
cp $TEMP/Scaffold/*.fasta $WORKDIR/reference.fa

echo -e "name\tpath\n$(echo $ISOLATE | cut -f2- -d "_")\t$WORKDIR/reference.fa" > $WORKDIR/fasta-txt

rm -r $TEMP

assembly-stats $WORKDIR/reference.fa

#This part of the script is under maintainance. In theory, it behaves different if you scaffold based on just one reference genome or based on several. But as I use the 
#if [[ $(basename -a $TEMP/*.fna | wc -l) == 1 ]]
#then
#    ragtag.py scaffold -o $TEMP/Scaffold $TEMP/*.fna $TEMP/$ISOLATE.fa;
#    cp $TEMP/Scaffold/*.fasta $WORKDIR/reference.fa
#
#elif [[ $(basename -a $TEMP/*.fna | wc -l) -gt 1 ]]
#then
#    declare -i x=1;
#    for i in $(basename -a $TEMP/*.fna);
#        do  ragtag.py scaffold -o $TEMP/multiple_$x $TEMP/$i $TEMP/$ISOLATE.fa;
#            x=$x+1;
#    done;
#    ragtag.py merge $TEMP/$ISOLATE.fa $TEMP/multiple_*/*.agp -o $TEMP/final_multiple;
#    cp $TEMP/final_multiple/ragtag.merge.fasta $WORKDIR/reference.fa
#elif [[ $(basename -a $TEMP/*.fna | wc -l) == 0 ]]
#then
#    echo "I could not download any reference genome from NCBI, therefore the un-scaffolded assembly will be used as reference";
#    cp $TEMP/$ISOLATE.fa $WORKDIR/reference.fa
#fi