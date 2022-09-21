#!/bin/bash
cd ~/PhD

eval "$(conda shell.bash hook)"
conda activate ragtag-2.1

echo 'Which bacteria do you want to scaffold? Type in the terminal the species with the format: "Genus_species"'
read

SPECIES=$(echo $REPLY | sed 's/_/ /')

if [[ ! -f Microbiome_pangenomic_analysis/data/$REPLY/samples-txt ]]
then
    echo 'You might want to run gen_samples_txt.sh first :/';
    exit
fi

mkdir Microbiome_pangenomic_analysis/data/temp/
WORKDIR=Microbiome_pangenomic_analysis/data/$REPLY
TEMP=Microbiome_pangenomic_analysis/data/temp

cat Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv > $TEMP/taxonomy.tsv
SAMPLES=$(awk -v s="$SPECIES" -F "\t" '$2 ~ s {print $1}' $TEMP/taxonomy.tsv | sort)

ncbi-genome-download bacteria \
    -g "$SPECIES" \
    -s refseq \
    -F fasta \
    -l "all"\
    -P \
    --flat-output \
    -o $TEMP 
gunzip $TEMP/*.gz

echo -e "We have downloaded "$(basename $TEMP/*.fna | wc -l)" genomes from the desired species"

echo 'Which isolate do you want to scaffold? I would not pick any "bin", but it is up to you. These are the available assemblies:'
echo $SAMPLES
read ISOLATE

cp Isolates_assembly/Pool_$(echo $ISOLATE | cut -f1 -d "_")/07.GTDB-Tk/Genomes/$ISOLATE.fa $TEMP

declare -i x=1
for i in $(basename $TEMP/*.fna)
do  ragtag.py scaffold -o $TEMP/multiple_$x $TEMP/$i $TEMP/$ISOLATE.fa;
    x=$x+1;
done
ragtag.py merge $TEMP/$ISOLATE.fa $TEMP/multiple_*/*.agp -o $TEMP/final_multiple

cp $TEMP/final_multiple/ragtag.merge.fasta $WORKDIR/reference.fa

echo -e "name\tpath\n$ISOLATE\t$WORKDIR/reference.fa" > $WORKDIR/fasta-txt

rm -r $TEMP

assembly-stats $WORKDIR/reference.fa