#!/bin/bash
# This script builds the pangenome of the four Lactobacillus isolates that we have in order to find specific genes that are either temperature-specific or isolate-specific. These genes will be used to design specific primers for each isolate.
# 
# Bosco Gracia Alvira, 2024

### VARIABLES
WORKDIR="/Users/bgracia/Dropbox (PopGen)/Bosco/PhD_Dropbox/Amplicon_primer_design/Lac_specific_primers"
TEMP="/Users/bgracia/PhD_local/temp"


### COMMANDS
conda activate

if [[ -d "$WORKDIR"/00_GENOMES ]]
then  
  mkdir -p \
    "$WORKDIR"/01_PROKKA \
    "$WORKDIR"/02_ROARY \
    "$WORKDIR"/03_PPANGGOLIN \
    "$WORKDIR"/04_SUBSETS \
    "$WORKDIR"/05_PRIMERS
fi


# We call prokka to annotate the genomes
for i in $(basename "$WORKDIR"/00_GENOMES/*.fasta)
do
    prokka  \
                --outdir "$WORKDIR"/01_PROKKA \
                --force \
                --prefix ${i%.fasta} \
                "$WORKDIR"/00_GENOMES/$i \
                --cpus 16
done

cat "$WORKDIR"/01_PROKKA/*.faa > "$WORKDIR"/04_SUBSETS/proteins.faa
cat "$WORKDIR"/01_PROKKA/*.ffn > "$WORKDIR"/04_SUBSETS/genes.ffn
cat "$WORKDIR"/01_PROKKA/*.tsv > "$WORKDIR"/04_SUBSETS/annotations.tsv

# We call ppanggolin to build the pangenome
conda activate ppanggolin

rm "$WORKDIR"/03_PPANGGOLIN/genomes.txt

for i in $(basename "$WORKDIR"/01_PROKKA/*.gff)
do
    echo -e "${i%.gff}\t01_PROKKA/${i}" >> "$WORKDIR"/03_PPANGGOLIN/genomes.txt
done

ppanggolin all --anno "$WORKDIR"/03_PPANGGOLIN/genomes.txt --kingdom bacteria --output "$WORKDIR"/03_PPANGGOLIN -f --cpu 16








