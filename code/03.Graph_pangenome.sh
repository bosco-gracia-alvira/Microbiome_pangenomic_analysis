#!/bin/bash
# This script creates a graph pangenome from all the isolates of a species
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data"
ISOLATES="$WORKDIR/Isolates"

### COMMANDS
IFS="
"

# Create the SuperPang directory
if [[ ! -f "$WORKDIR"/"$REPLY"/SuperPang ]]
then
    mkdir -p "$WORKDIR"/"$REPLY"/SuperPang
fi

# Ask the user which species to analyse
echo 'Which species do you want to analyse with Anvio? Type in the terminal the species with the format: "Genus_species"'
echo
cut -f2 "$WORKDIR"/taxonomy.tsv | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
read

# Set the variable SPECIES from the input and pick the samples that belong to the species of interest
SPECIES=$(echo "$REPLY" | sed 's/_/ /')
SAMPLES=$(awk -v s="$SPECIES" -F "\t" '$2 ~ s {print $1}' "$WORKDIR"/taxonomy.tsv | grep -v "user")

# Create the folder for the species
ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO

if [[ ! -d ~/Bosco/Microbial_pangenomic_analysis/"$REPLY"/SuperPang/Genomes ]]
then  
    mkdir -p ~/Bosco/Microbial_pangenomic_analysis/"$REPLY"/SuperPang/Genomes
fi

FOO

# Copy the genomes to the VetLinux server
for i in $SAMPLES
do
        scp -r "$ISOLATES/$i.fasta" \
        vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Microbial_pangenomic_analysis/"$REPLY"/SuperPang/Genomes
done


# In this long EOF I run the whole SuperPang pipeline remotely in VetLinux
ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO

eval "\$(conda shell.bash hook)"
conda activate SuperPang-0.9

cd ~/Bosco/Microbial_pangenomic_analysis/"$REPLY"/SuperPang/Genomes

SuperPang.py \
        --fasta ~/Bosco/Microbial_pangenomic_analysis/"$REPLY"/SuperPang/Genomes/*.fasta \
        --assume-complete \
        --output-dir ~/Bosco/Microbial_pangenomic_analysis/"$REPLY"/SuperPang/output \
        --force-overwrite \
        --threads 24

rm -r Genomes

FOO

# I copy the results back to my computer
scp -r vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Microbial_pangenomic_analysis/"$REPLY"/SuperPang/output/* \
"$WORKDIR/$REPLY/SuperPang"
