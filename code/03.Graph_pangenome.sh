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

# First argument of the script is the species to analyse
REPLY=$1
REPLY="Morganella_morganii_B"
if [[ -z "$REPLY" ]]
then
    echo "You need to provide the species you want to analyse as first argument"
    echo
    echo "These are the species that are available for the pangenomic analysis"
    echo "Copy the name of the species you want to analyse and paste it as first argument of the script"
    echo
    cut -f2 "$WORKDIR"/taxonomy.tsv | awk -F 's__' '{print $2}' | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
    echo
    exit
elif [[ "$REPLY" == "h" ]]
then
    echo "These are the species that are available for the pangenomic analysis"
    echo "Copy the name of the species you want to analyse and paste it as first argument of the script"
    echo
    cut -f2 "$WORKDIR"/taxonomy.tsv | awk -F 's__' '{print $2}' | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
    echo
    exit
fi

# Set the variable SPECIES from the input and pick the samples that belong to the species of interest
SPECIES=$(echo "$REPLY" | sed 's/_/ /')
SAMPLES=$(awk -v s="$SPECIES" -F "\t" '$2 ~ s {print $1}' "$WORKDIR"/taxonomy.tsv | grep -v "user")

species_count=$(cut -f2 "$WORKDIR"/taxonomy.tsv | awk -F 's__' '{print $2}' | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column)
count=$(echo "$species_count" | grep -w "$REPLY" | awk '{print $1}')

# Check if the count is equal to 1
if [ "$count" -eq 1 ]; then
  echo "There is only one genome available for the species $REPLY." There is no need to make a graph pangenome.
  exit
else
  echo "There are $count genomes available for the species $REPLY." Thus, we will make a pangenome.
fi

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
