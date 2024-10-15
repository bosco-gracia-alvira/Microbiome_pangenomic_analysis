#!/bin/bash
# This script builds the phylogeny of the species L. plantarum, including my isolates as well as publicly available ones from different isolation sources
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/Lactiplantibacillus_plantarum/Phylogeny"
ISOLATES="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/Isolates"
FASTAS="$WORKDIR/fastas"

### COMMANDS
IFS=$'\n'

# Check that you have the accessions file and the scripts
if [[ ! -f "$WORKDIR"/accessions.txt ]]
then
  echo "You need the list with accessions of publicly available genomes delimited by newlines"
  exit
fi

# Create the subfolders if they don't exist
mkdir -p "$FASTAS"
mkdir -p "$METADATA"

