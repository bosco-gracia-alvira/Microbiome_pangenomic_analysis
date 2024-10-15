#!/bin/bash
# This script builds the phylogeny of the species L. plantarum, including my isolates as well as publicly available ones from different isolation sources
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/Lactiplantibacillus_plantarum/Phylogeny"
ISOLATES="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/Isolates"

### COMMANDS
IFS=$'\n'

# Create the working directory
if [[ ! -f "$WORKDIR"/accessions.txt ]]
then
  echo "You need the list with accessions of publicly available genomes to compare to yours"
  exit
fi












LUXM00000000, LUWN00000000, LUXL00000000, LUXN00000000, LUXO00000000, LTAU00000000, LUWA00000000, LUWB00000000, LUWC00000000, LUWD00000000, LUWE00000000, LUWF00000000, LUWG00000000, LUWH00000000, LUWI00000000, LUWJ00000000, LUWK00000000, LUWL00000000, LUWM00000000, LUWO00000000, LUWP00000000, LUWQ00000000, LUWR00000000, LUWS00000000, LUWT00000000, LUWU00000000, LUWV00000000, LUWW00000000, LUWX00000000, LUWY00000000, LUWZ00000000, LUXA00000000, LUXB00000000, LUXC00000000, LUXD00000000, LUXE00000000, LUXF00000000, LUXG00000000, LUXH00000000, LUXI00000000, LUXJ00000000, LUXK00000000