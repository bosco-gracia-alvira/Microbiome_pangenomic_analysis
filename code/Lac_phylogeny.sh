#!/bin/bash
# This script builds the phylogeny of the species L. plantarum, including my isolates as well as publicly available ones from different isolation sources
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/Lactiplantibacillus_plantarum/Phylogeny"
ISOLATES="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/Isolates"
URL_FASTA="https://www.ebi.ac.uk/ena/browser/api/fasta"
SCRIPTS="$WORKDIR/ENABrowserTools/python3"
FASTAS="$WORKDIR/fastas"
URL_METADATA="https://www.ebi.ac.uk/ena/browser/api/xml"
METADATA="$WORKDIR/metadata"

### COMMANDS
IFS=$'\n'

# Check that you have the accessions file and the scripts
if [[ ! -f "$WORKDIR"/accessions.txt ]]
then
  echo "You need the list with accessions of publicly available genomes delimited by newlines"
  exit
fi

if [[ ! -d "$SCRIPTS" ]]
then
  echo "You need the enaBrowserTools scripts. Get them from:
  https://github.com/enasequence/enaBrowserTools?tab=readme-ov-file"
  exit
fi

# Create the subfolders if they don't exist
mkdir -p "$FASTAS"
mkdir -p "$METADATA"

 
python3 "$SCRIPTS"/enaDataGet.py -f fasta -d "$FASTAS" GCA_000203855


# Loop through each accession number and download the corresponding files
while IFS=$'\t' read -r strain accession; do
    echo "Downloading sequence data for $strain..."
    python3 "$SCRIPTS"/enaDataGet.py -m -f fasta -d "$FASTAS" "${accession}"
done < <(cut -f1,8 "$WORKDIR"/metadata.tsv | grep -v "Strain")

echo "Download completed."
