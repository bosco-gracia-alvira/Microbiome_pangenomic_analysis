#!/bin/bash
# This script builds the phylogeny of the species L. plantarum, including my isolates as well as publicly available ones from different isolation sources
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/Lactiplantibacillus_plantarum/Phylogeny"
TEMP="/Users/bgracia/PhD_local/temp"

### COMMANDS
IFS=$'\n'

# Check that you have the metadata file
if [[ ! -f "$WORKDIR/metadata.tsv" ]]
then
  echo "You need a metadata table with isolates"
  exit
elif [[ ! -d "$WORKDIR/genomes" ]]
then
  echo "You need a folder with the genomes that you want to compare"
  exit
fi

# In this environment I have the packages to do phylogenies (parsnp, fasttree, harvesttools...)
conda activate phylogeny

# Create the folders that we need
mkdir -p "$TEMP"/{roary,gffs,snp-sites}
cp -r "$WORKDIR/genomes" "$TEMP"

# We want to cluster the genomes based on their functions, for which we will use roary
# Get GFF files  for each gene
for i in "$TEMP"/genomes/*.fa
do
    # Get the base name of the file (without extension)
    base=$(basename "$i" .fa)
    
    # Run Prokka to generate the GFF file
    prokka --outdir "$TEMP/gffs/$base" --prefix "$base" "$i"
    mv "$TEMP/gffs/$base/$base.gff" "$TEMP/gffs/"
    rm -r "$TEMP/gffs/$base"
done

# Run roary to build a pangenome and move the results out of the _* folder
roary   -e -n \
        -p 16 \
        -f "$TEMP"/roary/ \
        "$TEMP"/gffs/*.gff

mv "$TEMP"/roary/_*/* "$TEMP"/roary/
rm -r "$TEMP"/roary/_*

# Run snp-sites to extract only the SNPs
snp-sites \
    -mvp \
    -o $TEMP/snp-sites/Lpla \
    $TEMP/roary/*.aln

# Build the phylogenetic tree using iqtree with 1000 bootstrap replicates and the GTR+ASC model
cd $TEMP/snp-sites || exit
iqtree \
    -s $TEMP/snp-sites/Lpla.phylip \
    --boot 1000 \
    -m GTR+ASC \
    -T 8

# Move the results to the working directory
mv "$TEMP/roary" "$WORKDIR/"
mv "$TEMP/gffs" "$WORKDIR/"
mv "$TEMP/snp-sites" "$WORKDIR/"
rm -r "${TEMP:?}/*"
