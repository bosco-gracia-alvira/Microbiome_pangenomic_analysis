#!/bin/bash
# This script builds the pangenome and phylogeny of the L. plantarum collection of isolates from our lab
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/Lactiplantibacillus_plantarum/Pangenomic_Florida"
GENOMES="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/Isolates"
TEMP="/Users/bgracia/PhD_local/temp"

### COMMANDS
IFS=$'\n'

# Check that you have the metadata file
if [[ ! -f "$WORKDIR/metadata.tsv" ]]
then
  echo "You need a metadata table with isolates"
  exit
fi

# Create the folder with genomes
mkdir -p "$WORKDIR/genomes"

# Copy the genomes into the genomes folder
while IFS=$'\t' read -r Strain Fasta rest
do  
    cp "$GENOMES/$Fasta" "$WORKDIR/genomes/$Strain.fa"
done < "$WORKDIR/metadata.tsv"

# In this environment I have the packages to do phylogenies (parsnp, fasttree, harvesttools...)
conda activate phylogeny

# Create the folders that we need
mkdir -p "$TEMP"/{roary,gffs,snp-sites,snp-dists}
cp -r "$WORKDIR/genomes" "$TEMP"

# We want to cluster the genomes based on their functions, for which we will use roary
# First we need the GFF file for each genome
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

# We run snp-dists to calculate the pairwise SNPs in the core genome
snp-dists -m "$TEMP"/roary/core_gene_alignment.aln > "$TEMP"/snp-dists/pairwise_snps_long.tsv
snp-dists "$TEMP"/roary/core_gene_alignment.aln > "$TEMP"/snp-dists/pairwise_snps.tsv

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

# Finally, we use pyANI v0.2.13 to compute the pairwise nucleotide identity between the genomes
average_nucleotide_identity.py \
    -i "$TEMP"/genomes \
    -o "$TEMP"/pyANI \
    --seed 999 \
    -m ANIb

# Move the results to the working directory
mv "$TEMP/roary" "$WORKDIR/"
mv "$TEMP/gffs" "$WORKDIR/"
mv "$TEMP/snp-dists" "$WORKDIR/"
mv "$TEMP/snp-sites" "$WORKDIR/"
mv "$TEMP/pyANI" "$WORKDIR/"
rm -r "${TEMP:?}/*"
