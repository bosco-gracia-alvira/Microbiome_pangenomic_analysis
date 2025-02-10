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

# Create the parsnp folder
mkdir -p "$TEMP"/{parsnp,gffs}
cp -r "$WORKDIR/genomes" "$TEMP"

# Download the reference genome in GenBank format using wget
wget -P "$TEMP" \
    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/203/855/GCA_000203855.3_ASM20385v3/GCA_000203855.3_ASM20385v3_genomic.gbff.gz

# Unzip the downloaded file
gunzip \
    --to-stdout "$TEMP/GCA_000203855.3_ASM20385v3_genomic.gbff.gz" > "$TEMP/Lpla_ref.gbff"
rm "$TEMP/GCA_000203855.3_ASM20385v3_genomic.gbff.gz"

# Download the reference genome in fna formats using wget
wget -P "$TEMP" \
    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/203/855/GCA_000203855.3_ASM20385v3/GCA_000203855.3_ASM20385v3_genomic.fna.gz

# Unzip the downloaded file
gunzip \
    --to-stdout "$TEMP/GCA_000203855.3_ASM20385v3_genomic.fna.gz" > "$TEMP/Lpla_ref.fna"
rm "$TEMP/GCA_000203855.3_ASM20385v3_genomic.fna.gz"

# Run parsnp, that extracts SNPs common to all the genomes
parsnp \
    -g "$TEMP/Lpla_ref.gbff" \
    -d "$TEMP/genomes" \
    -o "$TEMP/parsnp" \
    -p 12

# Build the phylogenetic tree using iqtree with 1000 bootstrap replicates and the GTR+ASC model
cd $TEMP/snp-sites || exit
iqtree -T 8 -s $TEMP/snp-sites/Lpla.phylip --boot 1000 -m MF

# Now we want to cluster the genomes based on their functions, for which we will use roary
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

roary   -e -n \
        -p 16 \
        -f "$TEMP"/roary/ \
        "$TEMP"/gffs/*.gff

mv "$TEMP/parsnp" "$WORKDIR/"
mv "$TEMP/roary" "$WORKDIR/"
mv "$TEMP/gffs" "$WORKDIR/"
rm -r "${TEMP:?}/*"
