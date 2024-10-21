#!/bin/bash
# This script builds the phylogeny of the species L. plantarum, including my isolates as well as publicly available ones from different isolation sources
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/Lactiplantibacillus_plantarum/Phylogeny"
FASTAS="$WORKDIR/ENAfastas"
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
mkdir -p "$TEMP"/parsnp
cp -r "$WORKDIR/genomes" "$TEMP"

# Download the reference genome in GenBank format using wget
wget -P "$TEMP/parsnp" \
    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/203/855/GCA_000203855.3_ASM20385v3/GCA_000203855.3_ASM20385v3_genomic.gbff.gz
# Unzip the downloaded file
gunzip \
    --to-stdout "$TEMP/parsnp/GCA_000203855.3_ASM20385v3_genomic.gbff.gz" > "$TEMP/parsnp/Lpla_ref.gbff"
rm "$TEMP/parsnp/GCA_000203855.3_ASM20385v3_genomic.gbff.gz"

wget -P "$TEMP/parsnp" \
    ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/203/855/GCA_000203855.3_ASM20385v3/GCA_000203855.3_ASM20385v3_genomic.fna.gz
# Unzip the downloaded file
gunzip \
    --to-stdout "$TEMP/parsnp/GCA_000203855.3_ASM20385v3_genomic.fna.gz" > "$TEMP/parsnp/Lpla_ref.fna"
rm "$TEMP/parsnp/GCA_000203855.3_ASM20385v3_genomic.fna.gz"

# Run parsnp, that extracts SNPs common to all the genomes
parsnp \
    -g "$TEMP/parsnp/Lpla_ref.gbff" \
    -d "$TEMP/genomes" \
    -o "$TEMP/parsnp" \
    -p 12

mv "$TEMP/parsnp" "$WORKDIR/"
rm -r "${TEMP:?}/*"

# # Loop through each accession number and download the corresponding files
# # while IFS=$'\t' read -r strain accession
# # do
# #     cp "$FASTAS"/"${accession%0000000}"1.fasta "$FASTAS"/$strain.fa
# # done < <(cut -f1,8 "$WORKDIR"/metadata.tsv | grep -v "Strain" | grep -v "GCA_" | sed "s/ //g")

# # while IFS=$'\t' read -r strain accession
# # do
# #     cp "$FASTAS"/"${accession}".fasta "$FASTAS"/$strain.fa
# # done < <(cut -f1,8 "$WORKDIR"/metadata.tsv | grep -v "Strain" | grep "GCA_" | sed "s/ //g")

# ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO
# rm -r ~/Bosco/Microbial_pangenomic_analysis/Lactiplantibacillus_plantarum/Phylogeny
# mkdir -p \
#     ~/Bosco/Microbial_pangenomic_analysis/Lactiplantibacillus_plantarum/Phylogeny/{genomes,identity,alignment,tree}

# FOO

# # Copy the ENA genomes to the server
# rsync -av \
#     "$FASTAS"/*.fa \
#     vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Microbial_pangenomic_analysis/Lactiplantibacillus_plantarum/Phylogeny/genomes


# # Copy my genomes to the server
# cd "$WORKDIR" || exit

# while IFS=$'\t' read -r path
# do
#     rsync -av \
#         "$path" \
#         vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Microbial_pangenomic_analysis/Lactiplantibacillus_plantarum/Phylogeny/genomes/$(basename ${path%sta})
# done < <(cut -f2 "$WORKDIR"/Isolates_path | grep -v "path")

# # We connect to the server again to run the first steps of GTDB-Tk
# ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO

# cd ~/Bosco/Microbial_pangenomic_analysis/Lactiplantibacillus_plantarum/Phylogeny

# eval \$(conda shell.bash hook)
# conda activate gtdbtk-2.1.1

# export GTDBTK_DATA_PATH="/home/vetlinux05/Bosco/db/gtdbtk_r214_database"

# # Identify the 
# gtdbtk identify \
#     -x fa \
#     --genome_dir genomes \
#     --out_dir identity \
#     --cpus 20

# gtdbtk align \
#     --taxa_filter "s__Lactiplantibacillus paraplantarum" \
#     --identify_dir identity \
#     --out_dir alignment \
#     --cpus 20

# gtdbtk infer \
#     --msa_file alignment/align/gtdbtk.bac120.msa.fasta.gz \
#     --out_dir tree \
#     --cpus 20

# FOO

# # I copy the results and the genomes back to my computer
# rsync -av \
#     vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Microbial_pangenomic_analysis/Lactiplantibacillus_plantarum/Phylogeny/* \
#     "$WORKDIR"

# # I can also do the phylogeny myself
# conda activate phylogeny

# # Create a new folder
# mkdir "$WORKDIR"/tree_manual

# gunzip -c "$WORKDIR"/alignment/align/gtdbtk.bac120.msa.fasta.gz \
#     > "$WORKDIR"/tree_manual/gtdbtk.bac120.msa.fasta



# # Run model test to fit the best model for the phylogeny
# modeltest-ng \
#     -d aa \
#     -i "$WORKDIR"/tree_manual/gtdbtk.bac120.msa.fasta \
#     -o "$WORKDIR"/tree_manual/test \
#     -p 12 \
#     -t ml

# # There is a software that does multilocus phylogeny of different strains
# mkdir "$WORKDIR"/parsnp

# cp -r "$WORKDIR"/genomes "$WORKDIR"

# # Download the reference genome in GenBank format using wget
# wget -P "$WORKDIR/parsnp" \
#     ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/203/855/GCA_000203855.3_ASM20385v3/GCA_000203855.3_ASM20385v3_genomic.gbk.gz
# # Unzip the downloaded file
# gunzip \
#     --to-stdout "$WORKDIR/parsnp/GCA_000203855.3_ASM20385v3_genomic.gbk.gz" > "$WORKDIR/parsnp/Lpla_ref.gbk"
# rm "$WORKDIR/parsnp/GCA_000203855.3_ASM20385v3_genomic.gbk.gz"

# wget -P "$WORKDIR/parsnp" \
#     ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/203/855/GCA_000203855.3_ASM20385v3/GCA_000203855.3_ASM20385v3_genomic.fna.gz
# # Unzip the downloaded file
# gunzip \
#     --to-stdout "$WORKDIR/parsnp/GCA_000203855.3_ASM20385v3_genomic.fna.gz" > "$WORKDIR/parsnp/Lpla_ref.fna"
# rm "$WORKDIR/parsnp/GCA_000203855.3_ASM20385v3_genomic.fna.gz"




