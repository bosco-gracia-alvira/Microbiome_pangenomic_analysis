#!/bin/bash
# This script builds the phylogeny of the species L. plantarum, including my isolates as well as publicly available ones from different isolation sources
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/Lactiplantibacillus_plantarum/Phylogeny"
FASTAS="$WORKDIR/ENAfastas"
METADATA="$WORKDIR/metadata.tsv"

### COMMANDS
IFS=$'\n'

# Check that you have the metadata file
if [[ ! -f "$WORKDIR/metadata.tsv" ]]
then
  echo "You need a metadata table with isolates"
  exit
fi

# Loop through each accession number and download the corresponding files
# while IFS=$'\t' read -r strain accession
# do
#     cp "$FASTAS"/"${accession%0000000}"1.fasta "$FASTAS"/$strain.fa
# done < <(cut -f1,8 "$WORKDIR"/metadata.tsv | grep -v "Strain" | grep -v "GCA_" | sed "s/ //g")

# while IFS=$'\t' read -r strain accession
# do
#     cp "$FASTAS"/"${accession}".fasta "$FASTAS"/$strain.fa
# done < <(cut -f1,8 "$WORKDIR"/metadata.tsv | grep -v "Strain" | grep "GCA_" | sed "s/ //g")

ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO
rm -r ~/Bosco/Microbial_pangenomic_analysis/Lactiplantibacillus_plantarum/Phylogeny
mkdir -p \
    ~/Bosco/Microbial_pangenomic_analysis/Lactiplantibacillus_plantarum/Phylogeny/{genomes,identity,alignment,tree}

FOO

# Copy the ENA genomes to the server
rsync -av \
    "$FASTAS"/*.fa \
    vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Microbial_pangenomic_analysis/Lactiplantibacillus_plantarum/Phylogeny/genomes


# Copy my genomes to the server
cd "$WORKDIR" || exit

while IFS=$'\t' read -r path
do
    rsync -av \
        "$path" \
        vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Microbial_pangenomic_analysis/Lactiplantibacillus_plantarum/Phylogeny/genomes/$(basename ${path%sta})
done < <(cut -f2 "$WORKDIR"/Isolates_path | grep -v "path")

# We connect to the server again to run the first steps of GTDB-Tk
ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO

cd ~/Bosco/Microbial_pangenomic_analysis/Lactiplantibacillus_plantarum/Phylogeny

eval \$(conda shell.bash hook)
conda activate gtdbtk-2.1.1

export GTDBTK_DATA_PATH="/home/vetlinux05/Bosco/db/gtdbtk_r214_database"

# Identify the 
gtdbtk identify \
    -x fa \
    --genome_dir genomes \
    --out_dir identity \
    --cpus 20

gtdbtk align \
    --taxa_filter "g__Lactiplantibacillus" \
    --identify_dir identity \
    --out_dir alignment \
    --cpus 20

gtdbtk infer \
    --msa_file alignment/align/gtdbtk.bac120.msa.fasta.gz \
    --out_dir tree \
    --cpus 20

FOO

# I copy the results and the genomes back to my computer
rsync -av \
    vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Microbial_pangenomic_analysis/Lactiplantibacillus_plantarum/Phylogeny/* \
    "$WORKDIR"
