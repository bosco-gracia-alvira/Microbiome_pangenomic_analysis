#!/bin/bash

cd ~/PhD

echo 'Which species do you want to analyse with Anvio? Type in the terminal the species with the format: "Genus_species"'
echo
cat Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv | cut -f2 | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
read

if [[ ! -d ~/PhD/Microbiome_pangenomic_analysis/data/$REPLY/Anvio ]]
then
        echo -e "The species $REPLY is not availabe :(";
        echo -e "Maybe you have forgotten to run gen_samples_txt.sh"
        exit
fi

if [[ ! -f ~/PhD/Microbiome_pangenomic_analysis/data/$REPLY/Anvio/fasta-txt ]]
then
        echo -e "The species $REPLY is not availabe :(";
        echo -e "Maybe you have forgotten to run gen_fasta_txt.sh"
        exit
fi

cd ~/PhD/Microbiome_pangenomic_analysis/data/$REPLY/Anvio
CONTIG=$(cat fasta-txt | awk '!/name/ {print $1}')


eval "$(conda shell.bash hook)"
conda activate anvio-7.1

anvi-run-workflow -w metagenomics \
        -c ../../../code/config.json

anvi-run-pfams -c 03_CONTIGS/$CONTIG-contigs.db -T 4
anvi-run-ncbi-cogs -c 03_CONTIGS/$CONTIG-contigs.db -T 4
anvi-run-kegg-kofams -c 03_CONTIGS/$CONTIG-contigs.db -T 4

anvi-script-add-default-collection -p 06_MERGED/$CONTIG/PROFILE.db

anvi-gen-variability-profile -c 03_CONTIGS/$CONTIG-contigs.db \
                             -p 06_MERGED/$CONTIG/PROFILE.db \
                             -C DEFAULT \
                             -b EVERYTHING \
                             --engine CDN \
                             --include-site-pnps \
                             --kiefl-mode \
                             -o $REPLY-SCVs.txt

anvi-gen-variability-profile -c 03_CONTIGS/$CONTIG-contigs.db \
                             -p 06_MERGED/$CONTIG/PROFILE.db \
                             -C DEFAULT \
                             -b EVERYTHING \
                             --engine AA \
                             --kiefl-mode \
                             -o $REPLY-SAAVs.txt


awk '{ if ((NR == 1) || ($73!=$74)) { print } }' $REPLY-SCVs.txt > consensus-SCVs.tsv

awk '{ if ((NR == 1) || ($30!=$31)) { print } }' $REPLY-SAAVs.txt > consensus-SAAVs.tsv