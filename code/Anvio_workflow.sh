#!/bin/bash

echo 'Which species do you want to analyse with Anvio? Type in the terminal the species with the format: "Genus_species"'
read

if [[ ! -d ~/PhD/Microbiome_pangenomic_analysis/data/$REPLY ]]
then
        echo -e "The species $REPLY is not availabe :(";
        echo -e "Maybe you have forgotten to run gen_samples_txt.sh"
        exit
fi

if [[ ! -f ~/PhD/Microbiome_pangenomic_analysis/data/$REPLY/fasta-txt ]]
then
        echo -e "The species $REPLY is not availabe :(";
        echo -e "Maybe you have forgotten to run gen_fasta_txt.sh"
        exit
fi

cd ~/PhD/Microbiome_pangenomic_analysis/data/$REPLY
CONTIG=$(cat fasta-txt | awk '!/name/ {print $1}')


eval "$(conda shell.bash hook)"
conda activate anvio-7.1

anvi-run-workflow -w metagenomics \
        -c ../../code/config.json

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

awk 'NR > 1 && $77!=0 && $78==1 {print $0}' $REPLY-SCVs.txt > filtered-SCVs.tsv

awk 'NR > 1 && $34!=0 && $35==1 {print $0}' $REPLY-SAAVs.txt > filtered-SAAVs.tsv