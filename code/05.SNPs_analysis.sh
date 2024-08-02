#!/bin/bash
# This script runs Anvio population genomics analysis for the desired species
# Bosco Gracia Alvira, 2024

### VARIABLES

# Ask the user which species to analyse

# Set the paths
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data"
LOCATION_COLD="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads_cold"
LOCATION_HOT="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads_hot"

echo 'Which species do you want to analyse with Anvio? Type in the terminal the species with the format: "Genus_species"'
echo
cut -f2 "$WORKDIR"/taxonomy.tsv | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
read
SAMPLES=$(awk -v s="$SPECIES" -F "\t" '$2 ~ s {print $1}' "$WORKDIR"/taxonomy.tsv | grep -v "user" | sed 's/-/_/g')

REFERENCE="$WORKDIR/${REPLY}/Anvio_popgen/${REPLY}_pangenome.fasta"
SNPS="$WORKDIR/${REPLY}/SNPs_analysis"
BAMS="$WORKDIR/${REPLY}/Anvio_popgen/04_MAPPING/${REPLY}"

### COMMANDS
IFS="
"

if [[ ! -d "$WORKDIR/${REPLY}/Anvio_popgen/04_MAPPING/${REPLY}" ]]
then
        echo -e "The species $REPLY is not availabe :("
        echo -e "Maybe you have forgotten to run 04.Anvio_popgen_wf.sh"
        exit
fi

if [[ ! -f "$WORKDIR"/samples-txt ]]
then
        echo -e "The species $REPLY is not availabe :("
        echo -e "I don't find the samples-txt file"
        echo -e "You need to run 01.Gen_tables.sh"
        exit
fi

if [[ ! -f "$SNPS" ]]
then
    mkdir -p "$SNPS"
fi

cp "$REFERENCE" "$SNPS"/ref.fa
samtools faidx "$SNPS"/ref.fa

# This command calculates allele frequencies in all the positions, then does gene calling and finally changes the sample IDs, that originally include the whole path to the bam file and the .bam extension
bcftools mpileup -Ou -f "$SNPS"/ref.fa "$BAMS"/*.bam | \
    bcftools call  --ploidy 1 -Ou -mv  | sed -e "s|${BAMS}/||g" | sed -e "s|.bam||g" > "$SNPS/$REPLY.vcf"

# This command filters the SNPs by quality and depth
bcftools view -i 'QUAL>20 && DP>20' "$SNPS/$REPLY.vcf" > "$SNPS/$REPLY.filtered.vcf"

# We make a BED file
plink2 --vcf "$SNPS/$REPLY.filtered.vcf" --double-id --allow-extra-chr --make-bed --out "$SNPS/$REPLY"
plink2 --bfile "$SNPS/$REPLY" --double-id --allow-extra-chr --pca --out "$SNPS/$REPLY"


vcftools --vcf "$SNPS/$REPLY.filtered.vcf" --window-pi 100 --out "$SNPS/${REPLY}_100bp"

vcftools --vcf "$SNPS/$REPLY.filtered.vcf" --TajimaD 100 --out "$SNPS/${REPLY}_100bp"

