#!/bin/bash
# This script runs Anvio population genomics analysis for the desired species
# Bosco Gracia Alvira, 2024

### VARIABLES

# Ask the user which species to analyse

# Set the paths
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data"
LOCATION_COLD="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads_cold"
LOCATION_HOT="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads_hot"

echo 'Which species do you want to analyse with Anvio? Type in the terminal the species with the format: "Genus_species"'
echo
cut -f2 "$WORKDIR"/taxonomy.tsv | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
read
SPECIES=$(echo "$REPLY" | sed 's/_/ /')
SAMPLES=$(awk -v s="$SPECIES" -F "\t" '$8 ~ s {print $1}' "$WORKDIR"/Genome_metadata.tsv | grep -v "user" | sed 's/-/_/g')

REFERENCE="$WORKDIR/${REPLY}/Anvio_popgen/${REPLY}_pangenome.fasta"
SNPS="$WORKDIR/${REPLY}/SNPs_analysis"
RAW_READS="$SNPS/reads"
LOGS="$SNPS/logs"
BAMS="$SNPS/bams"

### COMMANDS
IFS="
"

if [[ ! -f "$REFERENCE" ]]
then
        echo -e "The species $REPLY is not availabe :("
        echo -e "Maybe you have forgotten to run 04.Anvio_popgen_wf.sh"
        exit
fi

if [[ ! -f "$RAW_READS" ]]
then
    mkdir -p "$RAW_READS"
fi

if [[ ! -f "$BAMS" ]]
then
    mkdir -p "$BAMS"
fi

if [[ ! -f "$LOGS" ]]
then
    mkdir -p "$LOGS"
fi

cp "$REFERENCE" "$SNPS"/ref.fa
samtools faidx "$SNPS"/ref.fa

# Link the poolseq reads sets to the working directory
for i in $(basename "$LOCATION_HOT"/F*)
do
  ln -s "$LOCATION_HOT"/F*/Pooled_${i}_Clean_noCont_1.fq.gz "$RAW_READS"/h${i}_1.fq.gz
  ln -s "$LOCATION_HOT"/F*/Pooled_${i}_Clean_noCont_2.fq.gz "$RAW_READS"/h${i}_2.fq.gz
done

for i in $(basename "$LOCATION_COLD"/F*)
do
  ln -s "$LOCATION_COLD"/F*/Pooled_${i}_Clean_noCont_1.fq.gz "$RAW_READS"/c${i}_1.fq.gz
  ln -s "$LOCATION_COLD"/F*/Pooled_${i}_Clean_noCont_2.fq.gz "$RAW_READS"/c${i}_2.fq.gz
done

# Link the isolates reads to the working

# Create a file linking the samples of interest and the location of their reads

# This script creates a file that links the genomes to their reads path (absolute_path_reads.txt)
Rscript "$WORKDIR"/../code/05.Reformat_SNP.R

head -n 1 "$WORKDIR"/absolute_path_reads.txt > "$SNPS/genomes_reads.txt"

for i in $SAMPLES
do
        grep -w "$i" "$WORKDIR"/absolute_path_reads.txt >> "$SNPS/genomes_reads.txt"
done

# Iterate over the files to link the reads to the working directory
while IFS=$'\t' read -r sample r1 r2
do
        # Skip the header line
        if [[ "$sample" != "sample" ]]
        then
            ln -s "$r1" "$RAW_READS/${sample}_1.fq.gz"
            ln -s "$r2" "$RAW_READS/${sample}_2.fq.gz"
        fi
done < "$SNPS/genomes_reads.txt"


# Map the reads to the reference pangenome

# Create an index for the reference combined genome
bowtie2-build --threads 16 "$SNPS"/ref.fa "$SNPS"/ref

echo -e "sample\tcoverage" > "$SNPS"/coverage.txt

# Competitive mapping mapping against each of the reads sets
for i in $(basename -a "$RAW_READS"/*_1.fq.gz)
do
  sample="${i%_1.fq.gz}"
  r1="${i}"
  r2="${sample}_2.fq.gz"

  echo "Mapping ${sample}"
  # Map paired end reads using bowtie with stringent settings and output the result to a sam file
  bowtie2 \
    -x "$SNPS"/ref \
    -q --very-sensitive \
    --no-mixed \
    --no-discordant \
    -1 "$RAW_READS/${r1}" \
    -2 "$RAW_READS/${r2}" \
    -S "$BAMS/${sample}.sam" \
    --threads 16 \
    --rg-id "${sample}" \
    --rg "SM:${sample}" > "$LOGS/bowtie2_${sample}.log" 2>&1

  # Turn the sam into bam and sort it
    samtools view \
        -bS \
        -@ 16 \
        "$BAMS/${sample}.sam" |\
    samtools sort \
        -@ 16 \
        -O bam \
        -o "$BAMS/${sample}_sorted.bam" \
        -

    # Delete the sam to save space
    rm "$BAMS/${sample}.sam"

    # Calculate the mean coverage of the sample
    mean_coverage=$(samtools depth "$BAMS/${sample}_sorted.bam" | awk '{sum+=$3} END {if (NR>0) print sum/NR; else print 0}')

    # Append the mean coverage to the coverage file
    echo -e "${sample}\t${mean_coverage}" >> "$SNPS/coverage.txt"
done

cd "$SNPS" || exit

# List the samples with coverage above 10 and 5
awk '$2 >= 10 {print "./bams/"$1"_sorted.bam"}' "$SNPS/coverage.txt" | grep -v "sample" > "$SNPS/coverage_10.txt"
awk '$2 >= 5 {print "./bams/"$1"_sorted.bam"}' "$SNPS/coverage.txt" | grep -v "sample" > "$SNPS/coverage_5.txt"

# Remove the bams whose coverage is below 10
for i in "$BAMS"/*_sorted.bam; do
    mean_coverage=$(samtools depth "$bam" | awk '{sum+=$3} END {if (NR>0) print sum/NR; else print 0}')
    if (( $(echo "$mean_coverage > 10" | bc -l) )); then
        echo "Keeping $bam with mean coverage $mean_coverage"
    else
        echo "Discarding $bam with mean coverage $mean_coverage"
        rm "$bam"
    fi
done

# This command calculates allele frequencies in all the positions, then does gene calling and finally changes the sample IDs, that originally include the whole path to the bam file and the .bam extension
bcftools mpileup -Ou -f "$SNPS"/ref.fa -b "$SNPS/coverage_5.txt" --threads 8 | \
    bcftools call  --ploidy 1 -Ou -mv > "$SNPS/${REPLY}.vcf"
bcftools view -i 'QUAL>20 && DP>5' "$SNPS/$REPLY.vcf" > "$SNPS/$REPLY.filtered.vcf"

bcftools mpileup -Ou -f "$SNPS"/ref.fa -b "$SNPS/coverage_5.txt" -d 5 --threads 8 | \
    bcftools call  --threads 8 --ploidy 1 -Ou -mv | bcftools view -i 'QUAL>20 && DP>5' - > "$SNPS/${REPLY}_5x.vcf"

bcftools mpileup -Ou -f "$SNPS"/ref.fa -b "$SNPS/coverage_10.txt" -d 10 --threads 8 | \
    bcftools call  --threads 8 --ploidy 1 -Ou -mv | bcftools view -i 'QUAL>20 && DP>10' - > "$SNPS/${REPLY}_10x.vcf"

# We make a BED file that is used by PLINK to compute the PCA of the samples based on SNPs frequency
plink2 --vcf "$SNPS/${REPLY}.filtered.vcf" --double-id --allow-extra-chr --make-bed --out "$SNPS/${REPLY}"
plink2 --bfile "$SNPS/${REPLY}" --double-id --allow-extra-chr --pca --out "$SNPS/${REPLY}"
plink2 --bfile "$SNPS/${REPLY}" --read-freq "$SNPS/${REPLY}.afreq" --score "$SNPS/${REPLY}.eigenvec.var" 2 3 header-read --out "$SNPS/${REPLY}_loadings"

# We make a BED file that is used by PLINK to compute the PCA of the samples based on SNPs frequency
plink2 --vcf "$SNPS/${REPLY}_5x.vcf" --double-id --allow-extra-chr --make-bed --out "$SNPS/${REPLY}_5x"
plink2 --bfile "$SNPS/${REPLY}_5x" --double-id --allow-extra-chr --pca --out "$SNPS/${REPLY}_5x"

# We make a BED file that is used by PLINK to compute the PCA of the samples based on SNPs frequency
plink2 --vcf "$SNPS/${REPLY}_10x.vcf" --double-id --allow-extra-chr --make-bed --out "$SNPS/${REPLY}_10x"
plink2 --bfile "$SNPS/${REPLY}_10x" --double-id --allow-extra-chr --pca --out "$SNPS/${REPLY}_10x"


vcftools --vcf "$SNPS/$REPLY.filtered.vcf" --window-pi 100 --out "$SNPS/${REPLY}_100bp"

vcftools --vcf "$SNPS/$REPLY.filtered.vcf" --TajimaD 100 --out "$SNPS/${REPLY}_100bp"
