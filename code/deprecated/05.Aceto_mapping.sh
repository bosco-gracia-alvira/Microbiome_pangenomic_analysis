#!/bin/bash
# This script tests if Acetobacter reads map ambiguously between species
# Bosco Gracia Alvira, 2024

### VARIABLES

# Ask the user which species to analyse

# Set the paths
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data"
LOCATION_COLD="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads_cold"
LOCATION_HOT="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/poolseq_reads_hot"

# First argument of the script is the species to analyse
REPLY="Acetobacter_indonesiensis"

if [[ -z "$REPLY" ]]
then
    echo "You need to provide the species you want to analyse as first argument"
    echo
    echo "These are the species that are available for the pangenomic analysis"
    echo "Copy the name of the species you want to analyse and paste it as first argument of the script"
    echo
    cut -f2 "$WORKDIR"/taxonomy.tsv | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
    echo
    exit
elif [[ "$REPLY" == "h" ]]
then
    echo "These are the species that are available for the pangenomic analysis"
    echo "Copy the name of the species you want to analyse and paste it as first argument of the script"
    echo
    cut -f2 "$WORKDIR"/taxonomy.tsv | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
    echo
    exit
fi

SPECIES=$(echo "$REPLY" | sed 's/_/ /')
SAMPLES=$(awk -v s="$SPECIES" -F "\t" '$8 ~ s {print $1}' "$WORKDIR"/Genome_metadata.tsv | grep -v "user" | sed 's/-/_/g')

REFERENCE="$WORKDIR/${REPLY}/Anvio_popgen/${REPLY}_pangenome.fasta"
RESULTS="$WORKDIR/${REPLY}/mapping_tests"
RAW_READS="$SNPS/test_reads"
LOGS="$SNPS/test_logs"
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

cp "$REFERENCE" "$RESULTS"/ref.fa
samtools faidx "$RESULTS"/ref.fa

# Link the isolates reads to the working

# S131 A. indonesiensis
ln -s "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly/Pool_503/02.Rm_adapters/fastq_clean/S131.clean_1.fq.gz" \
    "$RAW_READS/S131_1.fq.gz"
ln -s "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly/Pool_503/02.Rm_adapters/fastq_clean/S131.clean_2.fq.gz" \
    "$RAW_READS/S131_2.fq.gz"

# S133 A. indonesiensis
ln -s "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly/Pool_503/02.Rm_adapters/fastq_clean/S133.clean_1.fq.gz" "$RAW_READS/S133_1.fq.gz"
ln -s "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly/Pool_503/02.Rm_adapters/fastq_clean/S133.clean_2.fq.gz" "$RAW_READS/S133_2.fq.gz"

# X715 A. sicerae
ln -s "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly/Pool_591/02.Rm_adapters/fastq_clean/X715.clean_1.fq.gz" "$RAW_READS/X715_1.fq.gz"
ln -s "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly/Pool_591/02.Rm_adapters/fastq_clean/X715.clean_2.fq.gz" "$RAW_READS/X715_2.fq.gz"

# X302 A. oryzifermentans
ln -s "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly/Pool_591/02.Rm_adapters/fastq_clean/X302.clean_1.fq.gz" "$RAW_READS/X302_1.fq.gz"
ln -s "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly/Pool_591/02.Rm_adapters/fastq_clean/X302.clean_2.fq.gz" "$RAW_READS/X302_2.fq.gz"

# X311 A. malorum
ln -s "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly/Pool_591/02.Rm_adapters/fastq_clean/X311.clean_1.fq.gz" "$RAW_READS/X311_1.fq.gz"
ln -s "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly/Pool_591/02.Rm_adapters/fastq_clean/X311.clean_2.fq.gz" "$RAW_READS/X311_2.fq.gz"

# Map the reads to the reference pangenome and calculate the % identity of the reads
for i in "$RAW_READS"/*_1.fq.gz
do
    name=$(basename "$i" | rev | cut -d "/" -f1 | rev | cut -d "." -f1)
    anir.rb -r "${i}" -g "$REFERENCE" -i 0.7 \
        -L "$SNPS/identity_distribution"/${name}.list \
        -H "$SNPS/identity_distribution"/${name}.hist \
        -T "$SNPS/identity_distribution"/${name}.tsv \
        --threads 16
done




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
    mean_coverage=$(samtools depth "$i" | awk '{sum+=$3} END {if (NR>0) print sum/NR; else print 0}')
    if (( $(echo "$mean_coverage > 10" | bc -l) )); then
        echo "Keeping $bam with mean coverage $mean_coverage"
    else
        echo "Discarding $bam with mean coverage $mean_coverage"
        rm "$bam"
    fi
done

mkdir -p "$SNPS/identity_distribution"
for i in "$BAMS"/*_sorted.bam; do
name=$(basename "$i" | rev | cut -d "/" -f1 | rev | cut -d "." -f1)
anir.rb -m "$i" --m-format bam --threads 16 -H "$SNPS/identity_distribution/${name}.hist" \

done


# This command calculates allele frequencies in all the positions, then does gene calling and finally changes the sample IDs, that originally include the whole path to the bam file and the .bam extension
bcftools mpileup -Ou -f "$SNPS"/ref.fa -b "$SNPS/coverage_5.txt" --threads 16 | \
    bcftools call  --ploidy 1 -Ou -mv | bcftools view -i 'QUAL>20 && DP>5' - > "$SNPS/$REPLY.vcf"

# bcftools mpileup -Ou -f "$SNPS"/ref.fa -b "$SNPS/coverage_5.txt" -d 5 --threads 16 | \
#     bcftools call  --threads 8 --ploidy 1 -Ou -mv | bcftools view -i 'QUAL>20 && DP>5' - > "$SNPS/${REPLY}_5x.vcf"

# bcftools mpileup -Ou -f "$SNPS"/ref.fa -b "$SNPS/coverage_10.txt" -d 10 --threads 16 | \
#     bcftools call  --threads 8 --ploidy 1 -Ou -mv | bcftools view -i 'QUAL>20 && DP>10' - > "$SNPS/${REPLY}_10x.vcf"

# We make a BED file that is used by PLINK to compute the PCA of the samples based on SNPs frequency
plink2 --vcf "$SNPS/${REPLY}.vcf" --double-id --allow-extra-chr --make-bed --out "$SNPS/${REPLY}"
plink2 --bfile "$SNPS/${REPLY}" --double-id --allow-extra-chr --pca --out "$SNPS/${REPLY}"
plink2 --bfile "$SNPS/${REPLY}" --read-freq "$SNPS/${REPLY}.afreq" --score "$SNPS/${REPLY}.eigenvec.var" 2 3 header-read --out "$SNPS/${REPLY}_loadings"

# # We make a BED file that is used by PLINK to compute the PCA of the samples based on SNPs frequency
# plink2 --vcf "$SNPS/${REPLY}_5x.vcf" --double-id --allow-extra-chr --make-bed --out "$SNPS/${REPLY}_5x"
# plink2 --bfile "$SNPS/${REPLY}_5x" --double-id --allow-extra-chr --pca --out "$SNPS/${REPLY}_5x"

# # We make a BED file that is used by PLINK to compute the PCA of the samples based on SNPs frequency
# plink2 --vcf "$SNPS/${REPLY}_10x.vcf" --double-id --allow-extra-chr --make-bed --out "$SNPS/${REPLY}_10x"
# plink2 --bfile "$SNPS/${REPLY}_10x" --double-id --allow-extra-chr --pca --out "$SNPS/${REPLY}_10x"


# vcftools --vcf "$SNPS/$REPLY.filtered.vcf" --window-pi 100 --out "$SNPS/${REPLY}_100bp"

# vcftools --vcf "$SNPS/$REPLY.filtered.vcf" --TajimaD 100 --out "$SNPS/${REPLY}_100bp"

# This script plots the PCA of the samples based on the SNPs frequency
Rscript "$WORKDIR"/../code/05.Plot_PCA.R "$REPLY"
