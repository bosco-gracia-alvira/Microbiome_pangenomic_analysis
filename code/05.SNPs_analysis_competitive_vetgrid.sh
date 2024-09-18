#!/bin/bash
# This script maps the reads to the species reference pangenome and calculates the SNPs
# It uses reads that were previously mapped competitively to the whole microbiome
# Bosco Gracia Alvira, 2024

### VARIABLES

# Ask the user which species to analyse

# Set the paths
WORKDIR="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data"
LOCAL="/Volumes/Temp/competitive_mapping"
READS="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Ancestral_microbiome/data/Competitive_mapping_microbiome/mapped"
BBMAP="/Volumes/Data/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly/Scripts/bbmap"
BBMAP_LOCAL="$LOCAL/bbmap"
# First argument of the script is the species to analyse
REPLY=$1

if [[ -z "$REPLY" ]]
then
    echo "You need to provide the species you want to analyse as first argument"
    echo
    echo "These are the species that are available for the pangenomic analysis"
    echo "Copy the name of the species you want to analyse and paste it as first argument of the script"
    echo
    cut -f2 "$WORKDIR"/taxonomy.tsv | awk -F 's__' '{print $2}' | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
    echo
    exit
elif [[ "$REPLY" == "h" ]]
then
    echo "These are the species that are available for the pangenomic analysis"
    echo "Copy the name of the species you want to analyse and paste it as first argument of the script"
    echo
    cut -f2 "$WORKDIR"/taxonomy.tsv | awk -F 's__' '{print $2}' | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
    echo
    exit
fi

SPECIES=$(echo "$REPLY" | sed 's/_/ /')
SAMPLES=$(awk -v s="$SPECIES" -F "\t" '$8 ~ s {print $1}' "$WORKDIR"/Genome_metadata.tsv | grep -v "user" | sed 's/-/_/g')

# Evaluate the number of genomes available for the species
species_count=$(cut -f2 "$WORKDIR"/taxonomy.tsv | awk -F 's__' '{print $2}' | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r)
count=$(echo "$species_count" | grep -w "$REPLY" | awk '{print $1}')
echo "$count genomes available for taxon $SPECIES"

# If there is only one genome available, we will use it as reference
if [ "$count" -eq 1 ]; then
  echo "There is only one genome available for the species $REPLY." This genome will be used as reference.
  genome=$(grep "$SPECIES" "$WORKDIR"/Genome_metadata.tsv | awk -F "\t" '{print $1}')
  REFERENCE="$WORKDIR/Isolates/${genome}.fasta"
else
  echo "There are $count genomes available for the species $REPLY." Thus, we will use the pangenome as reference.
  REFERENCE="$WORKDIR/${REPLY}/SuperPang/assembly.fasta"
fi

if [[ ! -d "$LOCAL" ]]
then  
    mkdir -p "$LOCAL"
fi

# Copy BBMAP to the temporary local directory
cp -r "$BBMAP" "$LOCAL"

SNPS="$LOCAL/${REPLY}/SNPs_competitive"
RAW_READS="$SNPS/reads"
LOGS="$SNPS/logs"
BAMS="$SNPS/bams"

### COMMANDS
IFS="
"

if [[ ! -f "$REFERENCE" ]]
then
        echo -e "The species $REPLY is not availabe :("
        echo -e "Maybe you have forgotten to run 03.Graph_pangenome.sh"
        exit
fi

if [[ ! -d "$READS" ]]
then
        echo -e "The folder with bams does not seem to exist"
        echo -e "Maybe you have forgotten to run the competitive_mapping script"
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

# We reformat the reference genome to match the names used in the Anvio_popgen workflow
eval "$(conda shell.bash hook)"
conda activate anvio-7.1

anvi-script-reformat-fasta \
        -l 250 \
        -o "$SNPS"/ref.fa \
        --simplify-names \
        "$REFERENCE"

samtools faidx "$SNPS"/ref.fa
conda deactivate

# Extract the reads from the bam files
# From the hot pools
for i in "$READS"/h*
do
    sample=$(basename "${i}")
    echo "Extracting reads from ${sample}"
    samtools fastq \
        -F 4 \
        -1 "$RAW_READS"/"${sample}"_1.fq.gz \
        -2 "$RAW_READS"/"${sample}"_2.fq.gz \
        "${READS}"/"${sample}"/*"${REPLY}"*.bam
done

# From the cold pools
for i in "$READS"/c*
do
    sample=$(basename "${i}")
    echo "Extracting reads from ${sample}"
    samtools fastq \
        -F 4 \
        -1 "$RAW_READS"/"${sample}"_1.fq.gz \
        -2 "$RAW_READS"/"${sample}"_2.fq.gz \
        "${READS}"/"${sample}"/*"${REPLY}"*.bam
done

# From the isolates
for i in $SAMPLES
do
    sample="i${i}"
    echo "Extracting reads from ${sample}"
    samtools fastq \
            -F 4 \
            -1 "$RAW_READS"/"${i}"_1.fq.gz \
            -2 "$RAW_READS"/"${i}"_2.fq.gz \
            "${READS}"/"${sample}"/*"${REPLY}"*.bam
done

# Map the reads to the reference pangenome

# Create an index for the reference combined genome
"$BBMAP_LOCAL"/bbmap.sh ref="$SNPS"/ref.fa path="$SNPS" -Xmx24g

echo -e "sample\tcoverage" > "$SNPS"/coverage.txt

# Competitive mapping mapping against each of the reads sets
for i in $(basename -a "$RAW_READS"/*_1.fq.gz)
do
  sample="${i%_1.fq.gz}"
  r1="${i}"
  r2="${sample}_2.fq.gz"

  echo "Mapping ${sample}"
  # Map paired end reads using bowtie with stringent settings and output the result to a sam file
  "$BBMAP_LOCAL"/bbmap.sh \
    ref="$SNPS"/ref.fa \
    path="$SNPS" \
    in="$RAW_READS/${r1}" \
    in2="$RAW_READS/${r2}" \
    out="$BAMS/${sample}.sam" \
    pairedonly=t \
    killbadpairs=t \
    rgid="${sample}" \
    threads=auto \
    idtag=t \
    minid=0.95 > "$LOGS/bbmap_${name}.log" 2>&1

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

    # Calculate the mean coverage of the sample (Truncated coverage 80%)
    mean_coverage=$(
        samtools depth "$BAMS/${sample}_sorted.bam" | \
        sort -k3,3nr | \
        awk '{ all[NR] = $3; sum+=$3 } END {if (NR==0) print 0; else { for (i=int(NR*0.1)+1; i<=int(NR*0.9); i++) s+=all[i]; print s/(int(NR*0.9)-int(NR*0.1)) } }')
        #samtools depth "$BAMS/${sample}_sorted.bam" | awk '{sum+=$3} END {if (NR>0) print sum/NR; else print 0}'
    # Append the mean coverage to the coverage file
    echo -e "${sample}\t${mean_coverage}" >> "$SNPS/coverage.txt"
done

cd "$SNPS" || exit

# List the samples with coverage above 10 and 5
awk '$2 >= 10 {print "./bams/"$1"_sorted.bam"}' "$SNPS/coverage.txt" | grep -v "sample" > "$SNPS/coverage_10.txt"
awk '$2 >= 5 {print "./bams/"$1"_sorted.bam"}' "$SNPS/coverage.txt" | grep -v "sample" > "$SNPS/coverage_5.txt"

# # Remove the bams whose coverage is below 10
# for i in "$BAMS"/*_sorted.bam; do
#     mean_coverage=$(samtools depth "$i" | awk '{sum+=$3} END {if (NR>0) print sum/NR; else print 0}')
#     if (( $(echo "$mean_coverage >= 5" | bc -l) )); then
#         echo "Keeping $i with mean coverage $mean_coverage"
#     else
#         echo "Discarding $i with mean coverage $mean_coverage"
#         rm "$i"
#     fi
# done

# This chunk counts the reference and alternative allele frequency in each position and in each sample (mpileup), then calls the SNPs (call) and filters the SNPs (no indels) with a quality above 20 and a depth above 5 
bcftools mpileup -f "$SNPS"/ref.fa -b "$SNPS/coverage_5.txt" -Q 20 -D -d 50 -a DP,AD,QS,SCR -Ou --threads 16 | \
    bcftools call  --ploidy 1 -Ou -cv --threads 16 | \
    bcftools view -i 'QUAL > 20' -v snps -m2 -M2 -Ov - > "$SNPS/temp_$REPLY.vcf"
    # bcftools view -e 'FORMAT/DP[:0] < 3' -Ov - Filtering can be done in the R analysis step

# We reformat the headers of the VCF file, that by default include the relative path to the bams
bcftools view -h "$SNPS/temp_$REPLY.vcf" > "$SNPS/headers.txt"
sed -i '' 's|./bams/||g; s|_sorted.bam||g' "$SNPS/headers.txt"
bcftools reheader -h "$SNPS/headers.txt" -o "$SNPS/$REPLY.vcf" "$SNPS/temp_$REPLY.vcf"

# We extract the frequency of the SNPs in each sample
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "$SNPS/$REPLY.vcf" > "$SNPS/$REPLY.freq"

# We make a BED file that is used by PLINK to compute the PCA of the samples based on SNPs frequency
plink2 --vcf "$SNPS/${REPLY}.vcf" --double-id --allow-extra-chr --make-bed --out "$SNPS/${REPLY}"
plink2 --bfile "$SNPS/${REPLY}" --double-id --allow-extra-chr --pca --out "$SNPS/${REPLY}"

rm "$SNPS/temp_$REPLY.vcf"

# Move the results to the working directory and remove the temporary directory
rsync -av --remove-source-files "$SNPS" "$WORKDIR/${REPLY}/"
rm -r "${LOCAL:?}/${REPLY:?}"

# This script plots the PCA of the samples based on the SNPs frequency
Rscript "$WORKDIR"/../code/05.Plot_PCA.R "$REPLY"
