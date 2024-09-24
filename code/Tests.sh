
for i in Acetobacter_indonesiensis Acetobacter_malorum Acetobacter_oryzifermentans Acetobacter_persici Acetobacter_sicerae Lacticaseibacillus_paracasei Lactiplantibacillus_plantarum Leuconostoc_pseudomesenteroides Levilactobacillus_brevis Morganella_morganii_B
do
    REPLY="$i"
    echo "$REPLY"
    cd "$REPLY/SNPs_analysis" || continue
    echo "Working in $REPLY"
    SNPS="."
    rm "$SNPS/coverage.txt"
    for j in $(basename -a reads/*_1.fq.gz)
    do
        sample="${j%_1.fq.gz}"
        mean_coverage=$(
            samtools depth "bams/${sample}_sorted.bam" | \
            sort -k3,3nr | \
            awk '{ all[NR] = $3; sum+=$3 } END {if (NR==0) print 0; else { for (j=int(NR*0.1)+1; j<=int(NR*0.9); j++) s+=all[j]; print s/(int(NR*0.9)-int(NR*0.1)) } }')
        # Append the mean coverage to the coverage file
        echo -e "${sample}\t${mean_coverage}" >> "$SNPS/coverage.txt"
    done

    awk '$2 >= 10 {print "./bams/"$1"_sorted.bam"}' "$SNPS/coverage.txt" | grep -v "sample" > "$SNPS/coverage_10.txt"

    awk '$2 >= 5 {print "./bams/"$1"_sorted.bam"}' "$SNPS/coverage.txt" | grep -v "sample" > "$SNPS/coverage_5.txt"

    # This chunk counts the reference and alternative allele frequency in each position and in each sample (mpileup), then calls the SNPs (call) and filters the SNPs (no indels) with a quality above 20 and a depth above 5 
    bcftools mpileup -f "$SNPS"/ref.fa -b "$SNPS/coverage_5.txt" -Q 20 -D -d 50 -a DP,AD,QS,SCR -Ou --threads 16 | \
        bcftools call  --ploidy 1 -Ou -cv --threads 16 | \
        bcftools view -i 'QUAL > 20' -v snps -m2 -M2 -Ov - > "$SNPS/temp_$REPLY.vcf"
        #bcftools view -e 'FORMAT/DP[:0] < 3' 

    # We reformat the headers of the VCF file, that by default include the relative path to the bams
    bcftools view -h "$SNPS/temp_$REPLY.vcf" > "$SNPS/headers.txt"
    sed -i '' 's|./bams/||g; s|_sorted.bam||g' "$SNPS/headers.txt"
    bcftools reheader -h "$SNPS/headers.txt" -o "$SNPS/$REPLY.vcf" "$SNPS/temp_$REPLY.vcf"

    # We extract the frequency of the SNPs in each sample
    bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "$SNPS/$REPLY.vcf" > "$SNPS/$REPLY.freq"

    cd ../../
done



for i in Lactiplantibacillus_plantarum Acetobacter_indonesiensis Acetobacter_oryzifermentans Lacticaseibacillus_paracasei Leuconostoc_pseudomesenteroides Levilactobacillus_brevis  
do
    if [[ $i == "Lactiplantibacillus_plantarum" ]]
    then 
        echo "Making $i graph pangenome"
        REPLY="$i"
        ./code/03.Graph_pangenome.sh "$REPLY"
    else
        echo "Whole pipeline for taxon $i"
        REPLY="$i"
        rm -r data/"$REPLY"/Anvio_pangen
        ./code/01.Gen_tables.sh "$REPLY"
        ./code/02.Anvio_pangen_wf.sh "$REPLY"
        echo "Finished pangenomic analysis for taxon $i"
        ./code/03.Graph_pangenome.sh "$REPLY"
        echo "Finished graph pangenome for taxon $i"
    fi
done
