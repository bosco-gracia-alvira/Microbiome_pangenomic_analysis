#!/bin/bash

cd ~/PhD

mkdir -p Microbiome_pangenomic_analysis/data/temp/genomes
TEMP=Microbiome_pangenomic_analysis/data/temp
cat Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv > $TEMP/taxonomy.tsv

echo 'For which of the following species do you want to build a SNP tree? Type it in the terminal with the format: "Genus_species"'
echo 
cut -f2 $TEMP/taxonomy.tsv | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
echo
read

if [[ ! -d Microbiome_pangenomic_analysis/data/$REPLY ]]
then
    mkdir Microbiome_pangenomic_analysis/data/$REPLY
fi

SPECIES=$(echo $REPLY | sed 's/_/ /')

WORKDIR=~/PhD/Microbiome_pangenomic_analysis/data/$REPLY
mkdir $TEMP/snp-sites

cat Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv > $TEMP/taxonomy.tsv
SAMPLES=($(awk -v s="$SPECIES" -F "\t" '$2 ~ s {print $1}' $TEMP/taxonomy.tsv | tr "\n" " "))

echo 'Would you like to include a representative genome as an outgroup in your analysis? Type yes or no.'
read OUTGROUP

if [[ $OUTGROUP = yes ]]
then
        OUT="ref"
        eval "$(conda shell.bash hook)"
        conda activate
        ncbi-genome-download bacteria \
                -g "$SPECIES" \
                -s refseq \
                -R "reference,representative" \
                -F fasta \
                -l all \
                -P \
                --flat-output \
                -o $TEMP/genomes

        gunzip $TEMP/genomes/*.gz
        mv $TEMP/genomes/*.fna $TEMP/genomes/outgroup.fa

        prokka  --outdir $TEMP/prokka \
                --force \
                --prefix outgroup \
                $TEMP/genomes/outgroup.fa 
else
        OUT="noref"
fi


for i in ${SAMPLES[@]};
do  
    cp Isolates_assembly/Pool_$(echo $i | cut -f1 -d "_")/07.GTDB-Tk/Genomes/$i.fa $TEMP/genomes;
    
    prokka  --outdir $TEMP/prokka \
            --force \
            --prefix $i \
            $TEMP/genomes/$i.fa;
done

roary -v -p 8 -e --mafft -f $TEMP/roary $TEMP/prokka/*.gff
snp-sites -mvp -o $TEMP/snp-sites/$REPLY $TEMP/roary/*.aln

mkdir $WORKDIR/Phylogeny_$OUT
mv $TEMP/snp-sites $WORKDIR/Phylogeny_$OUT/
mv $TEMP/roary $WORKDIR/Phylogeny_$OUT/
rm -r $TEMP

eval "$(conda shell.bash hook)"
conda activate base
iqtree -s $WORKDIR/Phylogeny_$OUT/snp-sites/$REPLY.phylip -m GTR+ASC
mkdir $WORKDIR/Phylogeny_$OUT/iqtree
mv $WORKDIR/Phylogeny_$OUT/snp-sites/$REPLY.phylip.* $WORKDIR/Phylogeny_$OUT/iqtree/

#This lines build a phylogenetic tree for each lineage.
#cd ~/PhD/Microbiome_pangenomic_analysis/data/
#for i in $(basename *);
#do      iqtree -s $i/Phylogeny_noref/snp-sites/$i.phylip -m GTR+ASC;
#        mkdir $i/Phylogeny_noref/iqtree;
#        mv $i/Phylogeny_noref/snp-sites/$i.phylip.* $i/Phylogeny_noref/iqtree;
#done