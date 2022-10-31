#!/bin/bash

#Clade specific phylogeny

cd ~/PhD

mkdir -p Microbiome_pangenomic_analysis/data/temp/genomes{1,2}
TEMP=Microbiome_pangenomic_analysis/data/temp
cat Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv > $TEMP/taxonomy.tsv

echo 'For which of the following species do you want to build a SNP tree? Type it in the terminal with the format: "Genus_species"'
echo 
cut -f2 $TEMP/taxonomy.tsv | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
echo
read
REPLY="Acetobacter_indonesiensis"
SPECIES=$(echo $REPLY | sed 's/_/ /')

WORKDIR=~/PhD/Microbiome_pangenomic_analysis/data/$REPLY
mkdir $TEMP/snp-sites{1,2}


CLADE1=($(awk '$2==0 {print $1}' $WORKDIR/PopCOGenT/$REPLY.cluster.tsv))
CLADE2=($(awk '$2==1 {print $1}' $WORKDIR/PopCOGenT/$REPLY.cluster.tsv))

for i in ${CLADE1[@]};
do  
    cp Isolates_assembly/Pool_$(echo $i | cut -f1 -d "_")/07.GTDB-Tk/Genomes/$i.fa $TEMP/genomes1;
    
    prokka  --outdir $TEMP/prokka1 \
            --force \
            --prefix $i \
            $TEMP/genomes1/$i.fa;
done

roary -v -p 8 -e --mafft -f $TEMP/roary1 $TEMP/prokka1/*.gff
snp-sites -mvp -o $TEMP/snp-sites1/$REPLY $TEMP/roary1/*.aln

mkdir $WORKDIR/Phylogeny_Clade1
mv $TEMP/snp-sites1 $WORKDIR/Phylogeny_Clade1/
mv $TEMP/roary1 $WORKDIR/Phylogeny_Clade1/

for i in ${CLADE2[@]};
do  
    cp Isolates_assembly/Pool_$(echo $i | cut -f1 -d "_")/07.GTDB-Tk/Genomes/$i.fa $TEMP/genomes2;
    
    prokka  --outdir $TEMP/prokka2 \
            --force \
            --prefix $i \
            $TEMP/genomes2/$i.fa;
done

roary -v -p 8 -e --mafft -f $TEMP/roary2 $TEMP/prokka2/*.gff
snp-sites -mvp -o $TEMP/snp-sites2/$REPLY $TEMP/roary2/*.aln

mkdir $WORKDIR/Phylogeny_Clade2
mv $TEMP/snp-sites2 $WORKDIR/Phylogeny_Clade2/
mv $TEMP/roary2 $WORKDIR/Phylogeny_Clade2/


eval "$(conda shell.bash hook)"
conda activate base
for i in {1,2};
do      iqtree -s $WORKDIR/Phylogeny_Clade$i/snp-sites$i/$REPLY.phylip -m GTR+ASC;
        mkdir $WORKDIR/Phylogeny_Clade$i/iqtree$i;
        mv $WORKDIR/Phylogeny_Clade$i/snp-sites$i/$REPLY.phylip.* $WORKDIR/Phylogeny_Clade$i/iqtree$i/;
done
