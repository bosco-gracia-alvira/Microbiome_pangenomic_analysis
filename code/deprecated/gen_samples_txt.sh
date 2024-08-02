#!/bin/bash
# This script creates the file samples-txt for Anvio.
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data"

### COMMANDS
IFS="
"

echo 'Which species do you want to analyse with Anvio? Type in the terminal the species with the format: "Genus_species"'
echo
cut -f2 "$WORKDIR"/taxonomy.tsv | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
read

SPECIES=$(echo $REPLY | sed 's/_/ /')

if [[ ! -f "$WORKDIR"/$REPLY/Anvio/fasta-txt ]]
then    
        echo
        echo
        echo -e "The species $REPLY is not availabe :("
        echo -e "Maybe you have forgotten to run 01.gen_fasta_txt.sh"
        echo
        exit
fi

# Pick the samples that belong to the species of interest
SAMPLES=$(awk -v s="$SPECIES" -F "\t" '$2 ~ s {print $1}' "$WORKDIR"/taxonomy.tsv | grep -v "user")

# Create the sample-txt file
echo -e 'sample\tr1\tr2' > "$WORKDIR"/"$REPLY"/Anvio/sample-txt

# Iterate through the samples to populate the fasta-txt file
# Note that the sample names cannot have hyphens because Anvi'o does not like it
# Also, I use relative paths because Anvi'o does not like spaces in the absolute pahts
for i in $SAMPLES
do
    # Stupid Anvio does not like hyphens, so I replace them with underscores
    name=$(echo "${i}" | sed "s/-/\_/g")
    echo -e "${name}\t../../Isolates/${i}.fasta" >> "$WORKDIR"/"$REPLY"/Anvio/fasta-txt
done




























mkdir -p Microbiome_pangenomic_analysis/data/temp

cat Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv > $TEMP/taxonomy.tsv

echo 'Which of the following species do you want to analyse? Type it in the terminal with the format: "Genus_species"'
echo 
cut -f2 $TEMP/taxonomy.tsv | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
echo
read





mkdir Microbiome_pangenomic_analysis/data/temp

WORKDIR=Microbiome_pangenomic_analysis/data/$REPLY/Anvio

cat Isolates_assembly/Pool_???/01.Raw_data/Demultiplexed/reads.names > $TEMP/names.csv

SAMPLES=($(awk -v s="$SPECIES" -F "\t" '$2 ~ s {print $1}' $TEMP/taxonomy.tsv | cut -f1-6 -d "_" | sort | uniq | tr "\n" " "))

echo -e 'sample\tr1\tr2' > $TEMP/header.tmp
touch $TEMP/body.tmp
for i in ${SAMPLES[@]};
do  awk -F "," -v i="$i" '$2==i {print $2,"../../../../Isolates_assembly/Pool_"substr($2,1,3)"/02.Rm_adapters/fastq_clean/"$1".clean_1.fq.gz","../../../../Isolates_assembly/Pool_"substr($2,1,3)"/02.Rm_adapters/fastq_clean/"$1".clean_2.fq.gz"}' $TEMP/names.csv >> $TEMP/body.tmp;
done
cut -d "_" -f2- $TEMP/body.tmp | sort | awk 'FNR==1{
              print
              next
            }
            {
              A[$1]=$1 in A ? A[$1]","$2:$2
            }
            {
              B[$1]=$1 in B ? B[$1]","$3:$3
            }
         END{
              for(i in A)
                   print i,A[i],B[i]
            }' | tr " " "\t" >> $TEMP/body_merged.tmp # This strange awk command that I don't fully understand merges together rows by sample name and merges together the location of the r1 and r2 with commas.

cat $TEMP/header.tmp $TEMP/body_merged.tmp > $WORKDIR/samples-txt

rm -r $TEMP