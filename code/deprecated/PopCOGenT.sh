#!/bin/bash

cd ~/PhD

mkdir Microbiome_pangenomic_analysis/data/temp
TEMP=Microbiome_pangenomic_analysis/data/temp
cat Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv > $TEMP/taxonomy.tsv

echo 'Which of the following species do you want to analyse? Type it in the terminal with the format: "Genus_species"'
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
SAMPLES=($(awk -v s="$SPECIES" -F "\t" '$2 ~ s {print $1}' $TEMP/taxonomy.tsv))

for i in ${SAMPLES[@]};
do scp Isolates_assembly/Pool_$(echo $i | cut -f1 -d "_")/07.GTDB-Tk/Genomes/$i.fa vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/PopCOGenT/src/PopCOGenT/Genomes/;
done

ssh -T vetlinux05@pgnsrv043.vu-wien.ac.at << FOO

REPLY=$REPLY
cd ~/Bosco/PopCOGenT/src/PopCOGenT/

export BASH_ENV=~/.bashrc
export MUGSY_INSTALL=~/.local/miniconda3/envs/PopCOGenT/bin

source config.sh

eval "\$(conda shell.bash hook)"
conda activate PopCOGenT
source \${mugsy_env}

python get_alignment_and_length_bias.py --genome_dir \${genome_dir} --genome_ext \${genome_ext} --alignment_dir \${alignment_dir} --mugsy_path \${mugsy_path} --mugsy_env \${mugsy_env} --base_name \${base_name} --final_output_dir \${final_output_dir} --num_threads \${num_threads} \${keep_alignments}

python cluster.py --base_name \${base_name} --length_bias_file \${final_output_dir}/\${base_name}.length_bias.txt --output_directory \${final_output_dir} --infomap_path \${infomap_path} \${single_cell}

rm Genomes/*
rm -r proc/ __pycache__/ *log infomap_out
FOO

if [[ ! -d $WORKDIR/PopCOGenT ]]
then
    mkdir $WORKDIR/PopCOGenT
fi

rsync -avz --remove-source-files -e ssh vetlinux05@pgnsrv043.vu-wien.ac.at:/home/vetlinux05/Bosco/PopCOGenT/src/PopCOGenT/output/$REPLY"*" $WORKDIR/PopCOGenT/
mv $WORKDIR/PopCOGenT/${REPLY}_*.txt.cluster.tab.txt $WORKDIR/PopCOGenT/${REPLY}.cluster.tsv
mv $WORKDIR/PopCOGenT/${REPLY}_*.txt.unclust.graphml $WORKDIR/PopCOGenT/${REPLY}.unclust.graphml
rm -r $TEMP