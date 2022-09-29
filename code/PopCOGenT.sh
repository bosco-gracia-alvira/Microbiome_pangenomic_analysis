#!/bin/bash

cd ~/PhD

echo 'Which bacteria do you want to analyse? Type it in the terminal with the format: "Genus_species"'
read

SPECIES=$(echo $REPLY | sed 's/_/ /')

mkdir Microbiome_pangenomic_analysis/data/$REPLY/temp
WORKDIR=~/PhD/Microbiome_pangenomic_analysis/data/$REPLY
TEMP=Microbiome_pangenomic_analysis/data/$REPLY/temp

cat Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv > $TEMP/taxonomy.tsv

SAMPLES=($(awk -v s="$SPECIES" -F "\t" '$2 ~ s {print $1}' $TEMP/taxonomy.tsv))

for i in ${SAMPLES[@]};
do scp Isolates_assembly/Pool_$(echo $i | cut -f1 -d "_")/07.GTDB-Tk/Genomes/$i.fa vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/PopCOGenT/src/PopCOGenT/Genomes/;
done

ssh vetlinux05@pgnsrv043.vu-wien.ac.at << FOO
#!/bin/bash
REPLY=$REPLY
cd ~/Bosco/PopCOGenT/src/PopCOGenT/

export BASH_ENV=~/.bashrc
export MUGSY_INSTALL=~/.local/miniconda3/envs/PopCOGenT/bin

source config.sh

eval "$(conda shell.bash hook)"
conda activate PopCOGenT
source ${mugsy_env}

python get_alignment_and_length_bias.py --genome_dir ${genome_dir} --genome_ext ${genome_ext} --alignment_dir ${alignment_dir} --mugsy_path ${mugsy_path} --mugsy_env ${mugsy_env} --base_name ${base_name} --final_output_dir ${final_output_dir} --num_threads ${num_threads} ${keep_alignments}

python cluster.py --base_name ${base_name} --length_bias_file ${final_output_dir}/${base_name}.length_bias.txt --output_directory ${final_output_dir} --infomap_path ${infomap_path} ${single_cell}

#rm Genomes/*
rm -r proc/ __pycache__/ *log infomap_out
FOO

rsync -avz --remove-source-files -e ssh vetlinux05@pgnsrv043.vu-wien.ac.at:/home/vetlinux05/Bosco/PopCOGenT/src/PopCOGenT/output/$REPLY"*" $WORKDIR
rm -r $TEMP
#scp -r vetlinux05@pgnsrv043.vu-wien.ac.at:/home/vetlinux05/Bosco/PopCOGenT/src/PopCOGenT/output/* $WORKDIR