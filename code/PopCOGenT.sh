#!/bin/bash



ACETOBACTERACEAE=($(awk -F'\t' 'BEGIN {ORS = " "} $8 ~ /Acetobacteraceae/ {print $1}' Genomes_info.txt))
LACTOBACILLACEAE=($(awk -F'\t' 'BEGIN {ORS = " "} $8 ~ /Lactobacillaceae/ {print $1}' Genomes_info.txt))
INDO=($(awk -F'\t' 'BEGIN {ORS = " "} $8 ~ /indonesiensis/ {print $1}' Genomes_info.txt))

for i in ${ACETOBACTERACEAE};
do scp 00.Genomes/${i}.fa vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Genomes/Acetobacteraceae ;
done

for i in ${LACTOBACILLACEAE};
do scp 00.Genomes/${i}.fa vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Genomes/Lactobacillaceae ;
done

for i in ${INDO};
do scp 00.Genomes/${i}.fa vetlinux05@pgnsrv043.vu-wien.ac.at:~/Bosco/Genomes/Indonesiensis ;
done

ssh vetlinux05@pgnsrv043.vu-wien.ac.at
cd ~/Bosco/PopCOGenT/src/PopCOGenT/

for i in $(ls *.fa);
do mv $i.renamed.mugsy $i;
done

rm -r proc/ __pycache__/ *log infomap_out


#export BASH_ENV=~/.bashrc
#export MUGSY_INSTALL=~/.local/miniconda3/envs/PopCOGenT/bin

source config_Aceto.sh
source config_Lacto.sh
source config_Indo.sh

eval "$(conda shell.bash hook)"
conda activate PopCOGenT
source ${mugsy_env}

python get_alignment_and_length_bias.py --genome_dir ${genome_dir} --genome_ext ${genome_ext} --alignment_dir ${alignment_dir} --mugsy_path ${mugsy_path} --mugsy_env ${mugsy_env} --base_name ${base_name} --final_output_dir ${final_output_dir} --num_threads ${num_threads} ${keep_alignments}

python cluster.py --base_name ${base_name} --length_bias_file ${final_output_dir}/${base_name}.length_bias.txt --output_directory ${final_output_dir} --infomap_path ${infomap_path} ${single_cell}

exit

scp -r vetlinux05@pgnsrv043.vu-wien.ac.at:/home/vetlinux05/Bosco/PopCOGenT/src/PopCOGenT/output 04.PopCOGenT/