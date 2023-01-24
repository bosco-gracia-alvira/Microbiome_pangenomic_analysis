#!/bin/bash

echo 'To which species do you want to generate a misc-table? Type in the terminal the species with the format: "Genus_species"'
echo
cat Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv | cut -f2 | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
read

if [[ ! -d ~/PhD/Microbiome_pangenomic_analysis/data/$REPLY/Anvio_pangen/03_PAN ]]
then    
        echo
        echo
        echo -e "The species $REPLY is not availabe :("
        echo -e "Maybe you have forgotten to run Anvio_pangen_wf.sh"
        echo
        exit
fi

cd ~/PhD/Microbiome_pangenomic_analysis/data/$REPLY/Anvio_pangen/

mkdir KEGG_metabolism && cd KEGG_metabolism

anvi-estimate-metabolism -e ../my-external-genomes.txt \
                         -O $REPLY \
                         --matrix-format

anvi-matrix-to-newick ${REPLY}-completeness-MATRIX.txt

# dry run to get the profile db:
anvi-interactive -d ${REPLY}-completeness-MATRIX.txt \
                 -p ${REPLY}_metabolism_PROFILE.db \
                 --manual-mode \
                 --dry-run

# import the state file:
anvi-import-state -s additional-files/state-files/state-metabolism.json \
                  -p ${REPLY}_metabolism_PROFILE.db \
                  -n default

# run for reals:
anvi-interactive --manual-mode \
                 -d ${REPLY}-completeness-MATRIX.txt \
                 -t ${REPLY}-completeness-MATRIX.txt.newick \
                 -p ${REPLY}_metabolism_PROFILE.db \
                 --title "Acetobacter_indonesiensis Metabolism Heatmap"

##I stop here
# learn where the MODULES.db is:
export ANVIO_MODULES_DB=`python3 -c "import anvio; import os; print(os.path.join(os.path.dirname(anvio.__file__), '/usr/local/Caskroom/miniconda/base/envs/anvio-7.1/lib/python3.6/site-packages/anvio/data/misc/KEGG/MODULES.db'))"`

# start an empty file:
echo -e "module\tclass\tcategory\tsubcategory\tname" > modules_info.txt

# get module classes:
sqlite3 $ANVIO_MODULES_DB "select module, data_value from modules where data_name='CLASS'" | \
    sed 's/; /|/g' | \
    tr '|' '\t' >> module_class.txt

# get module names:
sqlite3 $ANVIO_MODULES_DB "select module,data_value from modules where data_name='NAME'" | \
    tr '|' '\t' > module_names.txt

# join everything
paste module_class.txt <(cut -f 2 module_names.txt ) >> modules_info.txt

# empty the trash bin:
rm module_names.txt module_class.txt

anvi-import-misc-data additional-files/metabolism/modules_info.txt \
                      -p Enterococcus_metabolism_PROFILE.db \
                      -t items

anvi-import-state -s additional-files/metabolism/metabolism_state.json \
                  -p Enterococcus_metabolism_PROFILE.db \
                  -n default