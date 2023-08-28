#!/bin/bash
# This script is used to estimate the metabolism of the taxa pangenomes using Anvi'o. I am following the instructions of the Anvi'o tutorial:
# https://merenlab.org/tutorials/infant-gut/#chapter-v-metabolism-prediction
# Bosco Gracia Alvira, 2023

### VARIABLES
WORKDIR="/Users/bgracia/Dropbox (PopGen)/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data"
START="/Users/bgracia/Dropbox (PopGen)/Bosco/PhD_Dropbox"

### COMMANDS

# First we activate the Anvi'o environment:
eval "$(conda shell.bash hook)"
conda activate anvio-7.1

# This command outputs the available taxa that we can analyse and lets you choose one:
echo 'To which taxon do you want to generate a misc-table? Type in the terminal the species with the format "Genus_species" or the genus with the format "Genus"'
echo
cat "$START"/Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv | cut -f2 | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -k2 | column
read -r REPLY

# This lets you choose if you want to use incomplete genomes (bins) or not:
echo 'Would you like to include incomplete genomes (bins) in the analysis? Type "yes" or "no"'
echo
read -r BINS

if [[ $BINS = "no" ]]
then
    bins=_no_bins
elif [[ $BINS = "yes" ]]
then
    bins=
else
    echo "You have not typed a valid answer. Please, type \"yes\" or \"no\""
    read -r BINS
fi

# If the taxon exists and Anvio_pangen_wf.sh had been run, the script continues. Otherwise, it stops.
if [[ ! -d "$WORKDIR/${REPLY}/Anvio_pangen/03_PAN" ]]
then    
        echo
        echo
        echo -e "The species $REPLY is not availabe :("
        echo -e "Maybe you have forgotten to run Anvio_pangen_wf.sh"
        echo
        #exit
    else
        mkdir -p "$WORKDIR/${REPLY}/Anvio_pangen/04_METABOLISM"
fi

METABOLISM="$WORKDIR/${REPLY}/Anvio_pangen/04_METABOLISM"
PANGEN="$WORKDIR/${REPLY}/Anvio_pangen"

# This command runs the annotated genomes against the KEGG database and estimates the completeness of each KEGG module (i.e. pathway) in each genome. It outputs a matrix with the completeness of each module in each genome, another matrix with the KOs (each component of the module) in each genome, and a third matrix with the modules that are complete in at least one genome.
anvi-estimate-metabolism -e "$PANGEN"/my-external-genomes${bins}.txt \
                         --matrix-format \
                         -O "$METABOLISM"/"${REPLY}"${bins} \
                         --module-completion-threshold 1

anvi-matrix-to-newick "$METABOLISM"/"${REPLY}"${bins}-completeness-MATRIX.txt
anvi-matrix-to-newick "$METABOLISM"/"${REPLY}"${bins}-ko_hits-MATRIX.txt

# Here we create the metabolism profile db:
anvi-interactive -d "$METABOLISM"/"${REPLY}"${bins}-completeness-MATRIX.txt \
                 -p "$METABOLISM"/"${REPLY}"${bins}_metabolism_PROFILE.db \
                 --manual-mode \
                 --dry-run

# Here we create the KO-hits profile db, because we cannot store them in the same as the modules:
anvi-interactive -d "$METABOLISM"/"${REPLY}"${bins}-ko_hits-MATRIX.txt \
                 -p "$METABOLISM"/"${REPLY}"${bins}_KO_PROFILE.db \
                 --manual-mode \
                 --dry-run

# The following chunk makes the interactive interface prettier, including the names and pathways that the modules belong to.

# We learn where the KEGG/MODULES.db is:
export ANVIO_MODULES_DB=`python -c "import anvio; import os; print(os.path.join(os.path.dirname(anvio.__file__), 'data/misc/KEGG/MODULES.db'))"`

# Start an empty file:
echo -e "module\tclass\tcategory\tsubcategory\tname" > "$METABOLISM"/modules_info.txt

# Get module classes:
sqlite3 $ANVIO_MODULES_DB "select module, data_value from kegg_modules where data_name='CLASS'" | \
    sed 's/; /|/g' | \
    tr '|' '\t' >> "$METABOLISM"/module_class.txt

# Get module names:
sqlite3 $ANVIO_MODULES_DB "select module, data_value from kegg_modules where data_name='NAME'" | \
    tr '|' '\t' > "$METABOLISM"/module_names.txt

# Join everything
paste "$METABOLISM"/module_class.txt <(cut -f 2 "$METABOLISM"/module_names.txt ) >> "$METABOLISM"/modules_info.txt

# Empty the trash bin:
rm "$METABOLISM"/module_names.txt "$METABOLISM"/module_class.txt

# We import the KEGG information into the profile databases:
anvi-import-misc-data "$METABOLISM"/modules_info.txt \
                      -p "$METABOLISM"/"${REPLY}"${bins}_metabolism_PROFILE.db \
                      -t items

# We do the same but with the KO database
export ANVIO_KOs_DB=`python -c "import anvio; import os; print(os.path.join(os.path.dirname(anvio.__file__), 'data/misc/KEGG/ko_list.txt'))"`

cut -f1,12 "$ANVIO_KOs_DB" | sed "s/^knum/KO/" > "$METABOLISM"/ko_list.txt

anvi-import-misc-data "$METABOLISM"/ko_list.txt \
                      -p "$METABOLISM"/"${REPLY}"${bins}_KO_PROFILE.db \
                      -t items

# In this part we create a phylogenetic tree based on the concatenated proteins that are present in all the genomes of the collection but just once, i.e., core single copy genes. This tree is used to visualize the results in the interactive interface.

# We get the number of genomes in the collection:
if [[ $BINS = "no" ]]
then
    NUM=$(grep "${REPLY}" "$START"/Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv | grep -v "bin" | wc -l | cut -d " " -f7)
else
    NUM=$(grep "${REPLY}" "$START"/Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv | wc -l | cut -d " " -f7)
fi

# We extract the core single copy genes:
anvi-get-sequences-for-gene-clusters -g "$PANGEN"/03_PAN/MYPAN-GENOMES.db \
                                     -p "$PANGEN"/03_PAN/MYPAN-PAN.db \
                                     --min-num-genomes-gene-cluster-occurs $(($NUM)) \
                                     --max-num-genes-from-each-genome 1 \
                                     --min-num-genes-from-each-genome 1 \
                                     --concatenate-gene-clusters \
                                     -o "$PANGEN"/03_PAN/concatenated-proteins.fa

# If we don't want to include bins, we remove their genes from the concatenated proteins:
if [[ $BINS = "no" ]]
then
    seqkit grep -v -r -p "bin_" "$PANGEN"/03_PAN/concatenated-proteins.fa -o "$PANGEN"/03_PAN/concatenated-proteins${bins}.fa
fi

# We build the tree using the FastTree algorithm:
anvi-gen-phylogenomic-tree -f "$PANGEN"/03_PAN/concatenated-proteins${bins}.fa \
                           -o "$PANGEN"/03_PAN/phylogeny${bins}.txt

# We cannot include the tree as it is in the Anvi'o profile database, so we transform it into a misc-data file:
echo -e "item_name\tdata_type\tdata_value" > "$PANGEN"/03_PAN/header.tmp
echo -e "phylogeny\tnewick\t$(cat "$PANGEN"/03_PAN/phylogeny${bins}.txt)" > "$PANGEN"/03_PAN/body.tmp
cat "$PANGEN"/03_PAN/header.tmp "$PANGEN"/03_PAN/body.tmp > "$PANGEN"/03_PAN/phylogeny_order${bins}.txt
rm "$PANGEN"/03_PAN/*.tmp

# Finally we import it into the profile databases:
anvi-import-misc-data "$PANGEN"/03_PAN/phylogeny_order${bins}.txt \
                      -p "$METABOLISM"/"${REPLY}"${bins}_metabolism_PROFILE.db \
                      -t layer_orders \
                      --just-do-it

anvi-import-misc-data "$PANGEN"/03_PAN/phylogeny_order${bins}.txt \
                      -p "$METABOLISM"/"${REPLY}"${bins}_KO_PROFILE.db \
                      -t layer_orders \
                      --just-do-it

# We estimate the metabolism again, but this time the output includes the presence of each module and each KO in each genome. The utility of these files is learning more about a specific module that is differentially complete between the two clades.
anvi-estimate-metabolism -e "$PANGEN"/my-external-genomes${bins}.txt \
                         -O "$METABOLISM"/"${REPLY}"${bins}_metabolism \
                         --module-completion-threshold 1 \
                         --kegg-output-modes modules,kofam_hits

# If the input does not contain "_" means that you have given the genus name instead of species_genus. This command associates the genome of the selection to the taxonomy at species level.
if [[ ! $REPLY =~ "_" ]]
then
        cat "$START"/Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv | grep "${REPLY}" | cut -f1  | cut -d "_" -f2- > "$PANGEN"/Genomes.tmp
        cat "$START"/Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv | grep "${REPLY}" | cut -f2 | rev | cut -d "_" -f1 | rev | grep " " | sed 's/ /_/' > "$PANGEN"/Clade.tmp

        paste "$PANGEN"/Genomes.tmp "$PANGEN"/Clade.tmp > "$PANGEN"/Bins_to_taxonomy.tmp
        echo -e "name\tspecies" > "$PANGEN"/header.tmp
        cat "$PANGEN"/header.tmp "$PANGEN"/Bins_to_taxonomy.tmp > "$PANGEN"/Bins_to_taxonomy.txt

    if [[ $BINS = "no" ]]
    then
        grep -v "bin_" "$PANGEN"/Bins_to_taxonomy.txt > "$PANGEN"/Bins_to_taxonomy${bins}.txt
    fi

        sed "s/species/group/" "$PANGEN"/Bins_to_taxonomy${bins}.txt > "$PANGEN"/Enrichment.tmp

anvi-import-misc-data "$PANGEN"/Bins_to_taxonomy${bins}.txt \
                      -p "$METABOLISM"/"${REPLY}"${bins}_metabolism_PROFILE.db \
                      -t layers \
                      --just-do-it

anvi-import-misc-data "$PANGEN"/Bins_to_taxonomy${bins}.txt \
                      -p "$METABOLISM"/"${REPLY}"${bins}_KO_PROFILE.db \
                      -t layers \
                      --just-do-it

anvi-compute-metabolic-enrichment \
                                -M "$METABOLISM"/"${REPLY}"${bins}_metabolism_modules.txt \
                                -G "$PANGEN"/Enrichment.tmp \
                                --module-completion-threshold 1 \
                                -o "$METABOLISM"/"${REPLY}"${bins}_enriched_modules.txt

rm "$PANGEN"/*.tmp

fi

# run for reals:
anvi-interactive --manual-mode \
                 -d "$METABOLISM"/${REPLY}${bins}-completeness-MATRIX.txt \
                 -t "$METABOLISM"/${REPLY}${bins}-completeness-MATRIX.txt.newick \
                 -p "$METABOLISM"/${REPLY}${bins}_metabolism_PROFILE.db \
                 --title "${REPLY} Metabolism Heatmap (Modules)"

anvi-interactive --manual-mode \
                 -d "$METABOLISM"/${REPLY}${bins}-ko_hits-MATRIX.txt \
                 -t "$METABOLISM"/${REPLY}${bins}-ko_hits-MATRIX.txt.newick \
                 -p "$METABOLISM"/${REPLY}${bins}_KO_PROFILE.db \
                 --title "${REPLY} Metabolism Heatmap (KOs)"
