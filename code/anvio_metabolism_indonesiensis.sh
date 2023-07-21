#!/bin/bash
# This script is used to estimate the metabolism of A. indonesiensis pangenome using Anvi'o. I am following the instructions of the Anvi'o tutorial:
# https://merenlab.org/tutorials/infant-gut/#chapter-v-metabolism-prediction
# Bosco Gracia Alvira, 2023

### VARIABLES
WORKDIR="/Users/bgracia/Dropbox (PopGen)/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/Acetobacter_indonesiensis/Anvio_pangen"

### COMMANDS
eval "$(conda shell.bash hook)"
conda activate anvio-7.1

cd "$WORKDIR" || exit

if [[ ! -d 04_METABOLISM ]]
then    
        mkdir 04_METABOLISM
fi

# This command runs the annotated genomes against the KEGG database and estimates the completeness of each KEGG module (i.e. pathway) in each genome. It outputs a matrix with the completeness of each module in each genome, another matrix with the KOs (each component of the module) in each genome, and a third matrix with the modules that are complete in at least one genome.
anvi-estimate-metabolism -e "$WORKDIR"/my-external-genomes_no_bins.txt \
                         --matrix-format \
                         -O "$WORKDIR"/04_METABOLISM/A_indonesiensis \
                         --module-completion-threshold 1 
                         #--module-specific-matrices M00569,M00641,M00555,M00542,M00568,M00745 This option only computes the 6 modules that are specifically enriched in one clade.

# This command transforms the matrix into a newick tree, which is used to visualize the results in the interactive interface.
anvi-matrix-to-newick "$WORKDIR"/04_METABOLISM/A_indonesiensis-completeness-MATRIX.txt

# We create the profile database, which stores the information of the completeness matrix.
anvi-interactive -d "$WORKDIR"/04_METABOLISM/A_indonesiensis-completeness-MATRIX.txt \
                 -p "$WORKDIR"/04_METABOLISM/A_indonesiensis_metabolism_PROFILE.db \
                 --manual-mode \
                 --dry-run

# The following commands make the interactive interface prettier, including the names and pathways that the modules belong to.

# We learn where the KEGG/MODULES.db is:
export ANVIO_MODULES_DB=`python -c "import anvio; import os; print(os.path.join(os.path.dirname(anvio.__file__), 'data/misc/KEGG/MODULES.db'))"`

# Start an empty file:
echo -e "module\tclass\tcategory\tsubcategory\tname" > "$WORKDIR"/04_METABOLISM/modules_info.txt

# Get module classes:
sqlite3 $ANVIO_MODULES_DB "select module, data_value from kegg_modules where data_name='CLASS'" | \
    sed 's/; /|/g' | \
    tr '|' '\t' >> "$WORKDIR"/04_METABOLISM/module_class.txt

# Get module names:
sqlite3 $ANVIO_MODULES_DB "select module, data_value from kegg_modules where data_name='NAME'" | \
    tr '|' '\t' > "$WORKDIR"/04_METABOLISM/module_names.txt

# Join everything
paste "$WORKDIR"/04_METABOLISM/module_class.txt <(cut -f 2 "$WORKDIR"/04_METABOLISM/module_names.txt ) >> "$WORKDIR"/04_METABOLISM/modules_info.txt

# Empty the trash bin:
rm "$WORKDIR"/04_METABOLISM/module_names.txt "$WORKDIR"/04_METABOLISM/module_class.txt

# We import the KEGG information into the profile database:
anvi-import-misc-data "$WORKDIR"/04_METABOLISM/modules_info.txt \
                      -p "$WORKDIR"/04_METABOLISM/A_indonesiensis_metabolism_PROFILE.db \
                      -t items

# We import the information about the temperature and the clade of each genome. I have created the table Layers_no_bins.txt manually.
anvi-import-misc-data "$WORKDIR"/Layers_no_bins.txt \
                      -p "$WORKDIR"/04_METABOLISM/A_indonesiensis_metabolism_PROFILE.db \
                      -t layers \
                      --just-do-it

# We estimate the metabolism again, but this time the output includes the presence of each module and each KO in each genome. The utility of these files is learning more about a specific module that is differentially complete between the two clades.
anvi-estimate-metabolism -e "$WORKDIR"/my-external-genomes_no_bins.txt \
                         -O "$WORKDIR"/04_METABOLISM/A_indonesiensis_metabolism \
                         --module-completion-threshold 1 \
                         --kegg-output-modes modules,kofam_hits

# Finally we run the interactive interface, this time adding the newick tree.
anvi-interactive --manual-mode \
                 -d "$WORKDIR"/04_METABOLISM/A_indonesiensis-completeness-MATRIX.txt \
                 -t "$WORKDIR"/04_METABOLISM/A_indonesiensis-completeness-MATRIX.txt.newick \
                 -p "$WORKDIR"/04_METABOLISM/A_indonesiensis_metabolism_PROFILE.db \
                 --title "A. indonesiensis Metabolism Heatmap"

# Finally, we can run the enrichment analysis, that screens for modules that are significantly enriched in one clade. What it asks is: which proportion of the clades has this complete pathway? So if all the genomes of the clade have 0.8 completeness and the members of the other clade have 0, this won't be detected because the module is not complete. However, if ALL the members of a clade have a completeness of 1 and the genomes from the other clade have 0.8, then it will detected as enriched.
anvi-compute-metabolic-enrichment -M "$WORKDIR"/04_METABOLISM/A_indonesiensis_metabolism_modules.txt \
                                  -G "$WORKDIR"/clades.txt \
                                  --module-completion-threshold 1 \
                                  -o "$WORKDIR"/04_METABOLISM/A_indonesiensis_enriched_modules.txt
