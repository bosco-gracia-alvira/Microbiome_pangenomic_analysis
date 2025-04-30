#!/bin/bash
# This script runs Anvio pangenomic analysis for the desired species
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the paths
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data"

### COMMANDS
IFS="
"

# First argument of the script is the species to analyse
REPLY=$1

if [[ -z "$REPLY" ]]
then
    echo "You need to provide the species you want to analyse as first argument"
    echo
    echo "These are the species that are available for the pangenomic analysis"
    echo "Copy the name of the species you want to analyse and paste it as first argument of the script"
    echo
    cut -f2 "$WORKDIR"/taxonomy.tsv | awk -F 's__' '{print $2}' | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
    echo
    exit
elif [[ "$REPLY" == "h" ]]
then
    echo "These are the species that are available for the pangenomic analysis"
    echo "Copy the name of the species you want to analyse and paste it as first argument of the script"
    echo
    cut -f2 "$WORKDIR"/taxonomy.tsv | awk -F 's__' '{print $2}' | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
    echo
    exit
fi

# If the fasta-txt is not available we cannot do the analysis!
if [[ ! -f "$WORKDIR"/$REPLY/Anvio_pangen/fasta-txt ]]
then    
        echo
        echo
        echo -e "The species $REPLY is not availabe :("
        echo -e "Maybe you have forgotten to run 01.Gen_tables.sh"
        echo
        exit
fi

eval "$(conda shell.bash hook)"
conda activate anvio-7.1

cd "$WORKDIR"/"$REPLY"/Anvio_pangen || exit

anvi-run-workflow -w contigs \
                  -c ../../../code/02.config-contig.json

anvi-run-workflow -w pangenomics \
                  -c ../../../code/02.config-pangen.json

# Import the misc data that was created in 01.Gen_tables.sh
anvi-import-misc-data \
        "$WORKDIR"/$REPLY/Anvio_pangen/Anvio_misc.tsv \
        -p "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/MYPAN-PAN.db \
        --just-do-it \
        -t layers

# Phylogenomic tree of the species
# Extract the SCGs from the pangenome. SCGs must be found in all the genomes and only once
SAMPLES=$(cut -f1 "$WORKDIR"/$REPLY/Anvio_pangen/fasta-txt | grep -v "name")
NUM=$(echo $SAMPLES | wc -w)

anvi-get-sequences-for-gene-clusters \
        -g "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/MYPAN-GENOMES.db \
        -p "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/MYPAN-PAN.db \
        -o "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/SCGs.fa \
        --min-num-genes-from-each-genome 1 \
        --max-num-genes-from-each-genome 1 \
        --max-functional-homogeneity-index 0.99 \
        --min-num-genomes-gene-cluster-occurs "$NUM" \
        --concatenate \
        > "$WORKDIR"/"$REPLY"/Anvio_pangen/00_LOGS/MYPAN-anvi_phylogeny.log 2>&1

# Generate a phylogenomic tree
anvi-gen-phylogenomic-tree \
        -f "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/SCGs.fa \
        -o "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/phylogenomic-tree.txt \
        >> "$WORKDIR"/"$REPLY"/Anvio_pangen/00_LOGS/MYPAN-anvi_phylogeny.log 2>&1

# We cannot include the tree as it is in the Anvi'o profile database, so we transform it into a misc-data file:
echo -e "item_name\tdata_type\tdata_value" > "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/header.tmp
echo -e "phylogeny\tnewick\t$(cat "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/phylogenomic-tree.txt)" > "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/body.tmp
cat "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/header.tmp "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/body.tmp > "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/phylogeny_order.txt
rm "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/*.tmp

# Finally we import it into the profile databases:
anvi-import-misc-data \
        "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/phylogeny_order.txt \
        -p "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/MYPAN-PAN.db \
        -t layer_orders \
        --just-do-it \
        >> "$WORKDIR"/"$REPLY"/Anvio_pangen/00_LOGS/MYPAN-anvi_phylogeny.log 2>&1

# Test if there are functional differences between the genomes from each temperature regime
mkdir -p "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/Functional_enrichment

if [[ "$REPLY" != "Lactiplantibacillus_plantarum" ]]
then
    category="Temperature"
else
    category="Genotype"
fi

for i in {COG20_FUNCTION,Pfam,KEGG_Module,KOfam,COG20_PATHWAY,KEGG_Class,COG20_CATEGORY}
do      
        anvi-compute-functional-enrichment-in-pan \
        -p "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/MYPAN-PAN.db \
        -g "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/MYPAN-GENOMES.db \
        --category-variable $category \
        --annotation-source "${i}" \
        -o "$WORKDIR"/"$REPLY"/Anvio_pangen/03_PAN/Functional_enrichment/"${i}"-enrichment.txt \
        >> "$WORKDIR"/"$REPLY"/Anvio_pangen/00_LOGS/MYPAN-anvi_enrichment.log 2>&1
done
