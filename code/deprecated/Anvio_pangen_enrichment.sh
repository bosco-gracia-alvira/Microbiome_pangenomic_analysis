#!/bin/bash
# This script is used to estimate the metabolism of A. indonesiensis pangenome using Anvi'o. I am following the instructions of the Anvi'o tutorial:
# https://merenlab.org/tutorials/infant-gut/#chapter-v-metabolism-prediction
# Bosco Gracia Alvira, 2023

### VARIABLES
START="/Users/bgracia/Dropbox (PopGen)/Bosco/PhD_Dropbox"

### COMMANDS
eval "$(conda shell.bash hook)"
conda activate anvio-7.1

cd "$START" || exit

echo 'To which taxon do you want to generate a misc-table? Type in the terminal the species with the format "Genus_species" or the genus with the format "Genus"'
echo
cat "$START"/Isolates_assembly/Pool_???/07.GTDB-Tk/summary.tsv | cut -f2 | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -k2 | column
read

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

PANGEN="$START"/Microbiome_pangenomic_analysis/data/$REPLY/Anvio_pangen

if [[ ! -d "$PANGEN"/03_PAN ]]
then    
        echo
        echo
        echo -e "The species $REPLY is not availabe :("
        echo -e "Maybe you have forgotten to run Anvio_pangen_wf.sh"
        echo
        exit
else
        mkdir -p "$PANGEN"/05_ENRICHMENT
fi

# This chunk of code only works if the user is working with a whole genus.
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

        for i in {COG20_FUNCTION,COG20_CATEGORY,KOfam,KEGG_Module,COG20_PATHWAY,KEGG_Class}
        do      
                anvi-compute-functional-enrichment-across-genomes \
                                -e "$PANGEN"/my-external-genomes${bins}.txt \
                                -o "$PANGEN"/05_ENRICHMENT/Taxon_enrichment"${i}".txt \
                                -G "$PANGEN"/Enrichment.tmp \
                                --annotation-source "${i}"
        done

        rm "$PANGEN"/*.tmp

fi

if [[ -f "$PANGEN"/Groups.txt ]]
then
        for i in {COG20_FUNCTION,COG20_CATEGORY,KOfam,KEGG_Module,COG20_PATHWAY,KEGG_Class}
        do      
                anvi-compute-functional-enrichment-across-genomes \
                                                -e "$PANGEN"/my-external-genomes${bins}.txt \
                                                -o "$PANGEN"/05_ENRICHMENT/Enrichment"${i}".txt \
                                                -G "$PANGEN"/Groups.txt \
                                                --annotation-source "${i}"
        done
fi

# for i in $(cut -f6 "$WORKDIR"/05_ENRICHMENT/EnrichmentCOG20_FUNCTION.txt)
# do grep "$i" "$WORKDIR"/COG20.tsv
# done | sort > "$WORKDIR"/COG20-2.tsv

# for i in $(cut -f6 "$WORKDIR"/05_ENRICHMENT/EnrichmentCOG20_FUNCTION.txt)
# do echo "${i}";
# done | sort

# awk -F, 'FNR==NR {lines[$2]; next} $2 in lines' file1 file2

#anvi-compute-functional-enrichment-in-pan -p $WORKDIR/03_PAN/MYPAN-PAN.db \
                                          #-g $WORKDIR/03_PAN/MYPAN-GENOMES.db \
                                          #--category-variable clade \
                                          #--annotation-source $WORKDIR/COG20_FUNCTION \
                                          #-o $WORKDIR/05_ENRICHMENT/functional-enrichment.txt