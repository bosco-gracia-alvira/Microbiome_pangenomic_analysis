#!/bin/bash
# This script shows the taxa that are available for the pangenomic analysis.
# Bosco Gracia Alvira, 2024

### VARIABLES

# Set the path
WORKDIR="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data"
CODE="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/code"
ASSEMBLY="/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Isolates_assembly"

### COMMANDS
IFS="
"

# Create or update the taxonomy file
cat "$ASSEMBLY"/Pool_503/07.GTDB-Tk/gtdbtk.bac120.summary.tsv > "$WORKDIR"/taxonomy.tsv
tail -n +2 "$ASSEMBLY"/Pool_591/07.GTDB-Tk/gtdbtk.bac120.summary.tsv >> "$WORKDIR"/taxonomy.tsv
tail -n +2 "$ASSEMBLY"/Pool_643/07.GTDB-Tk/gtdbtk.bac120.summary.tsv >> "$WORKDIR"/taxonomy.tsv

# Show the available scripts
echo
echo "These are the scripts available for the pangenomic analysis"
echo "Run them in numerical order. The first argument is always the species you want to analyse"
echo
basename -a "$CODE"/0*.sh | sort | column
echo

# Show the available species
echo
echo "These are the species that are available for the pangenomic analysis"
echo "Copy the name of the species you want to analyse and paste it as first argument of the scripts"
echo
cut -f2 "$WORKDIR"/taxonomy.tsv | rev | cut -d "_" -f1 | rev | grep " "| sed 's/ /_/' | sort | uniq -c | sort -r | column
echo
