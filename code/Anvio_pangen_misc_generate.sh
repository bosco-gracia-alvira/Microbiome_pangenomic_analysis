#!/bin/bash

cd ~/PhD

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

WORKDIR=~/PhD/Microbiome_pangenomic_analysis/data/$REPLY
echo -e 'samples\tclade\ttemperature' > $WORKDIR/Anvio_pangen/header.tmp
awk 'FNR>1 {print substr($1,5),"\""$2"\""}' $WORKDIR/PopCOGenT/$REPLY.cluster.tsv > $WORKDIR/Anvio_pangen/clade.tmp
awk -F'_' 'FNR>1 {if (int(substr($3,2)) <= 10) print "HOT"; else print "COLD"}' $WORKDIR/PopCOGenT/$REPLY.cluster.tsv > $WORKDIR/Anvio_pangen/temperature.tmp
paste $WORKDIR/Anvio_pangen/clade.tmp $WORKDIR/Anvio_pangen/temperature.tmp > $WORKDIR/Anvio_pangen/body.tmp
cat $WORKDIR/Anvio_pangen/header.tmp $WORKDIR/Anvio_pangen/body.tmp | sed 's/ /\t/g' > $WORKDIR/Anvio_pangen/Layers.txt
rm $WORKDIR/Anvio_pangen/*.tmp

anvi-import-misc-data   $WORKDIR/Anvio_pangen/Layers.txt \
                        -p $WORKDIR/Anvio_pangen/03_PAN/MYPAN-PAN.db \
                        --target-data-table layers \
                        --just-do-it