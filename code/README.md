# Microbiome_pangenomic_analysis - code

This is the collection of scripts that I am writing to analyse each bacterial lineage. In principle each analysis can be ran independently, but it is advised to run all the scripts in the established order, even if they are from different analyses.

## Anvi'o pangenomics

This workflow compiles all the genomes from the desired taxon (species level, GTDB) and builds a pangenome. The scripts have to be run in the following order:

1. `01.Gen_tables.sh`: records which genomes from the assembly collection belong to the desired species and creates the tables that Anvi'o needs subsetting the assembly collection to the desired species. This script requires the presence of `01.Reformat_metadata.R` in the same folder.

2. `02.Anvio_pangen_wf.sh` runs the Anvi'o pangenomic workflow using as input the tables with the genome locations. It computes automatically the Average Nucleotide Identity (ANI) and imports the metadata (isolate's temperature regime, fly generation...). After that, the script textracts from the pangenome the single copy genes (SCGs) that are not identical across the replicates and uses them to build the phylogeny. Finally, it makes an enrichment analysis; it looks for functions (pfams, COGs, KEGGs...) that are significantly enriched in the genomes from each temperature regime. This script requires two extra config.json files: `02.config-contig.json` and `02.config-pangen.json`. The output can be found in "03_PAN/Functional_enrichment/". The pangenome web interface can be displayed with the script `02.Anvio_pangen_display.sh`. 


## Anvi'o population genomics

