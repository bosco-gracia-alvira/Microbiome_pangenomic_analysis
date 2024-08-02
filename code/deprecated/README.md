# Microbiome_pangenomic_analysis - code

This is the collection of scripts that I am writing to analyse each bacterial lineage. In principle each analysis can be ran independently, but some are composed of more than one script. All of the scripts start by asking you which species do you want to analyse.

## Anvi'o population genomics analysis

This workflow maps the reads of all the samples from the same lineage to a reference assembly that can be chosen. The scripts have to be run in the following order:

1. `gen_samples_txt.sh`
2. `gen_fasta_txt.sh`
3. `Anvio_workflow.sh`
4. `Anvio_display.sh`

## Average Nucleotide Identity calculation

The script `ANI.sh` calculates the pairwise ANI of all the samples from the same lineage and exports a .png file with the result.

## Evaluate for intra-specific microbial populations

The concept of species in bacteria is too broad. Individuals from the same species can belong to different populations, which poses a barrier to recombination. With `PopCOGenT.sh` script we can test if the bacteria assignated to the same species belong also to the same population.

## Phylogeny construction

The script `Phylogeny.sh` calls the genes of all the genomes from the desired lineage and builds the pangenome. Then, it extracts the SNPs from the core genes (those present in all the genomes), and concatenates them into a single alignment that is used to build the phylogeny with iqtree.

It lets you choose if you want to add the NCBI representative genome as outgroup to your alignment or not. On the one hand, with an outgroup you can root the tree, but at the same time less genes will be considered "core", because the SNPs in the genes that are not present in the outgroup will not be considered.
