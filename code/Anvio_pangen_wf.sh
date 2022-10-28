https://merenlab.org/2018/07/09/anvio-snakemake-workflows/#pangenomics-workflow


anvi-run-workflow -w pangenomics \
                  -c config-pangenomics.json \
                  --save-workflow-graph

{
    "workflow_name": "pangenomics",
    "config_version": "2",
    "project_name": "MYPAN",
    "external_genomes": "my-external-genomes.txt",
    "fasta_txt": "fasta.txt",
    "output_dirs": {
        "FASTA_DIR": "01_FASTA_contigs_workflow",
        "CONTIGS_DIR": "02_CONTIGS_contigs_workflow",
        "LOGS_DIR": "00_LOGS_pan_workflow"
    }
}