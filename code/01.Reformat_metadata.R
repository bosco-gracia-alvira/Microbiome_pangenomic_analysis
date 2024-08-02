# This script merges genome taxonomy and sample metadata tables. It is run inside of 01.gen_fasta_txt.sh

library(data.table)
setwd("/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data")

# Load the data
sample_metadata <- fread("metadata.tsv")
genome_info <- fread("taxonomy.tsv")

# Extract the sample name from the genome name
genome_info$sample <- sub("_.*", "", genome_info$user_genome)
genome_info <- genome_info[,c(1,2,21)]

# Merge the data frames
merged_data <- merge(sample_metadata, genome_info, by.x = "sample", by.y = "sample")
merged_data <- merged_data[, c("user_genome","pool","sample","Temperature","Generation","Population","Replicate","classification")]

# Save the genome metadata table in the directory
write.table(merged_data, "Genome_metadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


#Now I will create three files needed for Anvi'o

# Anvio hates hyphens, so I replace them with underscores
merged_data$name <- gsub("-", "_", merged_data$user_genome)

# fasta-txt
# File fasta-txt contains two columns: name and path
fasta_txt <- merged_data[, c("name", "user_genome")]
fasta_txt$path <- paste0("../../Isolates/",fasta_txt$user_genome,".fasta")
fasta_txt <- fasta_txt[, c("name", "path")]

# Save fasta-txt table in the directory
write.table(fasta_txt, "fasta-txt", sep = "\t", row.names = FALSE, quote = FALSE)


# samples-txt
# File samples-txt contains three columns: sample, r1 and r2
samples_txt <- merged_data[, c("name", "pool", "sample")]
samples_txt$r1 <- paste0("../../../../Isolates_assembly/Pool_",samples_txt$pool,"/02.Rm_adapters/fastq_clean/",samples_txt$sample,".clean_1.fq.gz")
samples_txt$r2 <- paste0("../../../../Isolates_assembly/Pool_",samples_txt$pool,"/02.Rm_adapters/fastq_clean/",samples_txt$sample,".clean_2.fq.gz")
samples_txt$sample <- samples_txt$name
samples_txt <- samples_txt[, c("sample", "r1", "r2")]

# Save sample-txt table in the directory
write.table(samples_txt, "samples-txt", sep = "\t", row.names = FALSE, quote = FALSE)


# Anvio_misc
# This file contains information about the temperature regime and generation of each isolate
Anvio_misc <- merged_data[, c("name", "Temperature", "Generation")]

# Save the table in the directory
write.table(Anvio_misc, "Anvio_misc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
