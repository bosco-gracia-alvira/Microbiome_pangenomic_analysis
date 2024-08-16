# This script plots a PCA of the eigenvectors calculated in 05.SNPs_analysis.sh
library(lattice)
library(data.table)
library(Biostrings)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(stringr)
library(gtools)

# Get the command line arguments
args <- commandArgs(trailingOnly = TRUE)

# The first argument will be the value of $REPLY
reply <- args[1]

# Print the value of reply (for debugging purposes)
print(paste("The value of REPLY is:", reply))

setwd(paste0("/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/",reply,"/SNPs_analysis"))

# Paths to the files that we want to import
data_path <- paste0("/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/",reply,"/SNPs_analysis/")
metadata_isolates_path <- paste0("/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/data/",reply,"/Anvio/Anvio_misc.tsv")
visuals_path <- "/Users/bgracia/PopGen Dropbox/Martin McFly/Bosco/PhD_Dropbox/Microbiome_pangenomic_analysis/visuals/"

pca <- read.table(paste0(data_path,reply,".eigenvec"),sep=" ",header=F)
colnames(pca)[1] <- "name"
pca[,2] <- NULL

# How much variance is explained by each PC?
eigenval <- scan(paste0(data_path,reply,".eigenval"))
variance_explained <- eigenval / sum(eigenval)
percentage_variance_explained <- variance_explained * 100

# Import the metadata
metadata_isolates <- fread(metadata_isolates_path, header=T)
metadata_isolates$Source <- "Isolate"

metadata_pools <- as.data.frame(pca$name[grepl("^cF|^hF", pca$name)])
colnames(metadata_pools) <- "name"
metadata_pools <- metadata_pools %>%
  mutate(
    Temperature = case_when(
      str_detect(name, "F0") ~ "BASE",
      str_detect(name, "^h") ~ "HOT",
      str_detect(name, "^c") ~ "COLD"),
    Generation = str_extract(name, "F\\d+"),
    Source = "Pool")

metadata <- rbind(metadata_isolates,metadata_pools)

# Create a colour palette for our temperatures
temp_palette <- c("HOT" = "red", "BASE" = "grey", "COLD" = "lightblue")


# Merge the eigenvectors and the metadata
pca_merged <- merge(pca, metadata, by = "name")
pca_pool <- subset(pca_merged, pca_merged$Source=="Pool")
pca_isolate <- subset(pca_merged, pca_merged$Source=="Isolate")

# Plot the PCA
ggplot() +
  geom_point(data = pca_isolate, aes(x = V3, y = V4, color = Temperature, shape = Source), size = 1, alpha = 0.4) +
  geom_point(data = pca_pool, aes(x = V3, y = V4, color = Temperature, shape = Source), size = 2) +
  geom_text(data = pca_pool, aes(x = V3, y = V4, label = name, color = Temperature), size = 2, hjust = -0.5) +
  labs(
    x = paste0("PC1 (", round(percentage_variance_explained[1], 2), "% variance)"),
    y = paste0("PC2 (", round(percentage_variance_explained[2], 2), "% variance)"),
    color = "Temperature", shape = "Source") +
  scale_colour_manual(values = c("HOT" = "red", "BASE" = "grey", "COLD" = "lightblue", "COLD-CONSTANT" = "blue")) +
  theme_minimal()
ggsave(filename = paste0(visuals_path,reply,"_SNPs_PCA.png"), width = 6, height = 4, dpi = 300)
