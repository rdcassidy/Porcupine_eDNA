#Sorting out the data from the output of Moncho's pipeline 

library(dplyr)
library(tidyverse)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
library(Biostrings)

#1) isolating the metadata (i.e., library info etc) for only the Porcupine Seabight sites 

metadata <- read.csv("metadata.csv")

Porcupine_metadata <- metadata %>% filter(Project %in% c("Porcupine21" ,"Porcupine23" , "Porcupine22"))

#2) filtering the ASV table to only include porcupine sites

ASVs <- read.csv("ASV_table.csv")

ASVs_filtered <- ASVs %>% filter(Sample %in% Porcupine_metadata$sample_id)


write.csv(ASVs_filtered, "filtered_ASV_table.csv")
#3) filtering the hash key so that it only includes the 
Hash_key <- read.csv("Hash_key.csv")

filtered_Hash_key <- Hash_key %>% filter(Hash %in% ASVs_filtered$Hash)

write.csv(filtered_Hash_key, "filtered_hash_key.csv", row.names = F)
#4) write the filtered Hash key as a fasta file


# Assuming 'sequence' is the column with the actual sequences
seq_set <- DNAStringSet(filtered_Hash_key$sequence)

# Set the names of the sequences (sequence identifiers) from your CSV file
names(seq_set) <- filtered_Hash_key$Hash
# Write to FASTA file
writeXStringSet(seq_set, filepath = "Filtered_Hash_key.fasta", format = "fasta")

