library(here)
library(dplyr)
library(rentrez)
library(tidyverse)

BLAST_results_multi_stax <- 
  read_table(file = here("CESGA/Resultados_blast/resultados_blast_07_02_2025.txt"), col_names = c("qseqid", "sseqid",  "pident", "length" ,"mismatch", "gapopen" ,"qstart","qend", "sstart", "send" ,"evalue" ,"bitscore", "qlen", "staxid")) 


BLAST_results_multi_stax <- BLAST_results_multi_stax %>% filter(grepl(";", staxid) ) %>% filter(pident >= 90, length >= 0.7*qlen)


BLAST_results_multi_stax <- BLAST_results_multi_stax %>%
  group_by(qseqid) %>%
  slice_max(pident, n = 1) %>%
  ungroup()


#this extracts only the accession number from the sseqid column. i have to change this still so that it works for the refseq database as well as genbank

# blast_results_acc <- BLAST_results_multi_stax %>%
#   mutate(accession = str_extract(sseqid, "(?<=\\|gb\\|)[^.]+"))
# 
# 
# blast_results_acc <- BLAST_results_multi_stax %>%
#   mutate(accession = str_extract(sseqid, "(?<=\\|(?:gb|ref)\\|)[^|]+"))

blast_results_acc <- BLAST_results_multi_stax %>%
  mutate(accession = str_extract(sseqid, "(?<=\\|(gb|ref|emb|dbj)\\|)[^|]+")) 

blast_results_acc<- as.data.frame(blast_results_acc)

accessions <- blast_results_acc %>% filter(accession != "NA")
accessions <- accessions$accession
write_lines(accessions, "accessions.txt")


######trying to get taxids



# Load accession numbers

install.packages("rentrez")
library(rentrez)

accessions <- readLines("accessions.txt")

# Function to get Taxonomy ID for each accession
get_taxid <- function(acc) {
  search_result <- entrez_search(db = "nuccore", term = acc)
  if (length(search_result$ids) > 0) {
    taxid_result <- entrez_link(dbfrom = "nuccore", id = search_result$ids, db = "taxonomy")
    return(taxid_result$links$nuccore_taxonomy[1])
  } else {
    return(NA)  # Return NA if no match found
  }
}

# Apply function to all accessions
tax_ids <- sapply(accessions, get_taxid)

# Create a data frame
results <- data.frame(accession = accessions, TaxID = tax_ids)

# Save to a CSV file
write.csv(results, "accession_taxids.csv", row.names = FALSE)

print("Taxonomy ID retrieval complete! Check accession_taxids.csv")

#adding the good staxids back in

results_distinct <- results %>% arrange() %>% distinct()
results_light <- results %>% select(accession, TaxID)

results_combined <- left_join(blast_results_acc, results_light, by = "accession")

BLAST_results_fixed <- results_combined %>% select(-c(staxid, accession)) 
