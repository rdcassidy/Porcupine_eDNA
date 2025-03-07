#fixing staxid again 


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


# blast_results_acc <- BLAST_results_multi_stax %>%
#   mutate(accession = str_extract(sseqid, "(?<=\\|gb\\|)[^.]+"))

blast_results_acc <- BLAST_results_multi_stax %>%
  mutate(accession = str_extract(sseqid, "(?<=\\|(gb|ref|emb|dbj)\\|)[^.]+"))

results <- read.csv("accession_taxids.csv")

results <- results %>% rename(accession = Accession)
joined_taxids <- left_join(blast_results_acc, results, by ="accession", relationship = "many-to-many" )

joined_taxids <- joined_taxids %>% distinct()


joined_taxids <- joined_taxids %>%  select(-c(accession,staxid)) %>% rename(staxid = TaxID)




BLAST_results <- 
  read_table(file = here("CESGA/Resultados_blast/resultados_blast_07_02_2025.txt"), col_names = c("qseqid", "sseqid",  "pident", "length" ,"mismatch", "gapopen" ,"qstart","qend", "sstart", "send" ,"evalue" ,"bitscore", "qlen", "staxid")) 

BLAST_results <- BLAST_results %>%
  group_by(qseqid) %>%
  slice_max(pident, n = 1) %>%
  ungroup()

BLAST_resultsxx <- BLAST_results %>% filter(grepl(";", staxid) ) %>% filter(pident >= 90, length >= 0.7*qlen)




joined_taxids <- joined_taxids %>% select(sseqid, staxid) %>% distinct()

df_x_to_update <- BLAST_results %>%
  filter(grepl(";", staxid)) %>%
  filter(pident >= 90, length >= 0.7*qlen) %>% 
  left_join(joined_taxids, by = "sseqid") %>%
  mutate(staxid = ifelse(!is.na(staxid.y), staxid.y, staxid.x)) %>%
  select(-staxid.x, -staxid.y)  # Remove temporary columns

# Keep rows where staxid does NOT contain ";"
df_x_unchanged <- BLAST_results %>%
  filter(!grepl(";", staxid))
df_x_unchanged$staxid <- as.numeric(df_x_unchanged$staxid)

# Combine updated and unchanged rows back into one dataframe
df_x_final <- bind_rows(df_x_to_update, df_x_unchanged)


