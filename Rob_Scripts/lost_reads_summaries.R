library(dplyr)
library(tidyverse)

setwd("~/Dropbox/Porcupain_project/Run_04_03_2025/")
summarytable <- read.csv("summary.csv")


summarydada <- read.csv(file = "summary_dada2.csv")



####going to check that the we are going to check how many read are lost after removing primers 
  ##by subtracting the fwd + rev reads from "all loci" for each barcdode/PCR combo 
  
results_eval <- summarytable %>%
  filter(locus != "Error") 

results_eval$nReads <- as.numeric(results_eval$nReads)
  
results_eval <- results_eval %>% group_by(fastq_header) %>%
  summarise(
    all_loci_reads = nReads[locus == "all_loci"],
    primer_sum = sum(nReads[locus %in% c("Locus_COI_Fwd", "Locus_COI_Rev")]),
    reads_lost_primer_filtering = all_loci_reads - primer_sum,
    .groups = "drop"
  ) %>% 
  mutate(perecent_reads_lost = reads_lost_primer_filtering/all_loci_reads)

####as we can see, some groups had really low numbers of reads in total. of the barcodes with a normal-seeming number of reads, the highest % lost reads from removing primers was like 8.5% 


###now going to check to see if anything important was lost in the dada2 step

loss_stats_summarydada <- summarydada %>% 
  mutate(pct_lost_filtering_Fwd = (inputF - filtered.F)/inputF) %>% 
  mutate(pct_lost_filtering_Rev = (inputR - filtered.R)/inputR) %>% 
  mutate(dif_dadaF1s_dadaF2s = abs(dadaF1s - dadaF2s)) %>% 
  mutate(dif_dadaR1s_dadaR2s = abs(dadaR1s - dadaR2s)) %>% 
  mutate(Pct_fwd_lost_dereplicating_merging = ((((dadaF1s + dadaF2s) / 2) - mergersF) / ((dadaF1s + dadaF2s) / 2))) %>% 
  mutate(Pct_rev_lost_dereplicating_merging = ((((dadaR1s + dadaR2s) / 2) - mergersR) / ((dadaR1s + dadaR2s) / 2))) %>% 
  mutate(pct_lost_chimeras = (joined - nochim) / joined) 



write.csv(loss_stats_summarydada, "lost_reads_table_run_04_03_2025.csv")


metadata_libraries <- read_csv("metadata.csv")
metadata_libraries <- metadata_libraries %>% filter(grepl("Porcupine", Project)) %>% select(sample_id, SampleType)

Porcupine_reads_lost <- left_join(loss_stats_summarydada ,metadata_libraries, by= c("Sample"= "sample_id")) 

Porcupine_reads_lost <- Porcupine_reads_lost %>% filter(Sample %in% metadata_libraries$sample_id)


Porcupine_reads_lost %>% 
  group_by(SampleType) %>% 
  ggplot(aes(x = SampleType, y =Pct_fwd_lost_dereplicating_merging)) +
  geom_violin() -> lostreadsplot

Porcupine_reads_lost <- Porcupine_reads_lost %>% mutate("Site" = Sample)
Porcupine_reads_lost <- Porcupine_reads_lost %>% separate(Site, into = c("Site", "trash"), sep = "_")                                          

Porcupine_reads_lost <- left_join(Porcupine_reads_lost, metadata2, by = "Site")                                        
#Porcupine_reads_lost <- Porcupine_reads_lost %>% select(-c(trash, Site))
Porcupine_reads_lost$LANCE[is.na(Porcupine_reads_lost$LANCE)] <- Porcupine_reads_lost$Site[is.na(Porcupine_reads_lost$LANCE)]


Porcupine_reads_lost %>% filter(SampleType =="sponge") %>% 
  group_by(LANCE) %>% 
  ggplot(aes(x = LANCE, y =Pct_fwd_lost_dereplicating_merging, fill = SampleType)) +
  geom_violin() +geom_point() -> lostreadsbylance_sponge

Porcupine_reads_lost %>% filter(SampleType =="eDNA") %>% 
  group_by(LANCE) %>% 
  ggplot(aes(x = LANCE, y =Pct_fwd_lost_dereplicating_merging, fill = SampleType)) +
  geom_violin() +geom_point() -> lostreadsbylance_eDNA






######with the original run ######

setwd("~/Dropbox/Porcupain_project/Run_04_03_2025/")
summarytable <- read.csv("summary.csv")


summarydada <- read.csv(file = "summary_dada2.csv")



####going to check that the we are going to check how many read are lost after removing primers 
##by subtracting the fwd + rev reads from "all loci" for each barcdode/PCR combo 

results_eval <- summarytable %>%
  filter(locus != "Error") 

results_eval$nReads <- as.numeric(results_eval$nReads)

results_eval <- results_eval %>% group_by(fastq_header) %>%
  summarise(
    all_loci_reads = nReads[locus == "all_loci"],
    primer_sum = sum(nReads[locus %in% c("Locus_COI_Fwd", "Locus_COI_Rev")]),
    reads_lost_primer_filtering = all_loci_reads - primer_sum,
    .groups = "drop"
  ) %>% 
  mutate(perecent_reads_lost = reads_lost_primer_filtering/all_loci_reads)

####as we can see, some groups had really low numbers of reads in total. of the barcodes with a normal-seeming number of reads, the highest % lost reads from removing primers was like 8.5% 


###now going to check to see if anything important was lost in the dada2 step

loss_stats_summarydada <- summarydada %>% 
  mutate(pct_lost_filtering_Fwd = (inputF - filtered.F)/inputF) %>% 
  mutate(pct_lost_filtering_Rev = (inputR - filtered.R)/inputR) %>% 
  mutate(dif_dadaF1s_dadaF2s = abs(dadaF1s - dadaF2s)) %>% 
  mutate(dif_dadaR1s_dadaR2s = abs(dadaR1s - dadaR2s)) %>% 
  mutate(Pct_fwd_lost_dereplicating_merging = ((((dadaF1s + dadaF2s) / 2) - mergersF) / ((dadaF1s + dadaF2s) / 2))) %>% 
  mutate(Pct_rev_lost_dereplicating_merging = ((((dadaR1s + dadaR2s) / 2) - mergersR) / ((dadaR1s + dadaR2s) / 2))) %>% 
  mutate(pct_lost_chimeras = (joined - nochim) / joined) 



write.csv(loss_stats_summarydada, "lost_reads_table_run_04_03_2025.csv")


metadata_libraries <- read_csv("metadata.csv")
metadata_libraries <- metadata_libraries %>% filter(grepl("Porcupine", Project)) %>% select(sample_id, SampleType)

Porcupine_reads_lost <- left_join(loss_stats_summarydada ,metadata_libraries, by= c("Sample"= "sample_id")) 

Porcupine_reads_lost <- Porcupine_reads_lost %>% filter(Sample %in% metadata_libraries$sample_id)


Porcupine_reads_lost %>% 
  group_by(SampleType) %>% 
  ggplot(aes(x = SampleType, y =Pct_fwd_lost_dereplicating_merging)) +
  geom_violin() -> lostreadsplot

Porcupine_reads_lost <- Porcupine_reads_lost %>% mutate("Site" = Sample)
Porcupine_reads_lost <- Porcupine_reads_lost %>% separate(Site, into = c("Site", "trash"), sep = "_")                                          

Porcupine_reads_lost <- left_join(Porcupine_reads_lost, metadata2, by = "Site")                                        
#Porcupine_reads_lost <- Porcupine_reads_lost %>% select(-c(trash, Site))
Porcupine_reads_lost$LANCE[is.na(Porcupine_reads_lost$LANCE)] <- Porcupine_reads_lost$Site[is.na(Porcupine_reads_lost$LANCE)]


Porcupine_reads_lost %>% filter(SampleType =="sponge") %>% 
  group_by(LANCE) %>% 
  ggplot(aes(x = LANCE, y =Pct_fwd_lost_dereplicating_merging, fill = SampleType)) +
  geom_violin() +geom_point() -> lostreadsbylance_sponge

Porcupine_reads_lost %>% filter(SampleType =="eDNA") %>% 
  group_by(LANCE) %>% 
  ggplot(aes(x = LANCE, y =Pct_fwd_lost_dereplicating_merging, fill = SampleType)) +
  geom_violin() +geom_point() -> lostreadsbylance_eDNA




#creating a multi plot like monchos 

library(ggplot2)
library(patchwork)

# Create individual plots

plot1 <- ggplot(Porcupine_reads_lost, aes(x = SampleType, y = pct_lost_filtering_Fwd)) +
  geom_violin() +geom_point()+ ylim(0,1) + ggtitle("Filtering_fwd")

plot2 <- ggplot(Porcupine_reads_lost, aes(x = SampleType, y = pct_lost_filtering_Rev)) +
  geom_violin() +geom_point()+ ylim(0,1) + ggtitle("Filtering_rev")

plot3 <- ggplot(Porcupine_reads_lost, aes(x = SampleType, y = Pct_fwd_lost_dereplicating_merging)) +
  geom_violin() +geom_point()+ ylim(0,1) + ggtitle("Lost merging and dereplicating forward")

plot4 <- ggplot(Porcupine_reads_lost, aes(x = SampleType, y = Pct_rev_lost_dereplicating_merging)) +
  geom_violin() +geom_point() + ylim(0,1) + ggtitle("Lost merging and dereplicating reverse ")

plot5 <- ggplot(Porcupine_reads_lost, aes(x = SampleType, y = pct_lost_chimeras)) +
  geom_violin() +geom_point() + ylim(0,1) + ggtitle("Chimeras")

# Combine plots using patchwork
(plot1 | plot2  ) / (plot3 | plot4 | plot5)   # Arranges two side by side, then one below

