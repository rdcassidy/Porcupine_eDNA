#comparing the way the sponges filter 
sponge_sampler_metadata <- read_excel("P23_Muestras_EsponjasFiltros_eDNA_AnaRiesgo_FULLMETADATA.xlsx", sheet = "Esponjas" )

sponge_sampler_metadata <- sponge_sampler_metadata %>% select(ID,`IDENTIFICACIÓN A BORDO`) 


A_beatrix <- sponge_sampler_metadata %>% filter(grepl("Aphrocallistes", `IDENTIFICACIÓN A BORDO`)) %>% select(ID)
Suberites <- sponge_sampler_metadata %>% filter(grepl("Suberites", `IDENTIFICACIÓN A BORDO`)) %>% select(ID)
Thenea_muricata <- sponge_sampler_metadata %>% filter(grepl("Thenea muricata", `IDENTIFICACIÓN A BORDO`)) %>% select(ID)
Thenea_sp <- sponge_sampler_metadata %>% filter(grepl("Thenea sp", `IDENTIFICACIÓN A BORDO`)) %>% select(ID)
Phakellia <- sponge_sampler_metadata %>% filter(grepl("Phakellia", `IDENTIFICACIÓN A BORDO`)) %>% select(ID)
AXinella <- sponge_sampler_metadata %>% filter(grepl("Axinella", `IDENTIFICACIÓN A BORDO`)) %>% select(ID)
Pheronema <- sponge_sampler_metadata %>% filter(grepl("Pheronema", `IDENTIFICACIÓN A BORDO`)) %>% select(ID)

ASV_A_beatrix <- ASV_by_method %>% filter(Site %in% A_beatrix$ID)%>% filter(!grepl("Aphrocallistes", genus))
ASV_Suberites <- ASV_by_method %>% filter(Site %in% Suberites$ID) %>% filter(!grepl("Suberites", genus))
ASV_Thenea_muricata <- ASV_by_method %>% filter(Site %in% Thenea_muricata$ID)%>% filter(!grepl("Thenea", genus))
ASV_Thenea_sp <-  ASV_by_method %>% filter(Site %in% Thenea_sp$ID) %>% filter(!grepl("Thenea", genus))
ASV_Phakellia <-  ASV_by_method %>% filter(Site %in% Phakellia$ID) %>% filter(!grepl("Phakellia", genus))
ASV_Axinella <- ASV_by_method %>% filter(Site %in% AXinella$ID) %>% filter(!grepl("Axinellidae", family))
ASV_Pheronema <- ASV_by_method %>% filter(Site %in% Pheronema$ID) %>% filter(!grepl("Pheronema", genus))


Axinella_tb <- ASV_Axinella %>% group_by(Site, phylum, class, order, family, genus, species, Method) %>%
  summarise(per_nReads = sum(nReads), .groups = "drop") %>%
  mutate(species = fct_inorder(species)) %>%
  arrange(phylum, species) %>% mutate("prop_nReads" = per_nReads/sum(per_nReads)) %>% mutate("Sponge" = "Axinella")

Aphrocallistes_tb <- ASV_A_beatrix %>% group_by(Site, phylum, class, order, family, genus, species, Method) %>%
  summarise(per_nReads = sum(nReads), .groups = "drop") %>%
  mutate(species = fct_inorder(species)) %>%
  arrange(phylum, species) %>% mutate("prop_nReads" = per_nReads/sum(per_nReads)) %>% mutate("Sponge" = "Aphrocallistes")

Suberites_tb <- ASV_Suberites %>% group_by(Site, phylum, class, order, family, genus, species, Method) %>%
  summarise(per_nReads = sum(nReads), .groups = "drop") %>%
  mutate(species = fct_inorder(species)) %>%
  arrange(phylum, species) %>% mutate("prop_nReads" = per_nReads/sum(per_nReads)) %>% mutate("Sponge" = "Suberites")

Thenea_muricata_tb <- ASV_Thenea_muricata %>% group_by(Site, phylum, class, order, family, genus, species, Method) %>%
  summarise(per_nReads = sum(nReads), .groups = "drop") %>%
  mutate(species = fct_inorder(species)) %>%
  arrange(phylum, species) %>% mutate("prop_nReads" = per_nReads/sum(per_nReads)) %>% mutate("Sponge" = "Thenea muricata")

Thenea_sp_tb <- ASV_Thenea_sp %>% group_by(Site, phylum, class, order, family, genus, species, Method) %>%
  summarise(per_nReads = sum(nReads), .groups = "drop") %>%
  mutate(species = fct_inorder(species)) %>%
  arrange(phylum, species) %>% mutate("prop_nReads" = per_nReads/sum(per_nReads)) %>% mutate("Sponge" = "Thenea sp.")

Phakellia_tb <- ASV_Phakellia %>% group_by(Site, phylum, class, order, family, genus, species, Method) %>%
  summarise(per_nReads = sum(nReads), .groups = "drop") %>%
  mutate(species = fct_inorder(species)) %>%
  arrange(phylum, species) %>% mutate("prop_nReads" = per_nReads/sum(per_nReads)) %>% mutate("Sponge" = "Phakellia")

Pheronema_tb <- ASV_Pheronema %>% group_by(Site, phylum, class, order, family, genus, species, Method) %>%
  summarise(per_nReads = sum(nReads), .groups = "drop") %>%
  mutate(species = fct_inorder(species)) %>%
  arrange(phylum, species) %>% mutate("prop_nReads" = per_nReads/sum(per_nReads)) %>% mutate("Sponge" = "Pheronema")


Comparing_samplers <- rbind(Axinella_tb, Aphrocallistes_tb, Suberites_tb, Thenea_sp_tb, Thenea_muricata_tb, Phakellia_tb, Pheronema_tb)


Comparing_samplers %>% 
  group_by(Sponge, phylum) %>% 
  ggplot(aes(x = Sponge, y = prop_nReads, fill = phylum)) +
  geom_col() -> sponges_plot_phyla

Comparing_samplers %>% 
  group_by(Sponge, genus) %>% 
  ggplot(aes(x = Sponge, y = prop_nReads, fill = genus)) +
  geom_col() + 
  theme(legend.position = "none")-> sponges_plot_genera

Comparing_samplers %>% 
  group_by(Sponge) %>% 
  ggplot(aes(x = Sponge, y = n_distinct(species)/n_distinct(Site))) +
  geom_col() -> sponges_species_count_per_sample


Comparing_samplers_with_lance <- Comparing_samplers %>% left_join(metadata2, by = "Site")

Comparing_samplers <- Comparing_samplers_with_lance


Comparing_samplers %>% 
  group_by(LANCE, Sponge) %>% 
  ggplot(aes(x= LANCE, y = n_distinct(species))) +
  geom_col() -> sponge_perf_per_lance



#########comparing samplers with eDNA 


ASV_filteronly <- ASV_bothlances %>% filter(Method == "eDNA") %>% mutate("LANCE" = Site) %>% mutate("Sponge" = "eDNA") 

Filters_tb <- ASV_filteronly %>% group_by(Site, phylum, class, order, family, genus, species, Method) %>%
  summarise(per_nReads = sum(nReads), .groups = "drop") %>%
  mutate(species = fct_inorder(species)) %>%
  arrange(phylum, species) %>% mutate("prop_nReads" = per_nReads/sum(per_nReads)) %>% mutate("Sponge" = "eDNA") %>% mutate("LANCE" = Site)

comparing_samplers_with_filters <- rbind(Comparing_samplers, Filters_tb)

comparing_samplers_with_filters %>% 
  group_by(LANCE, Sponge) %>% 
  ggplot(aes(x= LANCE, y = n_distinct(species))) +
  geom_col() -> all_perf_per_lance

comparing_samplers_with_filters %>% filter(LANCE %in% ASV_bothlances$Site) %>% filter(LANCE %in% Comparing_samplers$LANCE) %>% 
   group_by(Sponge, LANCE) %>% 
  ggplot(aes(x = LANCE, y = n_distinct(species), fill = Sponge)) +
  geom_col() -> methods_across_lances




######venn diagram time


library(tidyverse)
library(ggVennDiagram)

# Read in the data if not already loaded
# df <- read.csv("your_file.csv") 

# Create a list where each sponge category contains unique species
species_by_sponge <- comparing_samplers_with_filters %>% filter(LANCE %in% ASV_bothlances$Site) %>% filter(LANCE %in% Comparing_samplers$LANCE) %>% 
  group_by(Sponge) %>%
  summarise(species_list = list(unique(phylum))) %>%
  deframe() # Convert to named list

# Generate Venn diagram
ggVennDiagram(species_by_sponge) + 
  scale_fill_gradient(low = "white", high = "blue", ) +
  theme_minimal()

comparingL10 <- comparing_samplers_with_filters %>% filter(LANCE == "L10") %>% 
  group_by(Sponge) %>%
  summarise(species_list = list(unique(species))) %>%
  deframe() # Convert to named list

comparingL18 <- comparing_samplers_with_filters %>% filter(LANCE == "L18") %>% 
  group_by(Sponge) %>%
  summarise(species_list = list(unique(species))) %>%
  deframe() # Convert to named list

comparingL19 <- comparing_samplers_with_filters %>% filter(LANCE == "L19") %>% 
  group_by(Sponge) %>%
  summarise(species_list = list(unique(species))) %>%
  deframe() # Convert to named list


# Generate Venn diagram without percentages
ggVennDiagram(comparingL10, label = "count") +  # Removes percentages
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal()


comparing_sponges_only <- Comparing_samplers %>% 
  group_by(Sponge) %>%
  summarise(species_list = list(unique(species))) %>%
  deframe() # Convert to named list

ggVennDiagram(comparing_sponges_only, label = "count") +  # Removes percentages
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal()

#by sites
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
ggvenn(
  comparingL10 , 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
