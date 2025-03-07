

#making new ASV table 

ASV_nsDNA <- species.table.without.pos %>% filter(grepl("-", Sample))
ASV_eDNA <-species.table.without.pos %>% filter(!grepl("-", Sample))

ASV_eDNA <- ASV_eDNA %>% separate(col = Sample, into = c("Site", "Replicat"), sep = "_")
ASV_nsDNA <- ASV_nsDNA %>% separate(col = Sample, into = c("Site", "Replicat"), sep = "_")

ASV_eDNA <-ASV_eDNA %>% mutate(Method = "eDNA")
ASV_nsDNA <- ASV_nsDNA %>% mutate(Method = "nsDNA")

ASV_by_method <- rbind(ASV_eDNA, ASV_nsDNA)

write_csv(ASV_by_method, "ASV_by_method.csv")





#BRINGING IN MORE METADATA
library(readxl)

metadata2 <- read_excel("P23_Muestras_EsponjasFiltros_eDNA_AnaRiesgo_FULLMETADATA.xlsx", sheet = "Esponjas" )

metadata2 <- metadata2 %>% select(ID,LANCE) %>% mutate(LANCE = (paste("L",LANCE, sep = ""))) 

metadata2  <- metadata2 %>% as.data.frame()
metadata2 <- metadata2  %>% dplyr::rename(Site = ID)

#########################################
#                                       #
#   LECTURAS POR ESTACION Y LOCALIDAD   #
#                                       #
#########################################

ASV <-read.csv(here("Outputs/good_output/ASV_tables/ASV_table_samples_metazoa_without_pos_host.csv"))

ASV_sum <- ASV_by_method %>% 
  group_by(Site, Method) %>% 
  dplyr::summarise(sum(nReads)) %>% 
   mutate(Method=factor(Method,levels = c("eDNA", "nsDNA"))) 

ASV_sum$nReads <- ASV_sum$"sum(nReads)"
ASV_sum <- ASV_sum[, -3]

 color <- c('#bebada','#fb8072','#8dd3c7','#ffffb3','#bebada','#fb8072','#8dd3c7','#ffffb3','#bebada','#fb8072','#8dd3c7','#ffffb3','#bebada','#fb8072','#8dd3c7','#ffffb3','#bebada','#fb8072','#8dd3c7','#ffffb3','#bebada','#fb8072','#8dd3c7','#ffffb3','#bebada','#fb8072','#8dd3c7','#ffffb3','#bebada','#fb8072')

ASV_sum %>% 
  ggplot() +
  geom_bar(mapping = aes(x = Site, y = nReads, fill = Method), 
           stat = "identity", position = "stack") +
  scale_fill_manual(values = color)+
  labs(y = "Abundance of reads",
       x = ""
  ) +
  #  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# 1.1. Vamos a hacer una grafica de caja con bigotes para ambas medias

ASV_kruskal <- ASV_sum %>% 
  mutate(Loc = case_when(
    str_detect(Station, "Ri") ~ "Ri",
    str_detect(Station, "Urr") ~ "Ur",
    TRUE ~ NA_character_
  )) %>%
  mutate(Season = case_when(
    str_detect(Station, "_P") ~ "Spring",
    str_detect(Station, "_V") ~ "Summer",
    TRUE ~ NA_character_
  )) %>% 
  select(Loc,Season,Station,Tube,nReads)

#Ploteamos el histograma a ver como es la distribucion

ggplot(data = ASV_kruskal, aes(x = nReads)) +
  geom_histogram(binwidth = 1000, position = "dodge") +
  labs(x = "Number of Reads", y = "Frequency") +
  scale_fill_discrete(name = "Station") +
  theme_minimal()


ASV_kruskal <- ASV_kruskal %>% 
  mutate(lognReads = log(nReads)) %>% 
  mutate(sqrtnReads = sqrt(nReads))


ggplot(data = ASV_kruskal, aes(x = lognReads)) +
  geom_histogram(binwidth = 1000, position = "dodge") +
  labs(x = "Number of Reads", y = "Frequency") +
  scale_fill_discrete(name = "Station") +
  theme_minimal()


# 1.2. Nuestros datos siguen una distribucion normal?, Si es significativo es que NO sigue una distribucion normal


resultado_shapiro <- shapiro.test(ASV_sum$nReads) 
resultado_shapiro <- shapiro.test(ASV_kruskal$lognReads[ASV_kruskal$Loc == "Ur"]) #No es significativo por lo que ya sigue una distribucion normal
resultado_shapiro <- shapiro.test(ASV_kruskal$sqrtnReads[ASV_kruskal$Loc == "Ur"])


bartlett.test( lognReads$ASV_kruskal )

# 1.3. Comparamos medias a traves del test de Wilcoxon

wilcox.test(lognReads ~ Loc, data=ASV_kruskal) # SI hay diferencias significativas entre localidades
wilcox.test(lognReads ~ Season, data=ASV_kruskal) # NO hay diferneicas significativas entre estaciones

# 1.4. ¿Y entre la interaccion de estacion y localidad?

ASV_kruskal <- ASV_kruskal %>% 
  mutate(kruskal = paste(Loc,Season,sep="_"))

aggregate(lognReads ~ kruskal, data = ASV_kruskal, FUN = median)
aggregate(lognReads ~ kruskal, data = ASV_kruskal, FUN = sd)

ggplot(data = ASV_kruskal, mapping = aes(x = kruskal, y = nReads, colour = kruskal)) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = "none")

ggsave(here("Outputs/good_output/figuras/bubble_plot.pdf"), plot = mi_ggplot)

ggplot(data = ASV_kruskal, mapping = aes(x = lognReads, colour = kruskal)) +
  geom_histogram() +
  theme_bw() +
  facet_grid(. ~ kruskal) +
  theme(legend.position = "none")

library(car)
leveneTest(sqrtnReads ~ kruskal, data = ASV_kruskal, center = "median")

kruskal.test(sqrtnReads ~ kruskal, data = ASV_kruskal)

pairwise.wilcox.test(x = ASV_kruskal$lognReads, g = ASV_kruskal$kruskal, p.adjust.method = "holm" )





#######BUBBBLE PLOT FOR REAL


library(ggplot2)
library(dplyr)
library(forcats) # For fct_inorder()


ASV_by_method <- read.csv("ASV_by_method.csv")

# Transform the table
table_tb <- ASV_by_method %>%
  group_by(Site, species, phylum, Method) %>%
  summarise(per_nReads = sum(nReads), .groups = "drop") %>%
  mutate(species = fct_inorder(species)) %>%
  arrange(phylum, species)

table_tb_sponges <- table_tb %>% filter(grepl("-", Site)) %>% left_join(y = metadata2, by = "Site")
table_tb_nosponges <- table_tb %>% filter(!grepl("-", Site)) 

# table_tb_sponges <- table_tb_sponges %>% filter(!grepl("ROV", LANCE)) %>% filter(!grepl("NA", LANCE)) 

table_tb_sponges <- table_tb_sponges %>% select(!Site) %>% 
  dplyr::rename("Site" = "LANCE") %>% 
  select(Site, species, phylum, Method, per_nReads) 

table_tb_2 <- rbind(table_tb_nosponges, table_tb_sponges)

table_tb_2 <- table_tb_2 %>% filter(Site %in% intersect(table_tb_nosponges$Site, table_tb_sponges$Site))
table_tb_2 <- table_tb_nosponges %>% filter(!Site %in% table_tb_sponges$Site) %>% rbind (table_tb_sponges)

setdiff(unique(table_tb_nosponges$Site), unique(table_tb_sponges$Site)) 

setdiff(unique(table_tb_sponges$Site), unique(table_tb_nosponges$Site)) 

setequal(unique(table_tb_nosponges$Site), unique(table_tb_sponges$Site))  

shared_lances <- intersect(unique(table_tb_nosponges$Site), unique(table_tb_sponges$Site))

table_tb_3_onlyshared <- table_tb_2 %>% filter(Site %in% shared_lances)

library(tidyverse)


# Ensure correct ordering and visualization
table_tb %>%
  arrange(phylum) %>% 
  mutate(species = fct_inorder(species)) %>%
  mutate(Method = factor(Method, levels = unique(Method))) %>%  # Ensure correct factor levels
  arrange(Method) %>% 
  ggplot(aes(x = Site, y = species, size = per_nReads)) +  # Use 'Site' instead of 'Sample'
  geom_point(aes(fill = Method), alpha = 0.8, pch = 21, color = "black") +
  scale_fill_brewer(palette = "Accent") +
  scale_size(range = c(3, 11), name = "Proportion (%)") +
  facet_grid(phylum ~ Method, scales = "free", space = "free") +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0, size = 7),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 5),
    axis.text.y.left = element_text(angle = 0, size = 7)
  ) -> mi_ggplot



table_tb %>%
  arrange(phylum) %>% 
  mutate(species = fct_inorder(species)) %>%
  mutate(Method = factor(Method, levels = unique(Method))) %>%  # Ensure correct factor levels
  arrange(Method) %>% 
  ggplot(aes(x = Site, y = species, size = per_nReads)) +  # Use 'Site' instead of 'Sample'
  geom_point(aes(fill = Method), alpha = 0.8, pch = 21, color = "black") +
  scale_fill_brewer(palette = "Accent") +
  scale_size(range = c(3, 11), name = "Proportion (%)") +
  facet_grid(phylum ~ Method, scales = "free", space = "free") +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0, size = 7),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 5),
    axis.text.y.left = element_text(size = 5),  # Make species labels smaller
    panel.spacing = unit(0.8, "lines")  # Increase spacing between facets
  ) -> mi_ggplot2



table_tb %>%
  arrange(phylum) %>% 
  mutate(species = fct_inorder(species)) %>%
  mutate(Method = factor(Method, levels = unique(Method))) %>%  # Ensure correct factor levels
  arrange(Method) %>% 
  ggplot(aes(x = Site, y = species, size = per_nReads)) +  # Use 'Site' instead of 'Sample'
  geom_point(aes(fill = Method), alpha = 0.8, pch = 21, color = "black") +
  scale_fill_brewer(palette = "Accent") +
  scale_size(range = c(3, 11), name = "Proportion (%)") +
  facet_grid(phylum ~ Method, scales = "free", space = "free") +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0, size = 7),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 5),
    axis.text.y.left = element_text(size = 5),  # Make species labels smaller
    panel.spacing = unit(0.8, "lines")  # Increase spacing between facets
  ) -> mi_ggplot2



table_tb %>%
  arrange(phylum) %>% 
  mutate(species = fct_inorder(species)) %>%
  mutate(Method = factor(Method, levels = unique(Method))) %>%  # Ensure correct factor levels
  arrange(Method) %>% 
  ggplot(aes(x = Site, y = species, size = per_nReads)) +  
  geom_point(aes(fill = Method), alpha = 0.5, pch = 21, color = "black") +  # Increased transparency (alpha)
  scale_fill_brewer(palette = "Accent") +
  scale_size(range = c(3, 11), name = "Proportion (%)") +
  facet_grid(phylum ~ ., scales = "free", space = "free") +  # Remove Method from facetting
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0, size = 7),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 5),
    axis.text.y.left = element_text(size = 5),  
    panel.spacing = unit(0.8, "lines")  
  ) -> mi_ggplot3



table_tb_2 %>%
  group_by(phylum, Site, Method) %>%  # Aggregate at phylum level
  summarise(total_reads = sum(per_nReads, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = Site, y = phylum, size = total_reads)) +  
  geom_point(aes(fill = Method), alpha = 0.5, pch = 21, color = "black") +  
  scale_fill_brewer(palette = "Accent") +
  scale_size(range = c(3, 11), name = "Proportion (%)") +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0, size = 7),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 5),
    axis.text.y.left = element_text(size = 7),  # Adjusted for phylum labels
    panel.spacing = unit(0.8, "lines")  
  ) -> mi_ggplot3


table_tb_3_onlyshared %>%
  group_by(phylum, Method) %>%  # Aggregate only at phylum and method level
  summarise(total_reads = sum(per_nReads, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = Method, y = phylum, size = total_reads)) +  # Use Method instead of Site
  geom_point(aes(fill = Method), alpha = 0.5, pch = 21, color = "black") +  
  scale_fill_brewer(palette = "Accent") +
  scale_size(range = c(3, 11), name = "Proportion (%)") +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0, size = 7),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7),
    axis.text.y.left = element_text(size = 7),  # Adjusted for phylum labels
    panel.spacing = unit(0.8, "lines")  
  ) -> mi_ggplot4



# table_tb_3_onlyshared %>%
#   group_by(phylum, Method) %>%  # Aggregate at phylum and method level
#   summarise(species_count = n_distinct(species), .groups = "drop") %>%  # Count unique species
#   ggplot(aes(x = Method, y = phylum, size = species_count)) +  # Use species_count for bubble size
#   geom_point(aes(fill = Method), alpha = 0.5, pch = 21, color = "black") +  
#   scale_fill_brewer(palette = "Accent") +
#   scale_size(range = c(3, 11), name = "Species Count") +  # Change legend title
#   theme_bw() +
#   theme(
#     strip.text.y = element_text(angle = 0, size = 7),
#     axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7),
#     axis.text.y.left = element_text(size = 7),  # Adjusted for phylum labels
#     panel.spacing = unit(0.8, "lines")  
#   ) -> mi_ggplot5


library(ggplot2)
library(dplyr)

table_tb_3_onlyshared %>%
  group_by(phylum, Method) %>%  # Aggregate at phylum and method level
  summarise(species_count = n_distinct(species), .groups = "drop") %>%  # Count unique species
  ggplot(aes(x = Method, y = phylum, size = species_count)) +  # Use species_count for bubble size
  geom_point(aes(fill = species_count), alpha = 0.7, pch = 21, color = "black") +  
  scale_fill_gradient(low = "blue", high = "red") +  # Gradient color scale based on species_count
  scale_size(range = c(3, 15), name = "Species Count") +  # Change legend title
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0, size = 7),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7),
    axis.text.y.left = element_text(size = 7),  # Adjusted for phylum labels
    panel.spacing = unit(0.8, "lines")  
  ) -> mi_ggplot5

# To view the plot
print(mi_ggplot5)

table_tb_3_onlyshared %>%
  group_by(phylum, Method) %>%  # Aggregate only at phylum and method level
  summarise(total_reads = sum(per_nReads, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = Method, y = phylum, size = total_reads)) +  # Use species_count for bubble size
  geom_point(aes(fill = total_reads), alpha = 0.7, pch = 21, color = "black") +  
  scale_fill_gradient(low = "blue", high = "red") +  # Gradient color scale based on species_count
  scale_size(range = c(3, 15), name = "Proportion (%)") +  # Change legend title
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0, size = 7),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 7),
    axis.text.y.left = element_text(size = 7),  # Adjusted for phylum labels
    panel.spacing = unit(0.8, "lines")  
  ) -> mi_ggplot6



table_tb_3_onlyshared



library(ggplot2)
library(dplyr)




# Save the plot
ggsave(here("Figures/bubble_plot6.pdf"), plot = mi_ggplot, width = 7.5, height = 8.5)


ggsave(here("Outputs/bubble_plot_DIRTY2.pdf"), plot = mi_ggplot2, width = 15, height = 25)

ggsave(here("Outputs/bubble_plot_SPECIES_per_phylum_coinciding_lances.pdf"), plot = mi_ggplot5, width = 4, height = 7.5)

ggsave(here("Outputs/bubble_plot_READCOUNTS_proportion_per_phylum_coinciding_lances.pdf"), plot = mi_ggplot6, width = 4, height = 7.5)


ASV_by_method|> 
  group_by(Method) |> 
  summarise(across(c(class, family, species), n_distinct)) |>
  # inner_join(metadata, by="Sample") |>
  ggplot(aes(x = SIT, fill = phylum, y = family))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9)) +
  geom_col()

# Calculate average count and standard deviation for each site and method
summary_data <- table_tb_3_onlyshared %>%
  group_by(Method) %>%
  summarise(
    avg_count = mean(per_nReads, na.rm = TRUE),
    sd_count = sd(per_nReads, na.rm = TRUE),
    .groups = 'drop'
  )

# Create a bar plot with error bars
ggplot(summary_data, aes(x = Method, y = avg_count)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = avg_count - sd_count, ymax = avg_count + sd_count), 
                position = position_dodge(0.7), width = 0.25) +
  theme_minimal() +
  scale_fill_manual(values = c("eDNA" = "blue")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ASV_methods_Versus <- ASV_by_method %>% filter(intersect(table_tb_sponges$Site, table_tb_nosponges$Site)) %>% group_by(Method) %>% 
  summarise(avg_count = mean(nReads, na.rm = TRUE),
            sd_count = sd(nReads, na.rm = TRUE),
            .groups = 'drop'
  )


fixedASV <- left_join(ASV_by_method, metadata2, by= "Site")

fixedASV$Site[!is.na(fixedASV$LANCE)] <- fixedASV$LANCE[!is.na(fixedASV$LANCE)]

ASV_bothlances <- fixedASV %>% filter(Site %in% intersect(table_tb_sponges$Site, table_tb_nosponges$Site))

ASV_bothlances %>% filter(!grepl("Porifera", phylum)) %>% 
  group_by(Method, phylum) %>% 
  summarise(
    avg_count = mean(nReads, na.rm = TRUE),
    sd_count = sd(nReads, na.rm = TRUE),
    .groups = 'drop'
  ) %>% 
  ggplot(aes(x = Method, y = avg_count, fill = phylum)) +
  geom_col() ->no_por


ASV_bothlances %>% 
  group_by(Method, phylum) %>% 
  summarise(
    avg_count = mean(nReads, na.rm = TRUE),
    sd_count = sd(nReads, na.rm = TRUE),
    .groups = 'drop'
  ) %>% 
  ggplot(aes(x = Method, y = avg_count, fill = phylum)) +
  geom_col() ->with_por

ggsave("stacked_bar_phyla_eDNA_nsDNA_no_porifera.pdf", plot = no_por )

ggsave("stacked_bar_phyla_eDNA_nsDNA.pdf", plot = with_por )



ASV_bothlances %>% group_by(Site, Method) %>% 
  summarise(totalReads =sum(nReads), 
            .groups = 'drop') %>% 
  ggplot(aes(x = Method, y = totalReads)) +  
  geom_violin() + geom_point()-> comp_count

ASV_bothlances %>% filter(!grepl("Porifera", phylum)) %>%  group_by(Site, Method) %>% 
  summarise(totalReads =sum(nReads), 
            .groups = 'drop') %>% 
  ggplot(aes(x = Method, y = totalReads)) +  
  geom_violin() + geom_point()-> comp_count_nopor


# +  # Bar plot for average count
# geom_errorbar(aes(ymin = avg_count - sd_count, ymax = avg_count + sd_count), width = 0.2)

