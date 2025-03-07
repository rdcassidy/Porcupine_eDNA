


#Solo necesito el combined_dataset (todos los hashes al 90/90; creo que podriamos hacer antes el filtrado por Metazoa y de pident de 95; HABLAR CON MONCHO)
#Conteamos para cada hash el numeor de diferentes valores para cada rank

#ROB: filtering to include only metazoans 
combined.dataset <- read_csv("Outputs/combined_dataset_to_modify_20_2_2025.csv")
combined.dataset_met <- combined.dataset %>% filter(kingdom == "Metazoa")

rank.count <- aggregate(cbind(phylum, class, order, family, genus, species) ~ qseqid, data = combined.dataset_met, FUN = function(x) c(count = length(unique(x))))

#Eliminamos las que han sido bien asignadas

multiple.time.assigned <- subset(rank.count, !(phylum == 1 & class == 1 & order == 1 & family == 1 & genus == 1 & species == 1))

#Queremos saber en que rango da multiples asignaciones

multiples.phylum <- multiple.time.assigned[multiple.time.assigned$phylum > 1, "qseqid"]
multiples.phylum <- data.frame(qseqid = multiples.phylum)

multiples.class <- multiple.time.assigned[multiple.time.assigned$class > 1, "qseqid"]
multiples.class <- data.frame(qseqid = multiples.class)
multiples.class <- anti_join(multiples.class, multiples.phylum)

multiples.order <- multiple.time.assigned[multiple.time.assigned$order > 1, "qseqid"]
multiples.order <- data.frame(qseqid = multiples.order)
multiples.order <- anti_join(multiples.order, multiples.phylum) %>%
  anti_join(multiples.class)

multiples.family <- multiple.time.assigned[multiple.time.assigned$family > 1, "qseqid"]
multiples.family <- data.frame(qseqid = multiples.family)
multiples.family <- anti_join(multiples.family, multiples.phylum) %>%
  anti_join(multiples.class) %>%
  anti_join(multiples.order)

multiples.genus <- multiple.time.assigned[multiple.time.assigned$genus > 1, "qseqid"]
multiples.genus <- data.frame(qseqid = multiples.genus)
multiples.genus <- anti_join(multiples.genus, multiples.phylum) %>%
  anti_join(multiples.class) %>%
  anti_join(multiples.order) %>%
  anti_join(multiples.family)

multiples.species <- multiple.time.assigned[multiple.time.assigned$species > 1, "qseqid"]
multiples.species <- data.frame(qseqid = multiples.species)
multiples.species <- anti_join(multiples.species, multiples.phylum) %>%
  anti_join(multiples.class) %>%
  anti_join(multiples.order) %>%
  anti_join(multiples.family) %>%
  anti_join(multiples.genus)

#Cogemos los qseqid y los miramos en la tabla que esta toda la info

to.explore.multiples.phylum <- semi_join(combined.dataset_met, multiples.phylum, by= "qseqid")
to.explore.multiples.class <- semi_join(combined.dataset_met, multiples.class, by= "qseqid")
to.explore.multiples.order <- semi_join(combined.dataset_met, multiples.order, by= "qseqid")
to.explore.multiples.family <- semi_join(combined.dataset_met, multiples.family, by= "qseqid")
to.explore.multiples.genus <- semi_join(combined.dataset_met, multiples.genus, by= "qseqid")
to.explore.multiples.species <- semi_join(combined.dataset_met, multiples.species, by= "qseqid")

#Conteamos el numero de casos para cada rango asignado 

ocurrencias.casos.unicos.phylum <- to.explore.multiples.phylum %>% 
  group_by(qseqid, phylum) %>% 
  summarise(count = n())

ocurrencias.casos.unicos.class <- to.explore.multiples.class %>% 
  group_by(qseqid, phylum, class) %>% 
  summarise(count = n())

ocurrencias.casos.unicos.order <- to.explore.multiples.order %>% 
  group_by(qseqid, phylum, class, order) %>% 
  summarise(count = n())

ocurrencias.casos.unicos.family <- to.explore.multiples.family %>% 
  group_by(qseqid, phylum, class, order, family) %>% 
  summarise(count = n())

ocurrencias.casos.unicos.genus <- to.explore.multiples.genus %>% 
  group_by(qseqid, phylum, class, order, family, genus) %>% 
  summarise(count = n())

ocurrencias.casos.unicos.species <- to.explore.multiples.species %>% 
  group_by(qseqid, phylum, class, order, family, genus, species) %>% 
  summarise(count = n())



