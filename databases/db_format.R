library(ape)
library(stringr)
library(tidyverse)

seqs = read.FASTA("MIDORI2_DB_trimmed_derep.fasta")

fill_vector = function(x){
  superkingdom = grep("^superkingdom_",x)
  phylum = grep("^phylum_",x)
  class = grep("^class_",x)
  order = grep("^order_",x)
  family = grep("^family_",x)
  genus = grep("^genus_",x)
  species = grep("^species_",x)
  
  if (length(superkingdom) == 0){superkingdom = NA} 
  if (length(phylum) == 0){phylum = NA} 
  if (length(class) == 0){class = NA} 
  if (length(order) == 0){order = NA} 
  if (length(family) == 0){family = NA} 
  if (length(genus) == 0){genus = NA} 
  if (length(species) == 0){species = NA} 
  
  tax_pos = c(superkingdom, phylum, class, order, family, genus, species)
  tax_pos[is.na(tax_pos)] = tax_pos[which(is.na(tax_pos))-1] 
  tax_pos[is.na(tax_pos)] = tax_pos[which(is.na(tax_pos))-1] 
  tax_pos[is.na(tax_pos)] = tax_pos[which(is.na(tax_pos))-1] 
  tax_pos[is.na(tax_pos)] = tax_pos[which(is.na(tax_pos))-1] 
  tax_pos[is.na(tax_pos)] = tax_pos[which(is.na(tax_pos))-1] 
  tax_pos[is.na(tax_pos)] = tax_pos[which(is.na(tax_pos))-1] 
  tax_pos[is.na(tax_pos)] = tax_pos[which(is.na(tax_pos))-1] 
  
  new_x = x[tax_pos]
  return(new_x)
}


new_names = names(seqs) %>% 
  lapply(., FUN = function(x){str_split(x, pattern = "\t")[[1]][2]}) %>%
  unlist() %>%
  data.frame(names = .) %>%
  as_tibble() %>%
  mutate(ranks = map(names, ~ unlist(str_split(.x, pattern = ";")[[1]]) %>%
                       fill_vector()
                       )) %>%
  unnest_wider(ranks, names_sep = "_") %>%
  select(starts_with("rank")) %>%
  unite(everything,sep = ";") %>%
  as.list() %>% unlist()


names(seqs) = new_names
write.FASTA(seqs, "db_formatted.fas")

