# Pivoting literature search columns to get effect sizes for microbial effects meta-analysis 
# Josh Fowler and Gwen Pohlmann
# Mar 14, 2024

library(tidyverse)
library(readxl)


# read in the inital screening

lit_search <- read_excel("~/Downloads/Microbial Effects Literature Search.xlsx", sheet = "Lit_search") %>% 
  filter(include == "yes") %>% 
  separate_wider_delim(cols = host_taxa, delim = ",", names = paste0("host_species", 1:15), too_few = "align_start") %>% 
  separate_wider_delim(cols = symbiont_taxa, delim = ",", names = paste0("symbiont_species", 1:15), too_few = "align_start") %>% 
  pivot_longer(cols = host_species1:host_species15, names_to = "host_id", values_to = "host_species") %>% filter(!is.na(host_species)) %>% 
  pivot_longer(cols = symbiont_species1:symbiont_species15, names_to = "symbiont_id", values_to = "symbiont_species") %>% filter(!is.na(symbiont_species)) %>% 
  pivot_longer(cols = metric1:metric10, names_to = "metric_id", values_to = "metric_category") %>% 
  filter(!is.na(metric_category)) %>% 
  mutate(country = NA,
         metric_description = NA,
         life_stage_description = NA,
         mean_symbiotic = NA,
         mean_aposymbiotic = NA,
         sd_symbiotic = NA,
         sd_aposymbiotic = NA,
         samplesize_symbiotic = NA,
         samplesize_aposymbiotic = NA,
         drop_reason = NA) %>% 
  select(study_number, include_notes, drop_reason, host_id, host_species, symbiont_id, symbiont_species, multiple_stages, number_metrics, metrics_notes, metric_id, metric_category, metric_description, life_stage_description, mean_symbiotic, mean_aposymbiotic, sd_symbiotic, sd_aposymbiotic, samplesize_symbiotic, samplesize_aposymbiotic,country)


View(lit_search)

write_csv(lit_search, "MicrobialEffectsMetaAnalysis_effect_sizes_blank.csv")


write_csv(as_tibble(unique(lit_search$study_number)), "MicrobialEffectsMetaAnalysis_studies_as_of_2024-03-14.csv")
  
           


# Just for fun looking to see if species are in compadre database

lit_search_names <- lit_search %>% 
  mutate(host_taxa = str_replace(host_species, " ", "_")) %>% 
  select(host_species)%>% distinct() 


compadre_data <- load("~/Downloads/COMPADRE_v.6.23.5.0.RData")

compadre_meta <- compadre$metadata %>% 
  select(SpeciesAccepted) %>% distinct() %>% as.vector()


sapply(lit_search_names, function (y) sapply(compadre_meta, function (x) grepl(y, x)))

lit_search_names <- lit_search_names$host_species
string  <- paste0("(", paste(lit_search_names, collapse="|"), ")")

c("out", "in")[grepl(string, compadre_meta) + 1]

