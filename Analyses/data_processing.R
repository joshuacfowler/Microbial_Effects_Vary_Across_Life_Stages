# Title: How commonly do microbial effects vary across host life stages? #####
# Purpose: Reads in data from literature search, cleans up the data to keep only those that meet criteria, and assigns taxonomy
# Authors: Josh Fowler and Gwen Pohlmann #####
# Date: Jan 27, 2025 #####

# library(renv)
# renv::init()
# renv::snapshot()
# renv::deactivate()


library(tidyverse)


library(rotl)
library(ape)


####### Reading in the data   #######
# This data is stored in Teams; we have downloaded the most recent version to a local directory as of Sep 17, 2024


joshpath <- c("~/Dropbox/Microbial_Effects_Metaanalysis/")


# gwen wd
# setwd("~/Desktop/afkhami_lab/meta_analysis/R/raw_data")
gwenpath <- c("./")


path <- joshpath
# path <- gwenpath

raw_effects_df <- read_csv(file = paste0(path,("Microbial Effects Literature Search(Effect_sizes).csv"))) %>% 
  # raw_effects_df <- read_csv(file = paste0(path,("20241012_effect_sizes.csv"))) %>% 
  filter(!is.na(mean_symbiotic)) %>% 
  mutate(across(mean_symbiotic:n_aposymbiotic, as.numeric))



# find out how many distinct studies we have extracted data from
length(unique(raw_effects_df$study_number))
length(unique(filter(raw_effects_df, !is.na(se_symbiotic))$study_number)) # number of studies if we drop the ones that don't have SD (probably some of these have SE but not SD, but still)



variance_RII <- function(Bw, Bo, SDw, SDo, Nw, No){
  # calculating the variance of RII following formula from from Armas et al. 2004 supplement
  VARw <- (SDw^2)
  VARo <- (SDo^2)
  first_term = ((VARw/Nw) + (VARo/No))/((Bw+Bo)^2)
  second_term = (Bw-Bo)^2/(Bw+Bo)^2
  rho = ((VARw/Nw)-(VARo/No))/((VARw/Nw)+(VARo/No))
  V_rii = first_term*(1+second_term - ((2*rho*(Bw-Bo))/(Bw+Bo)))
  return(V_rii)
}

# calculate effects sizes and add to the data frame, also clean up the data frame
invalid_genera <- c("AMF","NAB", "Endophytic", "DAXY0016C","DYXY033","DYSH004","DYXY023","DYXY013C","DYXY003","DYXY004","DYXY001","DYXYY2","DYXYXY1","DYXY002","DYXY111","DYXY112")
# REVISIT: invalid genera ==== 
# come back to the "genera" above and manually determine the correct genus. for now, we omit them from analysis
effects_calc_df <- raw_effects_df %>% 
  mutate(
    calc_sd_symbiotic = case_when((is.na(sd_symbiotic)) & (!is.na(se_symbiotic)) & (!is.na(n_symbiotic)) ~ (se_symbiotic*(sqrt(n_symbiotic))), TRUE ~ sd_symbiotic),
    calc_sd_aposymbiotic = case_when((is.na(sd_aposymbiotic)) & (!is.na(se_aposymbiotic)) & (!is.na(n_aposymbiotic)) ~ (se_aposymbiotic*(sqrt(n_aposymbiotic))), TRUE ~ sd_aposymbiotic),
    calc_se_symbiotic = case_when((is.na(se_symbiotic)) & (!is.na(sd_symbiotic)) & (!is.na(n_symbiotic)) ~ (sd_symbiotic/(sqrt(n_symbiotic))), TRUE ~ se_symbiotic),
    calc_se_aposymbiotic = case_when((is.na(se_aposymbiotic)) & (!is.na(sd_aposymbiotic)) & (!is.na(n_aposymbiotic)) ~ (sd_aposymbiotic/(sqrt(n_aposymbiotic))), TRUE ~ se_aposymbiotic)
  ) %>%
  mutate(
    RII  = (mean_symbiotic - mean_aposymbiotic)/(mean_aposymbiotic + mean_symbiotic),
    cohensD = (mean_symbiotic - mean_aposymbiotic)/sqrt((calc_sd_aposymbiotic^2 + calc_sd_symbiotic^2)/2),
    lRR = log(mean_symbiotic/mean_aposymbiotic),
    var_RII =variance_RII(Bw = mean_symbiotic, Bo = mean_aposymbiotic, SDw = calc_sd_symbiotic, SDo = calc_sd_aposymbiotic, Nw = n_symbiotic, No = n_aposymbiotic),
    sd_RII = sqrt(var_RII),
    pooled_sd = sqrt((sd_aposymbiotic^2 + sd_symbiotic^2)/2),
    pooled_se = pooled_sd/sqrt(n_aposymbiotic + n_symbiotic)) %>% #I will used the studies pooled se as part of the measurement-error model in BRMS, rather than using hedge's G
  #hedgesG = (mean_aposymbiotic - mean_symbiotic)/sqrt((n_aposymbiotic-1)*sd_aposymbiotic)
  mutate(treatment_label = paste(study_number, experiment_id, treatment_id, sep = "-")) %>%
  mutate(experiment_label = paste(study_number, experiment_id, sep = "-")) %>% 
  mutate(lifestage_general = case_when(startsWith(lifestage_description, "juvenile") ~ "juvenile",
                                       startsWith(lifestage_description, "embryo") ~ "embryo",
                                       TRUE ~ lifestage_description)) %>% 
  mutate(symbiont_genus = word(symbiont_species, 1)) %>% 
  mutate(symbiont_genus_clean = 
           case_when((symbiont_genus == "E.") | (symbiont_genus == "Epichloe\xa8") | (symbiont_genus == "Epichlo\xeb") ~ "Epichloe",
                     (symbiont_genus %in% invalid_genera) ~ NA,
                     TRUE ~ symbiont_genus)) %>%
  mutate(symbiont_species_clean = paste(symbiont_genus_clean, word(symbiont_species, 2), sep = " ")) %>%
  mutate(host_genus = word(host_species, 1),
         host_genus = case_when(host_genus == "Schedonorus" ~ "Lolium", TRUE ~ host_genus)) %>% # having issues with schedonorus. it has flag "barren" which OTL says means there are only higher taxa at and below this node, no species or unranked tips
  mutate(host_species_clean = paste(host_genus, word(host_species, 2), sep = " ")) %>% 
  filter(!is.na(sd_symbiotic) | !is.na(calc_sd_symbiotic))


# summarizing how many experiments across the studies
length(unique(effects_calc_df$experiment_label))

symbiota_summary <- effects_calc_df %>% 
  mutate(symbiota_label = paste(symbiont_species_clean, host_species_clean, sep = "_"))
length(unique(symbiota_summary$symbiota_label))
length(unique(effects_calc_df$host_genus))
length(unique(effects_calc_df$symbiont_genus_clean))
######## trying out the rotl package #########

host_genera = array(unique(effects_calc_df$host_genus))
# having issues with schedonorus. it has flag "barren" which OTL says means there are only higher taxa at and below this node, no species or unranked tips
# host_genera[host_genera == "Schedonorus"] <- "Lolium"



# matching the names we have to the ott_ids in OTL
host_genera_names = tnrs_match_names(host_genera) 

# we can inspect and update those taxa which have multiple matches. Fortunately (as of Jan 28, 2025) the function chose the correct taxa in all cases. 
mult_matches = subset(host_genera_names, number_matches > 1)


inspect(host_genera_names, taxon_name = "phaseolus")


# host_genera_names <- update(host_genera_names,
#                          taxon_name = "setaria",
#                          new_row_number = 2
# )

# test_tree = tol_induced_subtree(ott_ids = host_genera_names$ott_id)
# plot(test_tree, show.tip.label = FALSE)
# 
# host_cov_matrix <- vcv(test_tree)

# get a slightly cleaner version of the genus info
host_genera_ott <- host_genera_names %>% 
  mutate(genus_ott_id = as.character(ott_id),
         host_genus_name = as.character(unique_name)) %>% select(genus_ott_id, host_genus_name, search_string)


# and then connect this to the full taxonomy
host_taxonomy <- tax_lineage(taxonomy_taxon_info(host_genera_names$ott_id, include_lineage = TRUE))


host_taxonomy_df <- bind_rows(host_taxonomy, .id = "genus_ott_id") %>% 
  left_join(host_genera_ott, by = join_by(genus_ott_id)) %>% 
  select(-unique_name, -ott_id) %>%
  # group_by(genus_ott_id) %>% 
  # mutate(rank_number = row_number()) %>% 
  filter(rank != "no rank") %>% 
  pivot_wider(names_from = c(rank), values_from = c(name), names_prefix = "host_", values_fn = list) %>% 
  select(genus_ott_id, host_genus_name, search_string, host_domain, host_kingdom, host_phylum, host_class, host_order, host_family ) %>% 
    unnest(cols = everything()) %>% 
  filter(host_kingdom != "Archaeplastida") # I don't know exactly why, but the lineage function pulls out a different set of kingdom/phylum/class for some of the species
  
unique(host_taxonomy_df$host_class)
unique(host_taxonomy_df$host_family)
unique(host_taxonomy_df$host_kingdom)

  



# Joining the effect sizes with taxonomy


effects_df <- effects_calc_df %>% 
  mutate(search_string = tolower(host_genus)) %>% 
  left_join(host_taxonomy_df)

write_csv(effects_df, "effects_df.csv")  



