# Title: How commonly do microbial effects vary across host life stages? #####
# Purpose: Reads in data from literature search and runs meta-analytic regression models
# Authors: Josh Fowler and Gwen Pohlmann #####
# test from gwen
# Date: May 14, 2024 #####

library(renv)
renv::init()


library(tidyverse)

library(effsize)

library(rstan)
library(brms)


#############################################################################
####### Reading in the data   #######
#############################################################################
# This data is stored in Teams; we have downloaded the most recent version to a local directory as of May 14, 2024
joshpath <- c("~/Dropbox/Microbial_Effects_Metaanalysis/")

path <- joshpath

raw_effects_df <- read_csv(file = paste0(path,("Microbial Effects Literature Search(Effect_sizes).csv"))) %>% 
  filter(!is.na(mean_symbiotic)) %>% 
  mutate(mean_symbiotic = as.numeric(mean_symbiotic),
         mean_aposymbiotic = as.numeric(mean_aposymbiotic))

length(unique(raw_effects_df$study_number))



# Calculating effects sizes

effects_df <- raw_effects_df %>% 
  mutate(
    RII  = (mean_aposymbiotic - mean_symbiotic)/(mean_aposymbiotic + mean_symbiotic),
    cohensD = (mean_symbiotic - mean_aposymbiotic)/sqrt((sd_aposymbiotic^2 + sd_symbiotic^2)/2))



ggplot(effects_df) +
  geom_histogram(aes(x = RII))+facet_wrap(~metric_category)





ggplot(effects_df) +
  geom_histogram(aes(x = cohensD))+facet_wrap(~metric_category)



#############################################################################
####### Fitting meta   #######
#############################################################################

