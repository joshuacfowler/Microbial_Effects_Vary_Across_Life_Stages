# Title: How commonly do microbial effects vary across host life stages? #####
# Purpose: Reads in data from literature search and runs meta-analytic regression models
# Authors: Josh Fowler and Gwen Pohlmann #####
# Date: May 14, 2024 #####

library(renv)
# renv::init()


library(tidyverse)

library(effsize) 

library(lme4)
library(rstan)
library(brms)


#############################################################################
####### Reading in the data   #######
#############################################################################
# This data is stored in Teams; we have downloaded the most recent version to a local directory as of Sep 6, 2024


# joshpath <- c("~/Dropbox/Microbial_Effects_Metaanalysis/")


# gwen wd
setwd("~/Desktop/afkhami_lab/meta_analysis/R/raw_data")
gwenpath <- c("./")
  
  
# path <- joshpath
path <- gwenpath


raw_effects_df <- read_csv(file = paste0(path,("20240906_effect_sizes.csv"))) %>% 
  filter(!is.na(mean_symbiotic)) %>% 
  mutate(mean_symbiotic = as.numeric(mean_symbiotic),
         mean_aposymbiotic = as.numeric(mean_aposymbiotic),
         n_symbiotic = as.numeric(n_symbiotic),
         n_aposymbiotic = as.numeric(n_aposymbiotic)
         )
  # separate_wider_delim(symbiont_species, delim = " ", names = c("symbiont_genus"), too_many = "align_start")


# find out how many distinct studies we have extracted data from
length(unique(raw_effects_df$study_number))


# calculate effects sizes and add to the data frame
# also calculate SD and SE where needed
# also create a new column to combine study_number and experiment_id and treatment_id
effects_df <- raw_effects_df %>% 
  mutate(
    RII  = (mean_symbiotic - mean_aposymbiotic)/(mean_aposymbiotic + mean_symbiotic),
    cohensD = (mean_symbiotic-mean_aposymbiotic)/sqrt((sd_aposymbiotic^2 + sd_symbiotic^2)/2)
    ) %>% #hedgesG = (mean_aposymbiotic - mean_symbiotic)/sqrt((n_aposymbiotic-1)*sd_aposymbiotic)
  mutate(
    calc_sd_symbiotic = case_when((is.na(sd_symbiotic)) & (!is.na(se_symbiotic)) & (!is.na(n_symbiotic)) ~ (se_symbiotic*(sqrt(n_symbiotic))), TRUE ~ NA),
    calc_sd_aposymbiotic = case_when((is.na(sd_aposymbiotic)) & (!is.na(se_aposymbiotic)) & (!is.na(n_aposymbiotic)) ~ (se_aposymbiotic*(sqrt(n_aposymbiotic))), TRUE ~ NA),
    calc_se_symbiotic = case_when((is.na(se_symbiotic)) & (!is.na(sd_symbiotic)) & (!is.na(n_symbiotic)) ~ (sd_symbiotic/(sqrt(n_symbiotic))), TRUE ~ NA),
    calc_se_aposymbiotic = case_when((is.na(se_aposymbiotic)) & (!is.na(sd_aposymbiotic)) & (!is.na(n_aposymbiotic)) ~ (sd_aposymbiotic/(sqrt(n_aposymbiotic))), TRUE ~ NA)
  ) %>%
  mutate(treatment_label = paste(study_number, experiment_id, treatment_id, sep = "-"))


# plotting prelim data
ggplot(effects_df) +
  geom_histogram(aes(x = RII))+facet_wrap(~metric_category, scales = "free")

ggplot(effects_df) +
  geom_histogram(aes(x = cohensD))+facet_wrap(~lifestage_description, scales = "free")



#############################################################################
####### Fitting meta-analytic model  #######
#############################################################################

# Testing out simpler models
fit <- lm(data = effects_df, formula = RII ~ 0 + metric_category)
fit <- lmer(data = effects_df, formula = RII ~ 0 + metric_category+ (1|study_number))
fit <- lmer(data = effects_df, formula = RII ~ 0 + metric_category + (1|study_number) + (1|experiment_id))
fit <- lmer(data = effects_df, formula = RII ~ 0 + metric_category + (1|treatment_label))

summary(fit)


# Getting predictions from the model

prediction_df <- expand.grid( metric_category = unique(effects_df$metric_category),
                              study_number = NA)

preds <- predict(fit, newdat = prediction_df, se.fit = TRUE, re.form =  NA)

prediction_df <- bind_cols(prediction_df, preds) %>% 
  mutate(upr = fit + 1.96*se.fit,
         lwr = fit - 1.96*se.fit)


ggplot(data = prediction_df)+
  geom_jitter(data = effects_df, aes( x= metric_category, y = RII), color = "blue", width = .1, alpha = .2)+
  geom_point(aes(x = metric_category, y = fit), size = 3, alpha = .6) +
  geom_linerange(aes(x = metric_category, ymin = lwr, ymax = upr)) + theme_bw()



  
  # with Cohens D
effects_df_filtered <- effects_df %>% filter(!is.na(cohensD))
  fit <- lm(data = effects_df_filtered, formula = cohensD ~ lifestage_description)
  summary(fit)
  
  
  # Getting predictions from the model
  
  prediction_df <- expand.grid( lifestage_description = unique(effects_df_filtered$lifestage_description))
  
  preds <- predict(fit, newdata = prediction_df, type = "response", se.fit = TRUE)
  
  prediction_df <- bind_cols(prediction_df, preds) %>% 
    mutate(upr = fit + 1.96*se.fit,
           lwr = fit - 1.96*se.fit)
  
  ggplot(data = prediction_df)+
    geom_jitter(data = effects_df_filtered, aes( x= lifestage_description, y = cohensD), color = "blue", width = .1, alpha = .2)+
    geom_point(aes(x = lifestage_description, y = fit), size = 3) +
    geom_linerange(aes(x = lifestage_description, ymin = lwr, ymax = upr)) + theme_bw()
  
  
  

# Version with random effects is giving a weird error

library(lme4)

fit <- lmer(data = effects_df, formula = RII ~ metric_category + (1|host_id))
summary(fit)






# fitting brms
## run this code to optimize computer system settings for MCMC
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
mcmc_pars <- list(
  iter = 5000, 
  warmup = 2500, 
  thin = 1, 
  chains = 3
)

fit <- brm(formula = RII ~ 0 + metric_category + (1|study_number),
           data = effects_df, 
           family = "gaussian",
           prior = c(set_prior("normal(0,5)", class = "b")))
            iter = mcmc_pars$iter
           chains = mcmc_pars$chains
           warmup = mcmc_pars$warmup
summary(fit)

prediction_df <- expand.grid( metric_category = unique(effects_df$metric_category),
                              study_number = 31)

preds <- fitted(fit, newdata = prediction_df, re_formula = NULL, allow_new_levels = TRUE)

prediction_df <- bind_cols(prediction_df, preds) 


ggplot(data = prediction_df)+
  geom_jitter(data = effects_df, aes( x= metric_category, y = RII), color = "blue", width = .1, alpha = .2)+
  geom_point(aes(x = metric_category, y = Estimate), size = 3) +
  geom_linerange(aes(x = metric_category, ymin = Q2.5, ymax = Q97.5)) + theme_bw()



