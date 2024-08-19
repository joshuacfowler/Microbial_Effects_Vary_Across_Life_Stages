# Title: How commonly do microbial effects vary across host life stages? #####
# Purpose: Reads in data from literature search and runs meta-analytic regression models
# Authors: Josh Fowler and Gwen Pohlmann #####
# Date: May 14, 2024 #####

library(renv)
# renv::init()


library(tidyverse)

library(effsize)

library(rstan)
library(brms)


#############################################################################
####### Reading in the data   #######
#############################################################################
# This data is stored in Teams; we have downloaded the most recent version to a local directory as of May 14, 2024
joshpath <- c("~/Dropbox/Microbial_Effects_Metaanalysis/")
# gwenpath <- 
  
  
path <- joshpath
# path <- gwenpath


raw_effects_df <- read_csv(file = paste0(path,("Microbial Effects Literature Search(Effect_sizes).csv"))) %>% 
  filter(!is.na(mean_symbiotic)) %>% 
  mutate(mean_symbiotic = as.numeric(mean_symbiotic),
         mean_aposymbiotic = as.numeric(mean_aposymbiotic))

length(unique(raw_effects_df$study_number))



# Calculating effects sizes

effects_df <- raw_effects_df %>% 
  mutate(
    RII  = (mean_symbiotic - mean_aposymbiotic)/(mean_aposymbiotic + mean_symbiotic),
    cohensD = ( mean_symbiotic-mean_aposymbiotic)/sqrt((sd_aposymbiotic^2 + sd_symbiotic^2)/2)) #hedgesG = (mean_aposymbiotic - mean_symbiotic)/sqrt((n_aposymbiotic-1)*sd_aposymbiotic)



ggplot(effects_df) +
  geom_histogram(aes(x = RII))+facet_wrap(~metric_category, scales = "free")





ggplot(effects_df) +
  geom_histogram(aes(x = cohensD))+facet_wrap(~lifestage_description, scales = "free")



#############################################################################
####### Fitting meta-analytic model  #######
#############################################################################

# Testing out simpler models
fit <- lm(data = effects_df, formula = RII ~ 0 + metric_category*host_category)
fit <- lm(data = effects_df, formula = RII ~ 0 + metric_category*host_category)

summary(fit)


# Getting predictions from the model

prediction_df <- expand.grid( metric_category = unique(effects_df$metric_category),
                              host_category = unique(effects_df$host_category))

preds <- predict(fit, newdata = prediction_df, type = "response", se.fit = TRUE)

prediction_df <- bind_cols(prediction_df, preds) %>% 
  mutate(upr = fit + 1.96*se.fit,
         lwr = fit - 1.96*se.fit)


ggplot(data = prediction_df)+
  geom_jitter(data = effects_df, aes( x= metric_category, y = RII), color = "blue", width = .1, alpha = .2)+
  geom_point(aes(x = metric_category, y = fit), size = 3, alpha = .6) +
  geom_linerange(aes(x = metric_category, ymin = lwr, ymax = upr)) + facet_wrap(~host_category) + theme_bw()



  
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

fit <- brm(formula = RII ~ metric_category,
           data = effects_df, 
            iter = mcmc_pars$iter,
           chains = mcmc_pars$chains,
           warmup = mcmc_pars$warmup)



