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

library(rotl)
library(metafor)


####### Reading in the data   #######
# This data is stored in Teams; we have downloaded the most recent version to a local directory as of Sep 17, 2024


<<<<<<< HEAD
joshpath <- c("~/Dropbox/Microbial_Effects_Metaanalysis/")


# gwen wd
setwd("~/Desktop/afkhami_lab/meta_analysis/R/raw_data")
gwenpath <- c("./")
  
  
path <- joshpath
path <- gwenpath

raw_effects_df <- read_csv(file = paste0(path,("Microbial Effects Literature Search(Effect_sizes).csv"))) %>% 
# raw_effects_df <- read_csv(file = paste0(path,("20240912_effect_sizes.csv"))) %>% 
  filter(!is.na(mean_symbiotic)) %>% 
  mutate(across(mean_symbiotic:n_aposymbiotic, as.numeric))
  # separate_wider_delim(symbiont_species, delim = " ", names = c("symbiont_genus"), too_many = "align_start")
=======
# joshpath <- c("~/Dropbox/Microbial_Effects_Metaanalysis/")
# path <- joshpath


# gwen wd
setwd("~/Desktop/afkhami_lab/meta_analysis/R")


# raw_effects_df <- read_csv(file = paste0(path,("20240912_effect_sizes.csv"))) %>% 
raw_effects_df <- read.csv("./raw_data/20240919_effect_sizes.csv") %>% 
  mutate(across(mean_symbiotic:n_aposymbiotic, as.numeric)) %>%
  filter(!(is.na(mean_symbiotic & (sd_symbiotic | se_symbiotic))))
# REVISIT: filtering ==== 
# study 24 is filtering incorrectly. two rows where the se_symbiotic = 0 were not included in the df. i'm not sure why
>>>>>>> ac5967bc68f8c024b57cc96550fe5dae145d05cb

# find out how many distinct studies we have extracted data from
length(unique(raw_effects_df$study_number))

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
effects_df <- raw_effects_df %>% 
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
  mutate(host_genus = word(host_species, 1)) %>%
  mutate(host_species_clean = paste(host_genus, word(host_species, 2), sep = " "))




######## trying out the rotl package #########

host_genera = array(unique(effects_df$host_genus))
# having issues with schedonorus. it has flag "barren" which OTL says means there are only higher taxa at and below this node, no species or unranked tips
host_genera = host_genera[host_genera != "Schedonorus"]

# matching the names we have to the ott_ids in OTL
host_genera_names = tnrs_match_names(host_genera)
mult_matches = subset(host_genera_names, number_matches > 1)
inspect(host_genera_names, taxon_name = "prunella")
# REVISIT: host_genera_names ==== 
# need to manually check that all genus search strings are matched to the correct ott_ids

# making taxon map
taxon_map = structure(host_genera_names$search_string, names = host_genera_names$unique_name)
taxon_map["Tolumnia (genus in kingdom Archaeplastida)"]

# testing some phylogeny building
test_tree = tol_induced_subtree(ott_ids = host_genera_names$ott_id)
plot(test_tree, show.tip.label = FALSE)

test_tree2 = tol_induced_subtree(ott_id(host_genera_names)[is_in_tree(ott_id(host_genera_names))])
plot(test_tree2, show.tip.label = FALSE)

setdiff(test_tree$tip.label, test_tree2$tip.label)

otl_tips2 = strip_ott_ids(test_tree2$tip.label, remove_underscores = TRUE)
test_tree2$tip.label = taxon_map[otl_tips2]
plot(test_tree2)


# using rotl to check/filter host_species_clean column
host_names = array(unique(effects_df$host_species_clean))
# had to use the code below when making test_tree3 b/c Error: node_id 'ott3915043' was not found!list(ott3915043 = "pruned_ott_id")
#host_names = host_names[!(host_names == "Populus euramericana")]

host_species_names = tnrs_match_names(host_names)

find_mismatch_host = host_species_names %>% 
  mutate(search_epithet = word(search_string, 2)) %>% 
  mutate(otl_epithet = word(unique_name, 2)) %>% 
  mutate(mismatch = case_when((search_epithet == otl_epithet) ~ NA, TRUE ~ search_epithet))
which(!(is.na(find_mismatch_host$mismatch)))
# make sure to use the resulting row numbers from above as the values in c below
mismatch_rows_host = c(6,19,34,35,51)
mismatches_host = find_mismatch_host[mismatch_rows_host, ]
print(mismatches_host$search_string)
# REVISIT: host_species_names ====
# need to manually check that ott_ids are correctly matched by looking at the synonyms for the species below & referencing the papers
synonyms(host_species_names, taxon_name = "fragaria x")
synonyms(host_species_names, taxon_name = "schedonorus arundinaceus")
synonyms(host_species_names, taxon_name = "achnatherum sibiricum")
synonyms(host_species_names, taxon_name = "andropogon gerardii")
synonyms(host_species_names, taxon_name = "populus euramericana")

test_tree3 = tol_induced_subtree(ott_ids = host_species_names$ott_id)
plot(test_tree3, show.tip.label = FALSE)

# using rotl to check/filter symbiont taxonomy columns
symbiont_genera = array(unique(effects_df$symbiont_genus_clean))
symbiont_genera_names = tnrs_match_names(symbiont_genera)

symbiont_names = array(unique(effects_df$symbiont_species_clean))
symbiont_species_names = tnrs_match_names(symbiont_names)
invalid_species = c("Wolbachia wAv", "Wolbachia wBv", "Mycorrhizal fungi", "NA fungi", "Tulasnella strain", "NA diazotrophic", "NA (non-rhizobial", " E.", "Bradyrhizobium strain", "Epichloe (coexisting", " AMF", "NA bacteria", "Ceratobasidium B-3C-1", "Ceratobasidium B-4D-3", "Ceratobasidium B-4D-2", "Ceratobasidium Z-3A-3-2", "NA isolate")
# REVISIT: symbiont_names ====
# need to relabel the species above so they can be identified
symbiont_names_filtered = symbiont_names[!(symbiont_names %in% invalid_species)]
symbiont_species_names = tnrs_match_names(symbiont_names_filtered)

find_mismatch_symbio = symbiont_species_names %>%
  mutate(search_epithet = word(search_string, 2)) %>%
  mutate(otl_epithet = word(unique_name, 2)) %>%
  mutate(mismatch = case_when((search_epithet == otl_epithet) ~ NA, TRUE ~ search_epithet))
which(!(is.na(find_mismatch_symbio$mismatch)))
# make sure to use the resulting row numbers from above as the values in c below
mismatch_rows_symbio = c(8, 9, 10, 11, 12, 13, 15, 18, 23, 37, 38, 46, 49, 50, 52, 54, 55, 60)
mismatches_symbio = find_mismatch_symbio[mismatch_rows_symbio, ]
print(mismatches_symbio$search_string)
# REVISIT: symbiont_names_filtered ====
# need to manually check the search strings printed above by referencing the papers




######### plotting prelim data ########
ggplot(effects_df) +
  geom_histogram(aes(x = RII))+facet_wrap(~metric_category, scales = "free")

ggplot(effects_df) +
  geom_histogram(aes(x = lRR))+facet_wrap(~metric_category, scales = "free")

ggplot(effects_df) +
  geom_histogram(aes(x = cohensD))+facet_wrap(~metric_category, scales = "free")

ggplot(effects_df) +
  geom_histogram(aes(x = RII))+facet_wrap(~lifestage_general, scales = "free")

# prelim funnel plot
funnel(effects_df$RII, effects_df$var_RII, yaxis = "vi")





####### Fitting meta-analytic model  #######

# Testing out simpler models
fit <- lm(data = effects_df, formula = RII ~ 0 + metric_category)
fit <- lmer(data = effects_df, formula = RII ~ 0 + metric_category + (1|study_number))
fit <- lmer(data = effects_df, formula = RII ~ 0 + metric_category + (1|study_number) + (1|study_number:experiment_id) )
fit <- lmer(data = effects_df, formula = RII ~ 0 + metric_category + (1|study_number) + (1|experiment_label) )

# fit <- lmer(data = effects_df, formula = RII ~ 0 + metric_category + (1|study_number) + (1|experiment_label) + (1|treatment_label) + (1|species_label) +(metric_category|study_number))

summary(fit)
# summary(fit)$r.squared



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


# Fitting a model across life stages
fit <- lm(data = effects_df, formula = RII ~ 0 + lifestage_general)
fit <- lmer(data = effects_df, formula = RII ~ 0 + lifestage_description+ (1|study_number))
fit <- lmer(data = effects_df, formula = RII ~ 0 + lifestage_description + (1|study_number) + (1|study_number:experiment_id) )

summary(fit)


# Getting predictions from the model

prediction_df <- expand.grid( lifestage_general = unique(effects_df$lifestage_general),
                              study_number = NA,
                              experiment_id = NA)

preds <- predict(fit, newdat = prediction_df, se.fit = TRUE, re.form =  NA)

prediction_df <- bind_cols(prediction_df, preds) %>% 
  mutate(upr = fit + 1.96*se.fit,
         lwr = fit - 1.96*se.fit)


ggplot(data = prediction_df)+
  geom_jitter(data = effects_df, aes( x= lifestage_general, y = RII), color = "blue", width = .1, alpha = .2)+
  geom_point(aes(x = lifestage_general, y = fit), size = 3, alpha = .6) +
  geom_linerange(aes(x = lifestage_general, ymin = lwr, ymax = upr)) + theme_bw()


  
  # with Cohens D
effects_df_filtered <- effects_df %>% filter(!is.na(cohensD))
fit <- lmer(data = effects_df_filtered, formula = cohensD ~ 0 + metric_category + (1|study_number) + (1|study_number:experiment_id) )
summary(fit)
  
  
  # Getting predictions from the model
  
prediction_df <- expand.grid( metric_category = unique(effects_df$metric_category),
                              study_number = NA,
                              experiment_id = NA)  

preds <- predict(fit, newdat = prediction_df, se.fit = TRUE, re.form =  NA)

  prediction_df <- bind_cols(prediction_df, preds) %>% 
    mutate(upr = fit + 1.96*se.fit,
           lwr = fit - 1.96*se.fit)
  
  ggplot(data = prediction_df)+
    geom_jitter(data = effects_df_filtered, aes( x= metric_category, y = cohensD), color = "blue", width = .1, alpha = .2)+
    geom_point(aes(x = metric_category, y = fit), size = 3) +
    geom_linerange(aes(x = metric_category, ymin = lwr, ymax = upr)) + theme_bw()
  
  
  







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

# Version that does not incorporate measurement error
fit <- brm(formula = lRR~ 0 + metric_category + (1|study_number) + (1|experiment_label),
           data = effects_df, 
           family = "gaussian",
           prior = c(set_prior("normal(0,1)", class = "b"),
                     set_prior("student_t(3, 0, 2.5)", class = "sd")),
           iter = mcmc_pars$iter,
           chains = mcmc_pars$chains,
           warmup = mcmc_pars$warmup)


# version incorporating measurement error
effects_df_filtered <- effects_df %>% filter(!is.na(sd_RII))
# fit <- brm(formula = RII|se(sd_RII, sigma = TRUE) ~ 0 + metric_category + (1+metric_category|study_number) + (1|study_number:experiment_id),
fit <- brm(formula = RII|se(sd_RII, sigma = TRUE) ~ 0 + metric_category + (1|study_number) + (1|study_number:experiment_id),
           data = effects_df_filtered,
           family = "gaussian",
           prior = c(set_prior("normal(0,1)", class = "b"),
                     set_prior("student_t(3, 0, 2.5)", class = "sd")),
           iter = mcmc_pars$iter,
           chains = mcmc_pars$chains,
           warmup = mcmc_pars$warmup)

summary(fit)




# Version rescaling the RII data and fitting to a beta distribution
# rescale <- function(x){(x+1)/2}
# unscale <- function(x){(x*2)-1}
# 
# effects_rescaled = effects_df %>% 
#   mutate(RII_rescaled = rescale(RII))
# fit <- brm(formula = RII_rescaled~ 0 + metric_category + (1|study_number) + (1|experiment_label),
#            data = effects_rescaled, 
#            family = "zero_one_inflated_beta",
#            prior = c(set_prior("normal(0,1)", class = "b"),
#                      set_prior("student_t(3, 0, 2.5)", class = "sd")),
#            iter = mcmc_pars$iter,
#            chains = mcmc_pars$chains,
#            warmup = mcmc_pars$warmup)
# 
# 
# 
# # Version incorporating the measuremenerror through distributional regression
# beta_formula <- bf(RII_rescaled ~ 0 + metric_category + (1|study_number) + (1|experiment_label),
#               phi ~ 0 + metric_category + (1|study_number) + (1|experiment_label))
# fit <- brm(formula = beta_formula,
#            data = effects_rescaled, 
#            family = zero_one_inflated_beta(),
#            prior = c(set_prior("normal(0,1)", class = "b"),
#                      set_prior("student_t(3, 0, 2.5)", class = "sd")),
#            iter = mcmc_pars$iter,
#            chains = mcmc_pars$chains,
#            warmup = mcmc_pars$warmup)
# 
# # version incorporating measurement error
# fit <- brm(formula = RII|se(pooled_se) ~ 0 + metric_category + (1+|study_number) + (1|study_number:experiment_id),
#            data = effects_df,
#            family = "gaussian",
#            prior = c(set_prior("normal(0,1)", class = "b"),
#                      set_prior("student_t(3, 0, 2.5)", class = "sd")),
#             iter = mcmc_pars$iter,
#            chains = mcmc_pars$chains,
#            warmup = mcmc_pars$warmup)
# 
# summary(fit)


get_prior(fit)
prediction_df <- expand.grid( metric_category = unique(effects_df$metric_category),
                              study_number = NA,
                              experiment_id = NA,
                              sd_RII = 0)

preds <- fitted(fit, newdata = prediction_df, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), re_formula = NA)

prediction_df <- bind_cols(prediction_df, preds) #%>% 
  # mutate(across(Estimate:Q97.5,~unscale(.)))


ggplot(data = prediction_df)+
  geom_jitter(data = effects_df_filtered, aes( x= metric_category, y = RII, color = metric_category), width = .1, alpha = .2)+
  # geom_linerange(aes(x = metric_category, ymin = Q25, ymax = Q75), lwd = 1.2) + 
  geom_linerange(aes(x = metric_category, ymin = Q2.5, ymax = Q97.5)) + 
  geom_point(aes(x = metric_category, y = Estimate), size = 3) +
  theme_bw()




prediction_perstudy_df <- expand.grid( metric_category = unique(effects_df$metric_category),
                                       study_number = unique(effects_df_filtered$study_number),                              experiment_id = NA,
                              sd_RII = 0)

preds <- fitted(fit, newdata = prediction_perstudy_df, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), re_formula = ~(1|study_number))

prediction_perstudy_df <- bind_cols(prediction_perstudy_df, preds) %>% 
  group_by(metric_category) %>% 
  arrange(Estimate) %>% 
  mutate(ordered_id = row_number())

ggplot(data = prediction_perstudy_df)+
  # geom_linerange(aes(x = metric_category, ymin = Q25, ymax = Q75), lwd = 1.2) +
  geom_vline(aes(xintercept = 0))+
  geom_rect(data = prediction_df, aes(xmin=Q2.5, xmax=Q97.5, ymin=-Inf, ymax=Inf, fill = metric_category), alpha = .3)+
  geom_vline(data = prediction_df, aes(xintercept = Estimate, color = metric_category))+
  geom_linerange(aes(y = ordered_id, xmin = Q2.5, xmax = Q97.5), alpha = .7) +
  geom_point(aes(x = Estimate, y = ordered_id), size = 1, alpha = .7) +
  facet_wrap(~metric_category, scale = "free_y")+
  theme_bw()



