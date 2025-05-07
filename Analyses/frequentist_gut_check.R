# Title: Frequentist gut check #####
# Purpose: Reads in data from literature search and does 'simpler' statistical analysis of the meta-analytic dataset
# Authors: Josh Fowler  #####
# Date: May 6, 2025 #####

library(tidyverse)

library(effsize) 

library(lme4)


library(patchwork)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }
####### Reading in the data   #######
# The raw data is stored in Teams; we have downloaded the most recent version to a local directory as of May 6, 2025
# Data is imported, goes through some organization/cleaning, and we connect each study to its taxonomic information in the data_processing.R script

####### Reading in the data   #######
effects_df <- read_csv("effects_df.csv") %>% 
  filter(metric_category!="population metric") %>%   
  filter(!is.na(sd_RII))


####### Making taxa specific datasets for the most samples familes #######



asparagales <- effects_df %>% filter(host_order == "Asparagales")

poales <- effects_df %>% filter(host_order == "Poales")
# factorial anova

#######  performing sub analyses ####### 
# across all vital rates
vr_anova <- lm(formula = RII ~ 0+metric_category, data = effects_df)
summary(vr_anova)
vr_anova <- aov(formula = RII ~ metric_category, data = effects_df)
summary(vr_anova)

TukeyHSD(vr_anova)

# across all life stages
ls_anova <- lm(formula = RII ~ 0+lifestage_general, data = effects_df)
summary(ls_anova)
ls_anova <- aov(formula = RII ~ lifestage_general, data = effects_df)
summary(ls_anova)

TukeyHSD(ls_anova)


# within Asparagales

# across all vital rates
vr_anova <- lm(formula = RII ~ 0+metric_category, data = asparagales)
summary(vr_anova)
vr_anova <- aov(formula = RII ~ metric_category, data = asparagales)
summary(vr_anova)

TukeyHSD(vr_anova)

# across all life stages
ls_anova <- lm(formula = RII ~ 0+lifestage_general, data = asparagales)
summary(ls_anova)
ls_anova <- aov(formula = RII ~ lifestage_general, data = asparagales)
summary(ls_anova)

TukeyHSD(ls_anova)



# within Poales

# across all vital rates
vr_anova <- lm(formula = RII ~ 0+metric_category, data = poales)
summary(vr_anova)
vr_anova <- aov(formula = RII ~ metric_category, data = poales)
summary(vr_anova)

TukeyHSD(vr_anova)

# across all life stages
ls_anova <- lmer(formula = RII ~ 0+lifestage_general + (1|study_number) + (1|experiment_label), data = poales)
summary(ls_anova)
anova(ls_anova)
anova(ls_anova, test = "Chisq")
ls_anova <- aov(formula = RII ~ lifestage_general, data = poales)
summary(ls_anova)

TukeyHSD(ls_anova)



