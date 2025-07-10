# Title: How commonly do microbial effects vary across host life stages? #####
# Purpose: Reads in data from literature search, evaluates publication bias, and visualized the data across different measurements
# Authors: Josh Fowler #####
# Date: March, 20, 2025 #####

# library(renv)
# renv::init()
# renv::snapshot()
# renv::deactivate()


library(tidyverse)

library(lme4)
library(rstan)
library(brms)


library(patchwork)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }
####### Reading in the data   #######
# This data is stored in Teams; we have downloaded the most recent version to a local directory as of Mar 18, 2025
# Data is imported, goes through some organization/cleaning, and we connect each study to its taxonomic information in the data_processing.R script

####### Reading in the data   #######
effects_df <- read_csv("effects_df.csv") %>% 
  filter(!is.na(sd_RII)) %>% 
  mutate(total_n = n_aposymbiotic+n_symbiotic,
         se_RII = sd_RII/sqrt(total_n))

# setting up color schemes
metric_colors <- c("#77AADD", "#EE8866", "#EEDD88", "#FFAABB")
stage_colors <- c("#44BB99", "#BBCC33", "#AAAA00", "#99DDFF")

####### Plotting relationship between RII and precision or sd_RII #######
ggplot(effects_df)+
  geom_point(aes(y = se_RII, x = RII), alpha = .5)+
  scale_y_reverse()+
  facet_wrap(~ metric_category)


####### Performing an Egger's hypothesis test in a Bayesian framework #######

# first need to calculate scaled effect sizes, then assess relationship between these and standard error of estimates
# test whether intercept of this relationship differs from 0


effects_scaled <- effects_df %>% 
  mutate(scaled_RII = RII/se_RII,
         precision_RII = 1/se_RII,
         samplesize_RII = 1/n_symbiotic) %>% 
  filter(scaled_RII> -Inf & scaled_RII < Inf) %>% 
  filter(precision_RII<30000) %>%
  mutate(metric_category_nice = case_when(grepl("growth", metric_category) ~ "Growth",
                                          grepl("reproduction", metric_category) ~ "Reproduction",
                                          grepl("recruitment", metric_category) ~ "Recruitment",
                                          grepl("survival", metric_category) ~ "Survival"))
# View(effects_scaled %>% filter(precision_RII>30000)) #In central analysis, I will drop experiment number 112-1


# looking at the relationship described in Egger's test

ggplot(effects_scaled)+
  geom_point(aes(x=scaled_RII, y = precision_RII))+
  facet_wrap(~ metric_category)

  

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

######## Bias test across vital rates #####
fit <- brm(formula = scaled_RII ~ 0 + metric_category + metric_category*precision_RII + (1|study_number),
           data = effects_scaled, 
           family = "gaussian",#"student",
           prior = c(#set_prior("gamma(1666, .1)", class = "nu"),
                     set_prior("normal(0,100)", class = "b"),
                     set_prior("normal(0,100)", class = "sigma")),
           iter = mcmc_pars$iter,
           chains = mcmc_pars$chains,
           warmup = mcmc_pars$warmup)
summary(fit)





# getting and plotting the model prediction
precision_RII_growth = seq(min((effects_scaled %>% filter(metric_category == "growth"))$precision_RII), max((effects_scaled %>% filter(metric_category == "growth"))$precision_RII), length = 100)
precision_RII_survival = seq(min((effects_scaled %>% filter(metric_category == "survival"))$precision_RII), max((effects_scaled %>% filter(metric_category == "survival"))$precision_RII), length = 100)
precision_RII_repro = seq(min((effects_scaled %>% filter(metric_category == "reproduction"))$precision_RII), max((effects_scaled %>% filter(metric_category == "reproduction"))$precision_RII), length = 100)
precision_RII_recruit = seq(min((effects_scaled %>% filter(metric_category == "recruitment"))$precision_RII), max((effects_scaled %>% filter(metric_category == "recruitment"))$precision_RII), length = 100)

prediction_growth_df <- tibble(metric_category = "growth",study_number = NA,precision_RII = precision_RII_growth)
prediction_survival_df <- tibble(metric_category = "survival",study_number = NA,precision_RII = precision_RII_survival)
prediction_repro_df <- tibble(metric_category = "reproduction",study_number = NA,precision_RII = precision_RII_repro)
prediction_recruit_df <- tibble(metric_category = "recruitment",study_number = NA,precision_RII = precision_RII_recruit)
prediction_df <- bind_rows(prediction_growth_df, prediction_survival_df, prediction_repro_df, prediction_recruit_df)


# plotting overall mean prediction
preds <- fitted(fit, newdata = prediction_df, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), re_formula = NA)
prediction_df <- bind_cols(prediction_df, preds) %>%
  mutate(metric_category_nice = case_when(grepl("growth", metric_category) ~ "Growth",
                                     grepl("reproduction", metric_category) ~ "Reproduction",
                                     grepl("recruitment", metric_category) ~ "Recruitment",
                                     grepl("survival", metric_category) ~ "Survival"))



Eggers_plot <- ggplot(prediction_df)+
  geom_point(data = effects_scaled, aes(x = precision_RII, y = scaled_RII))+
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, x = precision_RII), alpha = .2)+
  geom_ribbon(aes(ymin = Q25, ymax = Q75, x = precision_RII), alpha = .2, fill = "blue")+
  geom_line(aes(y = Estimate, x = precision_RII), color = "blue")+
  facet_wrap(~metric_category_nice, ncol = 1, scales = "free")+
  labs(y = expression(RII/SE[RII]), x = expression(1/SE[RII]))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "grey95"))
Eggers_plot

# plotting the posteriors


posts <- tidybayes::spread_draws(fit, b_metric_categorygrowth, b_metric_categoryreproduction, b_metric_categoryrecruitment, b_metric_categorysurvival, ndraws = 1000) %>% 
  pivot_longer(cols = b_metric_categorygrowth:b_metric_categorysurvival, names_to = "parameter", values_to = "value") %>% 
  mutate(metric_category = case_when(grepl("growth", parameter) ~ "Growth",
                                     grepl("reproduction", parameter) ~ "Reproduction",
                                     grepl("recruitment", parameter) ~ "Recruitment",
                                     grepl("survival", parameter) ~ "Survival"))

posts_summary <- posts %>%
  group_by(metric_category) %>% 
  summarize(n_draws = n(),
            mean = mean(value),
            prob_neg = sum(value<0)/n_draws*100,
            Q97.5 = quantile(value, .975),
            Q2.5 = quantile(value, .025))


posts_df <- posts %>% 
  mutate(CI = case_when(metric_category == "Growth" & value < posts_summary[posts_summary$metric_category=="Growth",]$Q97.5 & value > posts_summary[posts_summary$metric_category=="Growth",]$Q2.5 ~ "inside",
                        metric_category == "Reproduction" & value < posts_summary[posts_summary$metric_category=="Reproduction",]$Q97.5 & value > posts_summary[posts_summary$metric_category=="Reproduction",]$Q2.5 ~ "inside",
                        metric_category == "Recruitment" & value < posts_summary[posts_summary$metric_category=="Recruitment",]$Q97.5 & value > posts_summary[posts_summary$metric_category=="Recruitment",]$Q2.5 ~ "inside",
                        metric_category == "Survival" & value < posts_summary[posts_summary$metric_category=="Survival",]$Q97.5 & value > posts_summary[posts_summary$metric_category=="Survival",]$Q2.5 ~ "inside",
                        TRUE ~ "out"))

bias_post_plot <- ggplot(posts_df)+
  geom_vline(aes(xintercept = 0),  color = "black", lwd = .2)+
  geom_histogram(aes(x = value, fill = CI), position = "identity", bins = 50, alpha = .7)+
  facet_wrap(~metric_category, scales = "free", ncol = 1)+
  scale_fill_manual(values = c("black", "grey50"))+
  guides(fill = "none")+
  labs(x = "Intercept Posterior", y = "")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "grey95"))
bias_post_plot


bias_plot <- Eggers_plot + bias_post_plot + plot_annotation(tag_levels = "A")  
bias_plot
ggsave(bias_plot, filename = "Plots/bias_plot_metrics.png", width = 5, height = 5)




######## Bias test across life stages #####
fit <- brm(formula = scaled_RII ~ 0 + lifestage_general + lifestage_general*precision_RII + (1|study_number),
           data = effects_scaled, 
           family = "gaussian",#"student",
           prior = c(#set_prior("gamma(1666, .1)", class = "nu"),
             set_prior("normal(0,100)", class = "b"),
             set_prior("normal(0,100)", class = "sigma")),
           iter = mcmc_pars$iter,
           chains = mcmc_pars$chains,
           warmup = mcmc_pars$warmup)
summary(fit)





# getting and plotting the model prediction
precision_RII_embryo = seq(min((effects_scaled %>% filter(lifestage_general == "embryo"))$precision_RII), max((effects_scaled %>% filter(metric_category == "growth"))$precision_RII), length = 100)
precision_RII_juvenile = seq(min((effects_scaled %>% filter(lifestage_general == "juvenile"))$precision_RII), max((effects_scaled %>% filter(metric_category == "survival"))$precision_RII), length = 100)
precision_RII_adult = seq(min((effects_scaled %>% filter(lifestage_general == "adult"))$precision_RII), max((effects_scaled %>% filter(metric_category == "reproduction"))$precision_RII), length = 100)

prediction_embryo_df <- tibble(lifestage_general = "embryo",study_number = NA,precision_RII = precision_RII_growth)
prediction_juvenile_df <- tibble(lifestage_general = "juvenile",study_number = NA,precision_RII = precision_RII_survival)
prediction_adult_df <- tibble(lifestage_general = "adult",study_number = NA,precision_RII = precision_RII_repro)
prediction_df <- bind_rows(prediction_embryo_df, prediction_juvenile_df, prediction_adult_df)


# plotting overall mean prediction
preds <- fitted(fit, newdata = prediction_df, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), re_formula = NA)
prediction_df <- bind_cols(prediction_df, preds) %>%
  mutate(lifestage_category_nice = case_when(grepl("embryo", lifestage_general) ~ "Embryo",
                                          grepl("juvenile", lifestage_general) ~ "Juvenile",
                                          grepl("adult", lifestage_general) ~ "Adult"))



Eggers_plot <- ggplot(prediction_df)+
  geom_point(data = effects_scaled, aes(x = precision_RII, y = scaled_RII))+
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, x = precision_RII), alpha = .2)+
  geom_ribbon(aes(ymin = Q25, ymax = Q75, x = precision_RII), alpha = .2, fill = "blue")+
  geom_line(aes(y = Estimate, x = precision_RII), color = "blue")+
  facet_wrap(~lifestage_category_nice, ncol = 1, scales = "free")+
  labs(y = expression(RII/SE[RII]), x = expression(1/SE[RII]))+
  theme_bw()+
  theme(strip.background = element_rect(fill = "grey95"))
Eggers_plot

# plotting the posteriors


posts <- tidybayes::spread_draws(fit, b_lifestage_generalembryo, b_lifestage_generaljuvenile, b_lifestage_generaladult, ndraws = 1000) %>% 
  pivot_longer(cols = b_lifestage_generalembryo:b_lifestage_generaladult, names_to = "parameter", values_to = "value") %>% 
  mutate(lifestage_category = case_when(grepl("embryo", parameter) ~ "Embryo",
                                     grepl("juvenile", parameter) ~ "Juvenile",
                                     grepl("adult", parameter) ~ "Adult"))

posts_summary <- posts %>%
  group_by(lifestage_category) %>% 
  summarize(n_draws = n(),
            mean = mean(value),
            prob_neg = sum(value<0)/n_draws*100,
            Q97.5 = quantile(value, .975),
            Q2.5 = quantile(value, .025))


posts_df <- posts %>% 
  mutate(CI = case_when(lifestage_category == "Embryo" & value < posts_summary[posts_summary$lifestage_category=="Embryo",]$Q97.5 & value > posts_summary[posts_summary$lifestage_category=="Embryo",]$Q2.5 ~ "inside",
                        lifestage_category == "Juvenile" & value < posts_summary[posts_summary$lifestage_category=="Juvenile",]$Q97.5 & value > posts_summary[posts_summary$lifestage_category=="Juvenile",]$Q2.5 ~ "inside",
                        lifestage_category == "Adult" & value < posts_summary[posts_summary$lifestage_category=="Adult",]$Q97.5 & value > posts_summary[posts_summary$lifestage_category=="Adult",]$Q2.5 ~ "inside",
                        TRUE ~ "out"))

bias_post_plot <- ggplot(posts_df)+
  geom_vline(aes(xintercept = 0),  color = "black", lwd = .2)+
  geom_histogram(aes(x = value, fill = CI), position = "identity", bins = 50, alpha = .7)+
  facet_wrap(~lifestage_category, scales = "free", ncol = 1)+
  scale_fill_manual(values = c("black", "grey50"))+
  guides(fill = "none")+
  labs(x = "Intercept Posterior", y = "")+
  theme_bw()+
  theme(strip.background = element_rect(fill = "grey95"))
bias_post_plot


bias_plot <- Eggers_plot + bias_post_plot + plot_annotation(tag_levels = "A")  
bias_plot
ggsave(bias_plot, filename = "Plots/bias_plot_lifestage.png", width = 5, height = 5)






####### Visualizing the general categories across data and mapping ##########
effects_df_summary <- effects_scaled %>% 
  group_by(country) %>% 
  dplyr::summarize(article_count = length(unique(study_number)),
                   experiment_count = length(unique(experiment_label)))

world <- map_data("world") 
world <- subset(world, region != "Antarctica")

article_map <- ggplot(effects_df_summary)+
  geom_map(data = world, map = world, aes(map_id = region),
    fill = "white", color = "#7f7f7f", size = 0.25
  ) +
  geom_map(map = world, aes(map_id = country, fill = article_count), size = 0.25) +
  scale_fill_gradient(low = "#fff7bc", high = "#cc4c02", name = "# of Articles") +
  expand_limits(x = world$long, y = world$lat)+
  coord_sf()+
  labs(x = "Longitude", y = "Latitude")+
  theme_bw()

article_map


experiment_map <- ggplot(effects_df_summary)+
  geom_map(data = world, map = world, aes(map_id = region),
           fill = "white", color = "#7f7f7f", linewidth = 0.25
  ) +
  geom_map(map = world, aes(map_id = country, fill = experiment_count), size = 0.25) +
  scale_fill_gradient(low = "#fff7bc", high = "#cc4c02", name = "# of Experiments") +
  expand_limits(x = world$long, y = world$lat)+
  coord_sf()+
  labs(x = "Longitude", y = "Latitude")+
  theme_bw()


experiment_map


study_map <- article_map /experiment_map + plot_annotation(tag_levels = "A")
study_map

ggsave(study_map, filename = "Plots/study_map.png", width = 10, height = 8)




# making some histograms of taxonomy/VR coverage
ggplot(effects_df) + 
  geom_histogram(aes(x = host_order, fill = host_class), stat = "count")+
  labs(y = "# of measurements", x = "Order")+
  theme_bw()+
  theme(axis.text = element_text(hjust = 1, angle = 45))

effects_VR_count <- effects_scaled %>% 
  group_by(metric_category, host_order) %>% 
  summarize(count = n())

VR_count_plot <- ggplot(effects_VR_count)+
  geom_point(aes( x = factor(metric_category, levels = c("growth", "survival", "reproduction", "recruitment")), y = host_order, size = count, color = metric_category))+
  scale_color_manual(values = metric_colors)+
  labs(x = "", y = "Host Order", color = "Vital Rate")+
  theme_bw()+
  theme(axis.text = element_text(hjust = 1, angle = 45, size = rel(1)))
VR_count_plot


effects_LS_count <- effects_scaled %>% 
  group_by(lifestage_general, host_order) %>% 
  summarize(count = n())

LS_count_plot <- ggplot(effects_LS_count)+
  geom_point(aes( x= factor(lifestage_general, levels = c("embryo", "juvenile", "adult", "combines multiple")), y = host_order, size = count, color = lifestage_general))+
  scale_color_manual(values = c("embryo" = stage_colors[1],
                                "juvenile" = stage_colors[2],
                                "adult" = stage_colors[3],
                                "combines multiple" = stage_colors[4]))+
  labs(x = "", y = "", color = "Life Stage")+
  theme_bw()+
  theme(axis.text = element_text(hjust = 1, angle = 45,size = rel(1)),
        axis.text.y = element_blank())
LS_count_plot


counts_plot <- VR_count_plot + LS_count_plot + plot_annotation(tag_levels = "A")

# counts_plot
ggsave(counts_plot, filename = "Plots/counts_plot.png", width = 7, height = 6)


# summaries across the entire dataset
order_summary <- effects_scaled %>% 
  group_by(host_family) %>% 
  summarize(n())
taxonomy_summary <- effects_scaled %>% 
  group_by(.) %>% 
  dplyr::summarize(n_order = length(unique(host_order)),
                   n_family = length(unique(host_family)),
                   n_genus = length(unique(host_genus)))


paper_summary <- effects_scaled %>% 
  group_by(.) %>% 
  dplyr::summarize(n_articles = length(unique(study_number)),
                   n_experiment = length(unique(experiment_label)))

symbiota_counts <- effects_scaled %>% 
  mutate(symbiota_id = paste(host_species_clean, symbiont_species, sep = "_")) %>% 
  # summarize(length(unique(symbiota_id)))
  group_by(symbiota_id) %>% 
  summarize(n_symbiota = length(unique(metric_description)))


juv_summary <- effects_LS_count %>% 
  group_by(lifestage_general) %>%
  summarize(number_juv = sum(count),
            percent_juv = number_juv/sum(effects_LS_count$count))

growth_summary <- effects_VR_count %>% 
  group_by(metric_category) %>% 
  summarize(number = sum(count),
            percent = number/sum(effects_VR_count$count))


             
lab_summary <- effects_scaled %>% 
  group_by(experiment_setting) %>% 
  summarize(number = n(),
            percent = number/sum(effects_VR_count$count))



