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
  filter(metric_category!="population metric") %>% 
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
         precision_RII = 1/se_RII) %>% 
  filter(scaled_RII<5000, scaled_RII> -Inf) %>% 
  mutate(metric_category_nice = case_when(grepl("growth", metric_category) ~ "Growth",
                                          grepl("reproduction", metric_category) ~ "Reproduction",
                                          grepl("recruitment", metric_category) ~ "Recruitment",
                                          grepl("survival", metric_category) ~ "Survival"))



# looking at the relationship described in Egger's test

ggplot(effects_scaled)+
  geom_point(aes(scaled_RII, x = precision_RII))  
  

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


fit <- brm(formula = scaled_RII ~ 0 + metric_category + metric_category*precision_RII + (1|study_number),
           data = effects_scaled, 
           family = "gaussian",
           prior = c(set_prior("normal(0,100)", class = "b"),
                     set_prior("normal(0,100)", class = "sigma")),
           iter = mcmc_pars$iter,
           chains = mcmc_pars$chains,
           warmup = mcmc_pars$warmup)
summary(fit)





# getting and plotting the model prediction
prediction_df <- expand.grid( metric_category = unique(effects_scaled$metric_category),
                              study_number = NA,
                              precision_RII = seq(min(effects_scaled$precision_RII), max(effects_scaled$precision_RII), length = 100))

# plotting overall mean prediction
preds <- fitted(fit, newdata = prediction_df, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), re_formula = NA)
prediction_df <- bind_cols(prediction_df, preds) %>%
  mutate(metric_category_nice = case_when(grepl("growth", metric_category) ~ "Growth",
                                     grepl("reproduction", metric_category) ~ "Reproduction",
                                     grepl("recruitment", metric_category) ~ "Recruitment",
                                     grepl("survival", metric_category) ~ "Survival"))



Eggers_plot <- ggplot(prediction_df)+
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, x = precision_RII), alpha = .2)+
  geom_ribbon(aes(ymin = Q25, ymax = Q75, x = precision_RII), alpha = .2)+
  geom_line(aes(y = Estimate, x = precision_RII))+
  geom_point(data = effects_scaled, aes(x = precision_RII, y = scaled_RII))+
  facet_wrap(~metric_category_nice, ncol = 1)+
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
ggsave(bias_plot, filename = "Plots/bias_plot.png", width = 5, height = 5)






####### Visualizing the general categories across data and mapping ##########
effects_df_summary <- effects_df %>% 
  group_by(country) %>% 
  dplyr::summarize(article_count = length(unique(study_number)),
                   experiment_count = length(unique(experiment_label)))

world <- map_data("world") 
world <- subset(world, region != "Antarctica")

article_map <- ggplot(effects_df_summary)+
  geom_map(data = world, map = world, aes(map_id = region),
    fill = "white", color = "#7f7f7f", linewidth = 0.25
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
