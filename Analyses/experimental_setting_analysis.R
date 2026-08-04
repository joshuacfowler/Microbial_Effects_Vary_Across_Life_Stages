# Title: How do effects of lab and natural systems differ? #####
# Purpose: evaluates how effect sizes differ between lab and greenhouse environments and how experiments vs natural symbiont treatments compare to one another
# Authors: Josh Fowler  #####
# Date: Oct 13, 2025 #####

library(tidyverse)

library(lme4)
library(rstan)
library(brms)

library(patchwork)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }

# setting up color schemes
metric_colors <- c("#77AADD", "#FFAABB",  "#EEDD88", "#EE8866")
stage_colors <- c("#44BB99", "#BBCC33", "#AAAA00", "#99DDFF")


####### Reading in the data   #######
# The raw data is stored in Teams; we have downloaded the most recent version to a local directory as of May 6, 2025
# Data is imported, goes through some organization/cleaning, and we connect each study to its taxonomic information in the data_processing.R script

####### Reading in the data   #######
effects_df <- read_csv("effects_df.csv") %>% 
  filter(experiment_id != "112-1") %>% 
  filter(!is.na(sd_RII)) %>% 
  mutate(RII_01 = (RII+1)/2) %>% 
  mutate(exp_obs = case_when(grepl("natural", status_symbiotic) ~ "natural",
                              grepl("experimental", status_symbiotic) ~ "experimental"),
         setting = case_when(experiment_setting %in% c("greenhouse", "shadehouse", "glasshouse", "wire house") ~ "greenhouse",
                             experiment_setting %in% c("lab", "growth chamber") ~ "lab",
                              experiment_setting %in% c("field") ~ "field")
         ) 

# first just looking at the raw means in effect size and the variability in effect size across the categories
setting_df <- effects_df %>% 
  group_by(setting, metric_category) %>% 
  summarize(mean_RII = mean(RII, na.rm = T),
            sd_RII = sd(RII),
            CV_RII = sd_RII/mean_RII)
ggplot(setting_df)+
  geom_point(aes(x = setting, y =  mean_RII))+
  facet_wrap(~metric_category)
ggplot(setting_df)+
  geom_point(aes(x = setting, y =  CV_RII))+
  facet_wrap(~metric_category)


setting_df <- effects_df %>% 
  group_by(setting, lifestage_general) %>% 
  summarize(mean_RII = mean(RII, na.rm = T),
            sd_RII = sd(RII),
            CV_RII = sd_RII/mean_RII)
ggplot(setting_df)+
  geom_point(aes(x = setting, y =  mean_RII))+
  facet_wrap(~lifestage_general)
ggplot(setting_df)+
  geom_point(aes(x = setting, y =  CV_RII))+
  facet_wrap(~lifestage_general)


treatment_df <- effects_df %>% 
  group_by(exp_obs, metric_category) %>% 
  summarize(mean_RII = mean(RII, na.rm = T),
            sd_RII = sd(RII),
            CV_RII = sd_RII/mean_RII)
ggplot(treatment_df)+
  geom_point(aes(x = exp_obs, y =  mean_RII))+
  facet_wrap(~metric_category)
ggplot(treatment_df)+
  geom_point(aes(x = exp_obs, y =  CV_RII))+
  facet_wrap(~metric_category)


treatment_df <- effects_df %>% 
  group_by(exp_obs, lifestage_general) %>% 
  summarize(mean_RII = mean(RII, na.rm = T),
            sd_RII = sd(RII),
            CV_RII = sd_RII/mean_RII)
ggplot(treatment_df)+
  geom_point(aes(x = exp_obs, y =  mean_RII))+
  facet_wrap(~lifestage_general)
ggplot(treatment_df)+
  geom_point(aes(x = exp_obs, y =  CV_RII))+
  facet_wrap(~lifestage_general)


####### Fitting meta-analytic model  #######


# fitting brms
## run this code to optimize computer system settings for MCMC
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
mcmc_pars <- list(
  iter = 4000, 
  warmup = 2000, 
  thin = 1, 
  chains = 4
)

#re-try with different syntax for nested effects
VR_setting_fit <- brm(formula = RII|se(sd_RII, sigma = TRUE) ~ 0 + metric_category*setting*exp_obs + (1||study_number) + (1||experiment_label) + (1+metric_category||host_family) + (1+metric_category||host_genus),
                        #metric_category*setting + metric_category*exp_obs + setting*exp_obs + (1|study_number) + (1|experiment_label) + (1+metric_category|host_family) + (1+metric_category|host_genus),
           data = effects_df, 
           family = "gaussian",
           prior = c(set_prior("normal(0,.2)", class = "b"),
                     set_prior("normal(0,1)", class = "sd"),
                     set_prior("normal(0,1)", class = "sigma")),
           iter = mcmc_pars$iter,
           chains = mcmc_pars$chains,
           warmup = mcmc_pars$warmup, 
           control = list(adapt_delta = 0.99))
VR_ppc_hist_plot <- pp_check(VR_setting_fit, type = "stat_grouped", group = "exp_obs", stat = "mean", ndraws = 500)+
  labs(x = "RII", y = "Density")+
  theme_classic()
VR_ppc_hist_plot



# getting and plotting the model prediction
prediction_df <- expand.grid( metric_category = unique(effects_df$metric_category),
                              setting = unique(effects_df$setting),
                              exp_obs = unique(effects_df$exp_obs),
                              study_number = NA,
                              experiment_id = NA,
                              host_order =  NA,
                              sd_RII = 0)
# plotting overall mean prediction
preds <- fitted(VR_setting_fit, newdata = prediction_df, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), re_formula = NA)
prediction_df <- bind_cols(prediction_df, preds) #%>% 




VR_setting_effects_plot <- ggplot(data = prediction_df)+
  geom_hline(aes(yintercept = 0), color = "black", lwd = .1)+
  geom_jitter(data = effects_df, aes( x= setting, y = RII, fill = metric_category), shape = 21, color = "white", width = .1, alpha = .8)+
  # geom_linerange(aes(x = metric_category, ymin = Q25, ymax = Q75), lwd = 1.2) + 
  geom_linerange(aes(x = setting, ymin = Q2.5, ymax = Q97.5), color = "grey25") + 
  geom_point(aes(x = setting, y = Estimate), shape=21, fill = "grey25", color = "white", size = 1.5) +
  scale_fill_manual(values = metric_colors)+
  guides(fill = "none")+
  labs(x = "", color = "Vital Rate")+
  facet_grid(metric_category ~ exp_obs)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1.1)))
VR_setting_effects_plot

VR_setting_effects_plot <- ggplot(data = prediction_df)+
  geom_hline(aes(yintercept = 0), color = "black", lwd = .1)+
  geom_jitter(data = effects_df, aes( x= exp_obs, y = RII, fill = factor(metric_category, levels = c("growth", "survival", "reproduction", "recruitment"))), shape = 21, color = "white", width = .1, alpha = .8)+
  # geom_linerange(aes(x = metric_category, ymin = Q25, ymax = Q75), lwd = 1.2) + 
  geom_linerange(aes(x = exp_obs, ymin = Q2.5, ymax = Q97.5), color = "grey25") + 
  geom_point(aes(x = exp_obs, y = Estimate), shape=21, fill = "grey25", color = "white", size = 1.5) +
  scale_fill_manual(values = metric_colors)+
  guides(fill = "none")+
  labs(x = "", color = "Vital Rate")+
  facet_grid(factor(metric_category, levels = c("growth", "survival", "reproduction", "recruitment")) ~ setting)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1.1)))
VR_setting_effects_plot



ggsave(VR_setting_effects_plot, filename = "Plots/VR_setting_effects_plot.png", width = 7, height = 5)


# some posterior summaries
posts <- tidybayes::epred_draws(VR_setting_fit, newdata = prediction_df, re_formula = NA, ndraws = 500)

exp_obs_posts <- posts %>% 
  pivot_wider(id_cols = c(metric_category, setting, .draw), names_from=(exp_obs), values_from=(.epred)) %>% 
  mutate(exp_effect = experimental-natural)

summary_exp_obs <- exp_obs_posts %>% 
  group_by(setting, metric_category) %>% 
  summarize(n_draws = n(),
            mean = mean(exp_effect),
            prop_pos = (sum(exp_effect>0)/n_draws)*100)
posts_summary <- posts %>%
  group_by(setting, exp_obs, metric_category) %>% 
  summarize(n_draws = n(),
            mean = mean(.epred),
            prop_pos = (sum(.epred>0)/n_draws)*100)

VR_setting_post_plot <- ggplot(posts)+
  geom_histogram(aes(x = .epred, fill= exp_obs), alpha = .5)+
  scale_fill_manual(values = metric_colors)+
  facet_grid(metric_category ~ setting)+
  theme_bw()
  
VR_setting_post_plot

VR_setting_post_plot <- ggplot(exp_obs_posts)+
  geom_vline(xintercept = 0)+
  geom_histogram(aes(x = exp_effect, fill = metric_category),bins = 50, alpha = .6)+
  scale_fill_manual(values = metric_colors)+
  facet_grid(metric_category ~ setting)+
  labs(title = "Effect of Exp. vs. Obs. Symbiont Treatment", x = "Exp. - Obs.", y = "Count")+
  guides(fill = "none")+
  theme_bw()

VR_setting_post_plot

ggsave(VR_setting_post_plot, filename = "Plots/VR_setting_posterior_plot.png", width = 7, height = 5)




# fitting model for life stage differences
#re-try with different syntax for nested effects
LS_setting_fit <- brm(formula = RII|se(sd_RII, sigma = TRUE) ~ 0 + lifestage_general*setting*exp_obs + (1||study_number) + (1||experiment_label) + (1+lifestage_general||host_family) + (1+lifestage_general||host_genus),
                      #lifestage_general*setting + lifestage_general*exp_obs + setting*exp_obs + (1|study_number) + (1|experiment_label) + (1+lifestage_general|host_family) + (1+lifestage_general|host_genus),
                      data = effects_df, 
                      family = "gaussian",
                      prior = c(set_prior("normal(0,.2)", class = "b"),
                                set_prior("normal(0,1)", class = "sd"),
                                set_prior("normal(0,1)", class = "sigma")),
                      iter = mcmc_pars$iter,
                      chains = mcmc_pars$chains,
                      warmup = mcmc_pars$warmup, 
                      control = list(adapt_delta = 0.99))
LS_ppc_hist_plot <- pp_check(LS_setting_fit, type = "stat_grouped", group = "exp_obs", stat = "mean", ndraws = 500)+
  labs(x = "RII", y = "Density")+
  theme_classic()
LS_ppc_hist_plot



# getting and plotting the model prediction
prediction_df <- expand.grid( lifestage_general = unique(effects_df$lifestage_general),
                              setting = unique(effects_df$setting),
                              exp_obs = unique(effects_df$exp_obs),
                              study_number = NA,
                              experiment_id = NA,
                              host_order =  NA,
                              sd_RII = 0)
# plotting overall mean prediction
preds <- fitted(LS_setting_fit, newdata = prediction_df, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), re_formula = NA)
prediction_df <- bind_cols(prediction_df, preds) #%>% 




LS_setting_effects_plot <- ggplot(data = prediction_df)+
  geom_hline(aes(yintercept = 0), color = "black", lwd = .1)+
  geom_jitter(data = effects_df, aes( x= setting, y = RII, fill = lifestage_general), shape = 21, color = "white", width = .1, alpha = .8)+
  # geom_linerange(aes(x = lifestage_general, ymin = Q25, ymax = Q75), lwd = 1.2) + 
  geom_linerange(aes(x = setting, ymin = Q2.5, ymax = Q97.5), color = "grey25") + 
  geom_point(aes(x = setting, y = Estimate), shape=21, fill = "grey25", color = "white", size = 1.5) +
  scale_fill_manual(values = stage_colors)+
  guides(fill = "none")+
  labs(x = "", color = "Vital Rate")+
  facet_grid(lifestage_general ~ exp_obs)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1.1)))
LS_setting_effects_plot

LS_setting_effects_plot <- ggplot(data = prediction_df)+
  geom_hline(aes(yintercept = 0), color = "black", lwd = .1)+
  geom_jitter(data = effects_df, aes( x= exp_obs, y = RII, fill = factor(lifestage_general, levels = c("embryo", "juvenile", "adult"))), shape = 21, color = "white", width = .1, alpha = .8)+
  # geom_linerange(aes(x = lifestage_general, ymin = Q25, ymax = Q75), lwd = 1.2) + 
  geom_linerange(aes(x = exp_obs, ymin = Q2.5, ymax = Q97.5), color = "grey25") + 
  geom_point(aes(x = exp_obs, y = Estimate), shape=21, fill = "grey25", color = "white", size = 1.5) +
  scale_fill_manual(values = stage_colors)+
  guides(fill = "none")+
  labs(x = "", color = "Vital Rate")+
  facet_grid(factor(lifestage_general, levels = c("embryo", "juvenile", "adult")) ~ setting)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1.1)))
LS_setting_effects_plot



ggsave(LS_setting_effects_plot, filename = "Plots/LS_setting_effects_plot.png", width = 7, height = 5)


# some posterior summaries
posts <- tidybayes::epred_draws(LS_setting_fit, newdata = prediction_df, re_formula = NA, ndraws = 500)

exp_obs_posts <- posts %>% 
  pivot_wider(id_cols = c(lifestage_general, setting, .draw), names_from=(exp_obs), values_from=(.epred)) %>% 
  mutate(exp_effect = experimental-natural)

posts_summary <- posts %>%
  group_by(setting, exp_obs, lifestage_general) %>% 
  summarize(n_draws = n(),
            mean = mean(.epred),
            prop_pos = (sum(.epred>0)/n_draws)*100)

LS_setting_post_plot <- ggplot(posts)+
  geom_histogram(aes(x = .epred, fill= exp_obs), alpha = .5)+
  scale_fill_manual(values = stage_colors)+
  facet_grid(lifestage_general ~ setting)+
  theme_bw()

LS_setting_post_plot

LS_setting_post_plot <- ggplot(exp_obs_posts)+
  geom_vline(xintercept = 0)+
  geom_histogram(aes(x = exp_effect, fill = factor(lifestage_general, levels = c("embryo", "juvenile", "adult"))),bins = 50, alpha = .6)+
  scale_fill_manual(values = stage_colors)+
  facet_grid(factor(lifestage_general, levels = c("embryo", "juvenile", "adult")) ~ setting)+
  labs(title = "Effect of Exp. vs. Obs. Symbiont treatment", x = "Exp. - Obs.", y = "Count")+
  guides(fill = "none")+
  theme_bw()

LS_setting_post_plot

ggsave(LS_setting_post_plot, filename = "Plots/LS_setting_posterior_plot.png", width = 7, height = 5)

