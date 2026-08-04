# Title: How commonly do microbial effects vary across host life stages? #####
# Purpose: Reads in data from literature search and runs meta-analytic regression models
# Authors: Josh Fowler and Gwen Pohlmann #####
# Date: May 14, 2024 #####

# library(renv)
# renv::init()
# renv::snapshot()
# renv::deactivate()


library(tidyverse)


library(lme4)
#library(rstan)
library(cmdstanr)
library(brms)
library(bayesplot)

library(rotl)
library(patchwork)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }
####### Reading in the data   #######
# The raw data is stored in Teams; we have downloaded the most recent version to a local directory as of May 6, 2025
# Data is imported, goes through some organization/cleaning, and we connect each study to its taxonomic information in the data_processing.R script

####### Reading in the data   #######
effects_df <- read_csv("effects_df.csv") %>% 
  filter(experiment_id != "112-1") %>% 
  filter(!is.na(sd_RII)) %>% 
  mutate(RII = case_when(RII >= (1.0 - 1e-7) ~ 1.0,
                         RII <= (-1.0 + 1e-7) ~ -1.0, TRUE ~ RII),
         RII_01 = (RII+1)/2)




length(unique(effects_df$study_number))
length(unique(effects_df$experiment_label))
dim(effects_df)



######### plotting prelim data ########
# ggplot(effects_df) +
#   geom_histogram(aes(x = RII))+facet_wrap(~metric_category, scales = "free")+expand_limits(x = c(-1,1))
# ggsave("./20241012_metric_effects.png")
# 
# ggplot(effects_df) +
#   geom_histogram(aes(x = lRR))+facet_wrap(~metric_category, scales = "free")
# 
# ggplot(effects_df) +
#   geom_histogram(aes(x = cohensD))+facet_wrap(~metric_category, scales = "free")
# 
# ggplot(effects_df) +
#   geom_histogram(aes(x = RII))+facet_wrap(~lifestage_general, scales = "free")+expand_limits(x = c(-1,1))
# ggsave("./20241012_lifestage_effects.png")
# 
# 
# familyXorder_layout <- ggplot(effects_df)+
#   geom_tile(aes(x = host_order, y = host_family))+
#   theme_minimal()+theme(axis.text.x = element_text(hjust = 1, angle = 45))
# ggsave(familyXorder_layout, filename = "Plots/familyXorder_layout.png")

# prelim funnel plot
# funnel(effects_df$RII, effects_df$var_RII, yaxis = "vi")



####### Fitting meta-analytic model  #######


# fitting brms
## run this code to optimize computer system settings for MCMC
# rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())
# set.seed(123)

## MCMC settings
mcmc_pars <- list(
  iter = 4000, 
  warmup = 2000, 
  thin = 1, 
  chains = 4,
  seed = 123
)



##### Evaluating variation in RII across vital rate categories ###########
# Version that  incorporates measurement error
# zoib_model <- bf(
#   RII_01|se(sd_RII, sigma = TRUE) ~ 0 + metric_category,
#   phi ~ 1,
#   zoi ~ 1,
#   coi ~ 1,
#   family = zero_one_inflated_beta()
# )
# fit <- brm(formula = zoib_model,
#            #(1|study_number) + (1|experiment_label) + (1+metric_category|host_order) + (1+metric_category|host_family),
#            data = effects_df, 
#            family = zero_one_inflated_beta(),
#            prior = c(set_prior("normal(0,.2)", class = "b"),
#                      set_prior("normal(0,.25)", class = "sd"),
#                      set_prior("normal(0,.25)", class = "sigma")),
#            iter = mcmc_pars$iter,
#            chains = mcmc_pars$chains,
#            warmup = mcmc_pars$warmup, 
#            control = list(adapt_delta = 0.99))
# summary(fit)



# The following stan model encodes the equivalent model formula based on brms, but is sped up for truncation with some additional vectorization
# make a data list compatible with brms syntax
# mc_data_list <- make_standata(formula = bf(RII|trunc(lb = -1, ub = 1) +se(sd_RII, sigma = TRUE) ~ 0 + metric_category + (1|study_number) + (1|experiment_label) + (1+metric_category|host_family) + (1+metric_category|host_genus)), data = effects_df)
# mc_data_list <- make_standata(formula = bf(RII|trunc(lb = -1, ub = 1) +se(sd_RII, sigma = TRUE) ~ 0 + metric_category + (0+metric_category||study_number) + (0+metric_category||experiment_label) + (0+metric_category||host_family) + (0+metric_category||host_genus)), data = effects_df)
# mc_data_list$lb <- -1
# mc_data_list$ub <- 1
# mc_data_list$fitted <- 0
# mc_data_list$inc_re <- c(1,1,1,1)
# 
# # model <- cmdstan_model(stan_file = "Analyses/zoib.stan")
# model <- cmdstan_model(stan_file = "Analyses/ordered_beta.stan")
# 
# # model <- cmdstan_model(stan_file = "Analyses/trunc_01_normal.stan",
# #                        cpp_options = list(stan_threads = TRUE)) # this option is for within chain parallelization
# fit <- model$sample(data = mc_data_list,
#                     seed = mcmc_pars$seed,
#                     chains = mcmc_pars$chains,
#                     # parallel_chains = mcmc_pars$chains,
#                     # threads_per_chain = getOption("mc.cores", 1)/mcmc_pars$chains,
#                     iter_warmup = mcmc_pars$warmup,
#                     iter_sampling = mcmc_pars$iter,
#                     refresh = 10, 
#                     init = 0.1) # Forces initial values closer to zero/identity matrix)
# saveRDS(fit, file = "Analyses/zoib_fit.Rds")
# saveRDS(fit, file = "Analyses/zoib_fit_latent.Rds")

# fit <- stan(file = "Analyses/trunc_01_normal.stan", data = mc_data_list,
#                    iter = mcmc_pars$iter,
#                    warmup = mcmc_pars$warmup,
#                    chains = mcmc_pars$chains, 
#                    thin = mcmc_pars$thin)
fit <- brm(formula = bf(RII|se(sd_RII, sigma = TRUE) ~ 0 + metric_category + (1||study_number) + (1||experiment_label) + (1 + metric_category||host_family) + (1+metric_category||host_genus)),
             #(1|study_number) + (1|experiment_label) + (1+metric_category|host_order) + (1+metric_category|host_family),
           data = effects_df, 
           family = "gaussian",
           prior = c(set_prior("normal(0,.2)", class = "b"),
                     set_prior("normal(0,1)", class = "sd"),
                     set_prior("normal(0,1)", class = "sigma")),
           iter = mcmc_pars$iter,
           chains = mcmc_pars$chains,
           warmup = mcmc_pars$warmup, 
           control = list(adapt_delta = 0.99))

stan_model <- stancode(fit)

# summary(fit)



# fit2 <- brm(formula = RII|se(sd_RII, sigma = TRUE) ~ 0 + metric_category + (1|study_number) + (1|experiment_label) + (1+metric_category|host_order) + (1|experiment_label) + (1+metric_category|host_family) + (1+metric_category|host_genus),
#            #(1|study_number) + (1|experiment_label) + (1+metric_category|host_order) + (1+metric_category|host_family),
#            data = effects_df, 
#            family = "gaussian",
#            prior = c(set_prior("normal(0,.2)", class = "b"),
#                      set_prior("normal(0,1)", class = "sd"),
#                      set_prior("normal(0,1)", class = "sigma")),
#            iter = mcmc_pars$iter,
#            chains = mcmc_pars$chains,
#            warmup = mcmc_pars$warmup, 
#            control = list(adapt_delta = 0.99))

# loo(fit,fit2)
# loo_compare(fit,fit2)
# look into projpred for model comparison



# graphical posterior predictive check
VR_ppc_overlay_plot <- pp_check(fit, type = "dens_overlay", ndraws = 100)+
  labs(x = "RII", y = "Density")+
  theme_classic()
VR_ppc_overlay_plot
ggsave(VR_ppc_overlay_plot, filename = "Plots/VR_ppc_overlay_plot_just_family.png", width = 4, height = 4)

VR_ppc_hist_plot <- pp_check(fit, type = "stat_grouped", group = "metric_category", stat = "mean", ndraws = 500)+
  labs(x = "RII", y = "Density")+
  theme_classic()
VR_ppc_hist_plot
ggsave(VR_ppc_hist_plot, filename = "Plots/VR_ppc_hist_plot_just_family.png", width = 5, height = 4)


# getting and plotting the model prediction
prediction_df <- expand.grid( metric_category = unique(effects_df$metric_category),
                              study_number = NA,
                              experiment_id = NA,
                              host_order =  NA,
                              sd_RII = 0)

# plotting overall mean prediction
preds <- fitted(fit, newdata = prediction_df, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), re_formula = NA)
prediction_df <- bind_cols(prediction_df, preds) #%>% 



# setting up color schemes
metric_colors <- c("#77AADD", "#EE8866", "#EEDD88", "#FFAABB")
stage_colors <- c("#44BB99", "#BBCC33", "#AAAA00", "#99DDFF")

overall_VR_effects_plot <- ggplot(data = prediction_df)+
  geom_hline(aes(yintercept = 0), color = "black", lwd = .1)+
  geom_jitter(data = effects_df, aes( x= metric_category, y = RII, fill = metric_category), shape = 21, color = "white", width = .1, alpha = .8)+
  # geom_linerange(aes(x = metric_category, ymin = Q25, ymax = Q75), lwd = 1.2) + 
  geom_linerange(aes(x = metric_category, ymin = Q2.5, ymax = Q97.5), color = "grey25") + 
  geom_point(aes(x = metric_category, y = Estimate), shape=21, fill = "grey25", color = "white", size = 1.5) +
  scale_fill_manual(values = metric_colors)+
  guides(fill = "none")+
  labs(x = "", color = "Vital Rate")+
  # facet_wrap(~host_order)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1.1)))
overall_VR_effects_plot
ggsave(overall_VR_effects_plot, filename = "Plots/overall_VR_effects_plot_just_family.png", width = 7, height = 5)


# some posterior summaries
posts <- tidybayes::epred_draws(fit, newdata = prediction_df, re_formula = NA, ndraws = 500)

posts_summary <- posts %>%
  group_by(metric_category) %>% 
  summarize(n_draws = n(),
            mean = mean(.epred),
            median = median(.epred),
            Q2.5 = quantile(.epred, probs = c(0.025)),
            Q97.5 = quantile(.epred, probs = c(0.975)),
            prop_pos = (sum(.epred>0)/n_draws)*100)




# plotting prediction for specific orders
# what are the 6 most common orders?
top_families <- names(sort(table(effects_df$host_family), decreasing = TRUE)[1:6])


prediction_df <- expand.grid( metric_category = unique(effects_df$metric_category),
                              study_number = NA,
                              experiment_id = NA,
                              host_family = c("Poaceae","Orchidaceae"),
                              # host_genus = c("Dendrobium","Vanda"),
                              #host_order =  top_orders,
                              sd_RII = 0)
preds <- fitted(fit, newdata = prediction_df, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), re_formula = ~ (1+metric_category|host_family) )

prediction_df <- bind_cols(prediction_df, preds) #%>% 
# mutate(across(Estimate:Q97.5,~unscale(.)))

effects_df_filtered <- effects_df %>% filter(!is.na(RII), host_family %in% c("Poaceae","Orchidaceae" )) #host_genus == "Vanda"| host_genus == "Dendrobium")#host_order %in% top_orders)

family_VR_effects_plot <- ggplot(data = prediction_df)+
  geom_hline(aes(yintercept = 0), color = "black", lwd = .1)+
  geom_jitter(data = effects_df_filtered, aes( x= metric_category, y = RII, fill = metric_category), shape = 21, color = "white", width = .1, alpha = .8)+
  # geom_linerange(aes(x = metric_category, ymin = Q25, ymax = Q75), lwd = 1.2) + 
  geom_linerange(aes(x = metric_category, ymin = Q2.5, ymax = Q97.5)) + 
  geom_point(aes(x = metric_category, y = Estimate), shape = 21, fill = "black", color = "white", size = 1.5) +
  facet_wrap(~host_family, nrow = 2)+
  scale_fill_manual(values = metric_colors)+
  guides(fill = "none")+
  labs(x = "")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1.1)),
        strip.background = element_rect(fill = "grey95"))
family_VR_effects_plot
ggsave(family_VR_effects_plot, filename = "Plots/family_VR_effects_plot.png", width = 5, height = 5)



# some posterior summaries
posts <- tidybayes::epred_draws(fit, newdata = prediction_df, re_formula = ~ (1+metric_category|host_family), ndraws = 500)

posts_summary <- posts %>%
  group_by(host_family, metric_category) %>% 
  summarize(n_draws = n(),
            mean = mean(.epred),
            Q2.5 = quantile(.epred, probs = c(0.025)),
            Q97.5 = quantile(.epred, probs = c(0.975)),
            prop_pos = (sum(.epred>0)/n_draws)*100)



# plotting prediction for one specific genus
# prediction_df <- expand.grid( metric_category = unique(effects_df$metric_category),
#                               study_number = NA,
#                               experiment_id = NA,
#                               host_order =  "Asparagales",
#                               host_family = "Orchidaceae",
#                               host_genus = "Dendrobium",
#                               sd_RII = 0)
# preds <- fitted(fit, newdata = prediction_df, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), re_formula = ~ (metric_category|host_order))
# 
# prediction_df <- bind_cols(prediction_df, preds) #%>% 
# # mutate(across(Estimate:Q97.5,~unscale(.)))
# 
# effects_df_filtered <- effects_df %>% filter(!is.na(RII), host_genus == "Dendrobium")
# 
# genus_VR_effects_plot <- ggplot(data = prediction_df)+
#   geom_jitter(data = effects_df_filtered, aes( x= metric_category, y = RII, color = metric_category), width = .1, alpha = .2)+
#   # geom_linerange(aes(x = metric_category, ymin = Q25, ymax = Q75), lwd = 1.2) + 
#   geom_linerange(aes(x = metric_category, ymin = Q2.5, ymax = Q97.5)) + 
#   geom_point(aes(x = metric_category, y = Estimate), size = 1.5) +
#   facet_wrap(~host_genus)+
#   scale_color_manual(values = metric_colors)+
#   labs(x = "", color = "Vital Rate")+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(genus_VR_effects_plot, filename = "Plots/genus_VR_effects_plot.png", width = 7, height = 5)
# 









# making a combo plot for the paper

combo_VR_plot <- overall_VR_effects_plot + family_VR_effects_plot + plot_layout(widths = c(1.25,1)) + plot_annotation(tag_levels = "A")


combo_VR_plot
ggsave(combo_VR_plot, filename = "Plots/combo_VR_plot.png", width = 10, height = 5)





# Plotting prediction for all families

prediction_df <- expand.grid( metric_category = unique(effects_df$metric_category),
                              study_number = NA,
                              experiment_id = NA,
                              host_family = unique(effects_df$host_family),
                              # host_genus = c("Dendrobium","Vanda"),
                              #host_order =  top_orders,
                              sd_RII = 0)
preds <- fitted(fit, newdata = prediction_df, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), re_formula = ~ (1+metric_category|host_family) )

prediction_df <- bind_cols(prediction_df, preds) #%>% 
# mutate(across(Estimate:Q97.5,~unscale(.)))
families_rank <- names(sort(table(effects_df$host_family), decreasing = TRUE)[1:27])

all_family_VR_effects_plot <- ggplot(data = prediction_df)+
  geom_hline(aes(yintercept = 0), color = "black", lwd = .1)+
  geom_jitter(data = effects_df, aes( x= metric_category, y = RII, fill = metric_category), shape = 21, color = "white", width = .1, alpha = 1)+
  # geom_linerange(aes(x = metric_category, ymin = Q25, ymax = Q75), lwd = 1.2) + 
  geom_linerange(aes(x = metric_category, ymin = Q2.5, ymax = Q97.5),  alpha = .8) + 
  geom_point(aes(x = metric_category, y = Estimate), shape = 21, fill = "black", color = "white", size = 1.5, alpha = .8) +
  facet_wrap(~factor(host_family, levels = families_rank), nrow = 5)+
  scale_fill_manual(values = metric_colors)+
  guides(fill = "none")+
  labs(x = "")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1.1)),
        strip.background = element_rect(fill = "grey95"))
# all_family_VR_effects_plot
ggsave(all_family_VR_effects_plot, filename = "Plots/all_family_VR_effects_plot.png", width = 7, height = 9)

# some posterior summaries
posts <- tidybayes::epred_draws(fit, newdata = prediction_df, re_formula = ~ (1+metric_category|host_family) , ndraws = 500)

posts_summary <- posts %>%
  group_by(host_family, metric_category) %>% 
  summarize(n_draws = n(),
            mean = mean(.epred),
            Q2.5 = quantile(.epred, probs = c(0.025)),
            Q97.5 = quantile(.epred, probs = c(0.975)),
            prop_pos = (sum(.epred>0)/n_draws)*100)



##### Fitting a model for variation in microbial effect across life stage #######
# Version that does  incorporate measurement error
# for now filtering out unknown lifestages, but this is just an oversite in data extraction

fit <- brm(formula = bf(RII|se(sd_RII, sigma = TRUE) ~ 0 + lifestage_general + (1||study_number) + (1||experiment_label) + (1+lifestage_general||host_family) + (1+lifestage_general||host_genus)),
           #(1|study_number) + (1|experiment_label) + (1+metric_category|host_order) + (1+metric_category|host_family),
           data = effects_df, 
           family = "gaussian",
           prior = c(set_prior("normal(0,.2)", class = "b"),
                     set_prior("normal(0,1)", class = "sd"),
                     set_prior("normal(0,1)", class = "sigma")),
           iter = mcmc_pars$iter,
           chains = mcmc_pars$chains,
           warmup = mcmc_pars$warmup, 
           control = list(adapt_delta = 0.99))
summary(fit)




# graphical posterior predictive check
LS_ppc_overlay_plot <- pp_check(fit, type = "dens_overlay", ndraws = 100)+
  labs(x = "RII", y = "Density")+
  theme_classic()
LS_ppc_overlay_plot
ggsave(LS_ppc_overlay_plot, filename = "Plots/LS_ppc_overlay_plot.png", width = 4, height = 4)

LS_ppc_hist_plot <- pp_check(fit, type = "stat_grouped", group = "lifestage_general", stat = "mean", ndraws = 500)+
  labs(x = "RII", y = "Density")+
  theme_classic()
LS_ppc_hist_plot
ggsave(LS_ppc_hist_plot, filename = "Plots/LS_ppc_hist_plot.png", width = 8, height = 4)


# getting and plotting the model prediction
prediction_df <- expand.grid( lifestage_general = na.omit(unique(effects_df$lifestage_general)),
                              study_number = NA,
                              experiment_id = NA,
                              host_family =  NA,
                              sd_RII = 0)

# plotting overall mean prediction
preds <- fitted(fit, newdata = prediction_df, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), re_formula = NA)

prediction_df <- bind_cols(prediction_df, preds) #%>% 
# mutate(across(Estimate:Q97.5,~unscale(.)))


overall_LS_effects_plot <- ggplot(data = prediction_df)+
  geom_hline(aes(yintercept = 0), color = "black", lwd = .1)+
  geom_jitter(data = effects_df, aes( x= factor(lifestage_general, levels = c("embryo", "juvenile", "adult")), y = RII, fill = lifestage_general), shape = 21, color = "white", width = .1, alpha = .8)+
  # geom_linerange(aes(x = metric_category, ymin = Q25, ymax = Q75), lwd = 1.2) + 
  geom_linerange(aes(x = lifestage_general, ymin = Q2.5, ymax = Q97.5), color = "grey25") + 
  geom_point(aes(x = lifestage_general, y = Estimate), shape=21, fill = "grey25", color = "white", size = 1.5) +
  scale_fill_manual(values = c("embryo" = stage_colors[1],
                               "juvenile" = stage_colors[2],
                               "adult" = stage_colors[3],
                               "combines multiple" = stage_colors[4]))+
  guides(fill = "none")+
  labs(x = "", color = "Vital Rate")+
  # facet_wrap(~host_order)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,  size = rel(1.1)))
overall_LS_effects_plot
ggsave(overall_LS_effects_plot, filename = "Plots/overall_LS_effects_plot.png", width = 7, height = 5)


# some posterior summaries
posts <- tidybayes::epred_draws(fit, newdata = prediction_df, re_formula = NA, ndraws = 500)

posts_summary <- posts %>%
  group_by(lifestage_general) %>% 
  summarize(n_draws = n(),
            mean = mean(.epred),
            Q2.5 = quantile(.epred, probs = c(0.025)),
            Q97.5 = quantile(.epred, probs = c(0.975)),
            prop_pos = (sum(.epred>0)/n_draws)*100)


prob_juv_over_emb <-sum(posts[posts$lifestage_general == "juvenile",]$.epred > posts[posts$lifestage_general == "embryo",]$.epred)/500

prob_adult_over_emb <-sum(posts[posts$lifestage_general == "adult",]$.epred > posts[posts$lifestage_general == "embryo",]$.epred)/500

prob_juv_over_adult <-sum(posts[posts$lifestage_general == "juvenile",]$.epred > posts[posts$lifestage_general == "adult",]$.epred)/500


# plotting prediction for specific families
# what are the 6 most common families?

# plotting prediction for specific families
# what are the 6 most common families?
top_families <- names(sort(table(effects_df$host_family), decreasing = TRUE)[1:6])


prediction_df <- expand.grid( lifestage_general = unique(effects_df$lifestage_general),
                              study_number = NA,
                              experiment_id = NA,
                              host_family = c("Poaceae","Orchidaceae"),
                              # host_genus = c("Dendrobium","Vanda"),
                              #host_order =  top_orders,
                              sd_RII = 0)

preds <- fitted(fit, newdata = prediction_df, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), re_formula = ~ (1+ lifestage_general|host_family))

prediction_df <- bind_cols(prediction_df, preds) #%>% 
# mutate(across(Estimate:Q97.5,~unscale(.)))

effects_df_filtered <- effects_df %>% filter(!is.na(RII),host_family %in% c("Poaceae", "Orchidaceae"))
                                                #host_order %in% top_orders)

family_LS_effects_plot <- ggplot(data = prediction_df)+
  geom_hline(aes(yintercept = 0), color = "black", lwd = .1)+
  geom_jitter(data = effects_df_filtered, aes( x= factor(lifestage_general, levels = c("embryo", "juvenile", "adult")), y = RII, fill = lifestage_general), shape = 21, color = "white", width = .1, alpha = .8)+
  # geom_linerange(aes(x = metric_category, ymin = Q25, ymax = Q75), lwd = 1.2) + 
  geom_linerange(aes(x = lifestage_general, ymin = Q2.5, ymax = Q97.5)) + 
  geom_point(aes(x = lifestage_general, y = Estimate), shape = 21, fill = "black", color = "white", size = 1.5) +
  facet_wrap(~host_family, nrow = 2)+
  scale_fill_manual(values = c("embryo" = stage_colors[1],
                               "juvenile" = stage_colors[2],
                               "adult" = stage_colors[3]))+  guides(fill = "none")+
  labs(x = "")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1.2)),strip.background = element_rect(fill = "grey95"))
family_LS_effects_plot
ggsave(family_LS_effects_plot, filename = "Plots/family_LS_effects_plot.png", width = 7, height = 5)



# some posterior summaries
posts <- tidybayes::epred_draws(fit, newdata = prediction_df, re_formula = ~ (1+ lifestage_general|host_family), ndraws = 500)

posts_summary <- posts %>%
  group_by(host_family, lifestage_general) %>% 
  summarize(n_draws = n(),
            mean = mean(.epred),
            Q2.5 = quantile(.epred, probs = c(0.025)),
            Q97.5 = quantile(.epred, probs = c(0.975)),
            prop_pos = (sum(.epred>0)/n_draws)*100)


prob_juv_over_emb <-sum(posts[posts$lifestage_general == "juvenile" & posts$host_family == "Poaceae",]$.epred > posts[posts$lifestage_general == "embryo" & posts$host_family == "Poaceae",]$.epred)/500

prob_adult_over_emb <-sum(posts[posts$lifestage_general == "adult" & posts$host_family == "Poaceae",]$.epred > posts[posts$lifestage_general == "embryo" & posts$host_family == "Poaceae",]$.epred)/500

prob_juv_over_emb <-sum(posts[posts$lifestage_general == "juvenile" & posts$host_family == "Orchidaceae",]$.epred > posts[posts$lifestage_general == "embryo" & posts$host_family == "Orchidaceae",]$.epred)/500

prob_adult_over_emb <-sum(posts[posts$lifestage_general == "adult" & posts$host_family == "Orchidaceae",]$.epred > posts[posts$lifestage_general == "embryo" & posts$host_family == "Orchidaceae",]$.epred)/500

prob_juv_over_adult <-sum(posts[posts$lifestage_general == "juvenile" & posts$host_family == "Orchidaceae",]$.epred > posts[posts$lifestage_general == "adult" & posts$host_family == "Orchidaceae",]$.epred)/500


# making a combo plot for the paper

combo_LS_plot <- overall_LS_effects_plot + family_LS_effects_plot + plot_layout(widths = c(1.25,1)) + plot_annotation(tag_levels = "A")


# combo_LS_plot
ggsave(combo_LS_plot, filename = "Plots/combo_LS_plot.png", width = 10, height = 5)







# Plotting prediction for all families

prediction_df <- expand.grid( lifestage_general = unique(effects_df$lifestage_general),
                              study_number = NA,
                              experiment_id = NA,
                              host_family = unique(effects_df$host_family),
                              # host_genus = c("Dendrobium","Vanda"),
                              #host_order =  top_orders,
                              sd_RII = 0)
preds <- fitted(fit, newdata = prediction_df, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), re_formula = ~ (1+lifestage_general|host_family) )

prediction_df <- bind_cols(prediction_df, preds) #%>% 
# mutate(across(Estimate:Q97.5,~unscale(.)))
families_rank <- names(sort(table(effects_df$host_family), decreasing = TRUE)[1:27])

all_family_LS_effects_plot <- ggplot(data = prediction_df)+
  geom_hline(aes(yintercept = 0), color = "black", lwd = .1)+
  geom_jitter(data = effects_df, aes( x= factor(lifestage_general, levels = c("embryo", "juvenile", "adult")), y = RII, fill = lifestage_general), shape = 21, color = "white", width = .1, alpha = 1)+
  # geom_linerange(aes(x = metric_category, ymin = Q25, ymax = Q75), lwd = 1.2) + 
  geom_linerange(aes(x = lifestage_general, ymin = Q2.5, ymax = Q97.5),  alpha = .8) + 
  geom_point(aes(x = lifestage_general, y = Estimate), shape = 21, fill = "black", color = "white", size = 1.5, alpha = .8) +
  facet_wrap(~factor(host_family, levels = families_rank), nrow = 5)+
  scale_fill_manual(values = c("embryo" = stage_colors[1],
                               "juvenile" = stage_colors[2],
                               "adult" = stage_colors[3]))+  guides(fill = "none")+ 
  labs(x = "")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1.1)),
        strip.background = element_rect(fill = "grey95"))
# all_family_LS_effects_plot
ggsave(all_family_LS_effects_plot, filename = "Plots/all_family_LS_effects_plot.png", width = 8, height = 9)

# some posterior summaries
posts <- tidybayes::epred_draws(fit, newdata = prediction_df, re_formula = ~ (1+lifestage_general|host_family) , ndraws = 500)

posts_summary <- posts %>%
  group_by(host_family, lifestage_general) %>% 
  summarize(n_draws = n(),
            mean = mean(.epred),
            prop_pos = (sum(.epred>0)/n_draws)*100)



prob_juv_orchid_over_juvenile_everyone <-sum(posts[posts$lifestage_general == "juvenile" & posts$host_family == "Orchidaceae",]$.epred > posts[posts$lifestage_general == "juvenile" & posts$host_family != "Orchidaceae",]$.epred)/13000








