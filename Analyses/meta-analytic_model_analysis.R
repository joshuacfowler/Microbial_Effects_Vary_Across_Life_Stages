# Title: How commonly do microbial effects vary across host life stages? #####
# Purpose: Reads in data from literature search and runs meta-analytic regression models
# Authors: Josh Fowler and Gwen Pohlmann #####
# Date: May 14, 2024 #####

# library(renv)
# renv::init()
# renv::snapshot()
# renv::deactivate()


library(tidyverse)

library(effsize) 

library(lme4)
library(rstan)
library(brms)

library(rotl)
library(metafor)
library(patchwork)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }
####### Reading in the data   #######
# The raw data is stored in Teams; we have downloaded the most recent version to a local directory as of May 6, 2025
# Data is imported, goes through some organization/cleaning, and we connect each study to its taxonomic information in the data_processing.R script

####### Reading in the data   #######
effects_df <- read_csv("effects_df.csv") %>% 
  filter(metric_category!="population metric") %>%   
  filter(!is.na(sd_RII)) %>% 
  mutate(RII_01 = (RII+1)/2)



length(unique(effects_df$study_number))
length(unique(effects_df$experiment_label))
dim(effects_df)

######### plotting prelim data ########
ggplot(effects_df) +
  geom_histogram(aes(x = RII))+facet_wrap(~metric_category, scales = "free")+expand_limits(x = c(-1,1))
ggsave("./20241012_metric_effects.png")

ggplot(effects_df) +
  geom_histogram(aes(x = lRR))+facet_wrap(~metric_category, scales = "free")

ggplot(effects_df) +
  geom_histogram(aes(x = cohensD))+facet_wrap(~metric_category, scales = "free")

ggplot(effects_df) +
  geom_histogram(aes(x = RII))+facet_wrap(~lifestage_general, scales = "free")+expand_limits(x = c(-1,1))
ggsave("./20241012_lifestage_effects.png")


familyXorder_layout <- ggplot(effects_df)+
  geom_tile(aes(x = host_order, y = host_family))+
  theme_minimal()+theme(axis.text.x = element_text(hjust = 1, angle = 45))
ggsave(familyXorder_layout, filename = "Plots/familyXorder_layout.png")

# prelim funnel plot
# funnel(effects_df$RII, effects_df$var_RII, yaxis = "vi")



####### Fitting meta-analytic model  #######


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





fit <- brm(formula = RII|se(sd_RII, sigma = TRUE) ~ 0 + metric_category + (1|study_number) + (1|experiment_label) + (1+metric_category|host_order) + (1+metric_category|host_family) + (1+metric_category|host_genus),
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



# look into projpred for model comparison


  
# graphical posterior predictive check
VR_ppc_overlay_plot <- pp_check(fit, type = "dens_overlay", ndraws = 100)+
  labs(x = "RII", y = "Density")+
  theme_classic()
VR_ppc_overlay_plot
ggsave(VR_ppc_overlay_plot, filename = "Plots/VR_ppc_overlay_plot.png", width = 4, height = 4)

VR_ppc_hist_plot <- pp_check(fit, type = "stat_grouped", group = "metric_category", stat = "mean", ndraws = 500)+
  labs(x = "RII", y = "Density")+
  theme_classic()
VR_ppc_hist_plot
ggsave(VR_ppc_hist_plot, filename = "Plots/VR_ppc_hist_plot.png", width = 4, height = 4)

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
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
overall_VR_effects_plot
ggsave(overall_VR_effects_plot, filename = "Plots/overall_VR_effects_plot.png", width = 7, height = 5)


# some posterior summaries
posts <- tidybayes::epred_draws(fit, newdata = prediction_df, re_formula = NA, ndraws = 500)

posts_summary <- posts %>%
  group_by(metric_category) %>% 
  summarize(n_draws = n(),
            mean = mean(.epred),
            prop_pos = (sum(.epred>0)/n_draws)*100)
  



# plotting prediction for specific orders
# what are the 6 most common orders?
top_orders <- names(sort(table(effects_df$host_order), decreasing = TRUE)[1:6])


prediction_df <- expand.grid( metric_category = unique(effects_df$metric_category),
                              study_number = "22",
                              experiment_id = "22-1",
                              host_order =  "Asparagales",
                              host_family = "Orchidaceae",
                              host_genus = c("Dendrobium","Vanda"),
                              #host_order =  top_orders,
                              sd_RII = 0)
preds <- fitted(fit, newdata = prediction_df, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), re_formula = ~ (1 + metric_category|host_order) + (1+metric_category|host_family) + (1+metric_category|host_genus) )

prediction_df <- bind_cols(prediction_df, preds) #%>% 
# mutate(across(Estimate:Q97.5,~unscale(.)))

effects_df_filtered <- effects_df %>% filter(!is.na(RII), host_genus == "Vanda"| host_genus == "Dendrobium")#host_order %in% top_orders)

order_VR_effects_plot <- ggplot(data = prediction_df)+
  geom_hline(aes(yintercept = 0), color = "black", lwd = .1)+
  geom_jitter(data = effects_df_filtered, aes( x= metric_category, y = RII, fill = metric_category), shape = 21, color = "white", width = .1, alpha = .8)+
  # geom_linerange(aes(x = metric_category, ymin = Q25, ymax = Q75), lwd = 1.2) + 
  geom_linerange(aes(x = metric_category, ymin = Q2.5, ymax = Q97.5)) + 
  geom_point(aes(x = metric_category, y = Estimate), shape = 21, fill = "black", color = "white", size = 1.5) +
  facet_wrap(~host_order+host_family + host_genus, nrow = 2)+
  scale_fill_manual(values = metric_colors)+
  guides(fill = "none")+
  labs(x = "")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey95"))
order_VR_effects_plot
ggsave(order_VR_effects_plot, filename = "Plots/order_VR_effects_plot.png", width = 5, height = 5)



# some posterior summaries
posts <- tidybayes::epred_draws(fit, newdata = prediction_df, re_formula = NA, ndraws = 500)

posts_summary <- posts %>%
  group_by(host_order, metric_category) %>% 
  summarize(n_draws = n(),
            mean = mean(.epred),
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

combo_VR_plot <- overall_VR_effects_plot + order_VR_effects_plot + plot_layout(widths = c(1,2)) + plot_annotation(tag_levels = "A")


combo_VR_plot
ggsave(combo_VR_plot, filename = "Plots/combo_VR_plot.png", width = 10, height = 5)








##### Fitting a model for variation in microbial effect across life stage #######
# Version that does  incorporate measurement error
# for now filtering out unknown lifestages, but this is just an oversite in data extraction
LS_effects_df <- effects_df %>% filter(!is.na(lifestage_general))
fit <- brm(formula = RII|se(sd_RII, sigma = TRUE) ~ 0 + lifestage_general + (1|study_number) + (1|experiment_label) + (1+lifestage_general|host_order) + (1+lifestage_general|host_family) + (1+lifestage_general|host_genus),
           #(1|study_number) + (1|experiment_label) + (1+metric_category|host_order) + (1+metric_category|host_family),
           data = LS_effects_df, 
           family = "gaussian",
           prior = c(set_prior("normal(0,.25)", class = "b"),
                     set_prior("normal(0,.25)", class = "sd"),
                     set_prior("normal(0,.25)", class = "sigma")),
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
ggsave(LS_ppc_hist_plot, filename = "Plots/LS_ppc_hist_plot.png", width = 4, height = 4)


# getting and plotting the model prediction
prediction_df <- expand.grid( lifestage_general = na.omit(unique(LS_effects_df$lifestage_general)),
                              study_number = NA,
                              experiment_id = NA,
                              host_order =  NA,
                              sd_RII = 0)

# plotting overall mean prediction
preds <- fitted(fit, newdata = prediction_df, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), re_formula = NA)

prediction_df <- bind_cols(prediction_df, preds) #%>% 
# mutate(across(Estimate:Q97.5,~unscale(.)))






overall_LS_effects_plot <- ggplot(data = prediction_df)+
  geom_hline(aes(yintercept = 0), color = "black", lwd = .1)+
  geom_jitter(data = LS_effects_df, aes( x= factor(lifestage_general, levels = c("embryo", "juvenile", "adult", "combines multiple")), y = RII, fill = lifestage_general), shape = 21, color = "white", width = .1, alpha = .8)+
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
overall_LS_effects_plot
ggsave(overall_VR_effects_plot, filename = "Plots/overall_LS_effects_plot.png", width = 7, height = 5)


# some posterior summaries
posts <- tidybayes::epred_draws(fit, newdata = prediction_df, re_formula = NA, ndraws = 500)

posts_summary <- posts %>%
  group_by(lifestage_general) %>% 
  summarize(n_draws = n(),
            mean = mean(.epred),
            prop_pos = (sum(.epred>0)/n_draws)*100)




# plotting prediction for specific orders
# what are the 6 most common orders?
top_orders <- names(sort(table(effects_df$host_order), decreasing = TRUE)[1:6])


prediction_df <- expand.grid( lifestage_general = unique(LS_effects_df$lifestage_general),
                              study_number = NA,
                              experiment_id = NA,
                              host_order =  "Asparagales",
                              host_family = "Orchidaceae",
                              host_genus = "Paphiopedilum",#top_orders,
                              sd_RII = 0)
preds <- fitted(fit, newdata = prediction_df, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), re_formula = ~ (metric_category|host_order))

prediction_df <- bind_cols(prediction_df, preds) #%>% 
# mutate(across(Estimate:Q97.5,~unscale(.)))

effects_df_filtered <- LS_effects_df %>% filter(!is.na(RII),host_genus == "Paphiopedilum")
                                                #host_order %in% top_orders)

order_LS_effects_plot <- ggplot(data = prediction_df)+
  geom_hline(aes(yintercept = 0), color = "black", lwd = .1)+
  geom_jitter(data = effects_df_filtered, aes( x= factor(lifestage_general, levels = c("embryo", "juvenile", "adult", "combines multiple")), y = RII, fill = lifestage_general), shape = 21, color = "white", width = .1, alpha = .8)+
  # geom_linerange(aes(x = metric_category, ymin = Q25, ymax = Q75), lwd = 1.2) + 
  geom_linerange(aes(x = lifestage_general, ymin = Q2.5, ymax = Q97.5)) + 
  geom_point(aes(x = lifestage_general, y = Estimate), shape = 21, fill = "black", color = "white", size = 1.5) +
  facet_wrap(~host_order, nrow = 2)+
  scale_fill_manual(values = c("embryo" = stage_colors[1],
                               "juvenile" = stage_colors[2],
                               "adult" = stage_colors[3],
                               "combines multiple" = stage_colors[4]))+  guides(fill = "none")+
  labs(x = "")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "grey95"))
order_LS_effects_plot
ggsave(order_LS_effects_plot, filename = "Plots/order_LS_effects_plot.png", width = 7, height = 5)



# some posterior summaries
posts <- tidybayes::epred_draws(fit, newdata = prediction_df, re_formula = NA, ndraws = 500)

posts_summary <- posts %>%
  group_by(host_order, lifestage_general) %>% 
  summarize(n_draws = n(),
            mean = mean(.epred),
            prop_pos = (sum(.epred>0)/n_draws)*100)






# making a combo plot for the paper

combo_LS_plot <- overall_LS_effects_plot + order_LS_effects_plot + plot_layout(widths = c(1,2)) + plot_annotation(tag_levels = "A")


combo_LS_plot
ggsave(combo_LS_plot, filename = "Plots/combo_LS_plot.png", width = 10, height = 5)






















#######################################################################################################################
############ Evaluating how microbial symbionts effect correlations between vital rates/ stages     ###################
#######################################################################################################################

metrics_per_study <- effects_df %>% 
  group_by(host_genus, symbiont_genus) %>% 
  summarize(length(unique(metric_category)))

vital_rate_estimates <- effects_df %>% 
  filter(metric_description != "fluorescence intensity (CFU/ml)") %>% 
  dplyr::select(-RII, -cohensD, -lRR, -var_RII, -sd_RII) %>% 
  group_by(metric_category) %>% 
  # standardizing the vital rate estimates
  mutate(avg_mean_symbiotic = mean(mean_symbiotic),
         sd_mean_symbiotic = sd(mean_symbiotic),
         avg_mean_aposymbiotic = mean(mean_aposymbiotic),
         sd_mean_aposymbiotic = sd(mean_aposymbiotic),
         std_mean_symbiotic = (mean_symbiotic - avg_mean_symbiotic)/sd_mean_symbiotic,
         std_mean_aposymbiotic = (mean_aposymbiotic - avg_mean_aposymbiotic)/sd_mean_aposymbiotic) %>% 
  group_by(metric_category, experiment_label) %>% 
  pivot_wider(names_from = c(metric_category), values_from = c(mean_symbiotic, mean_aposymbiotic, std_mean_symbiotic, std_mean_aposymbiotic,sd_symbiotic, sd_aposymbiotic, se_symbiotic, se_aposymbiotic, n_symbiotic, n_aposymbiotic)) 

stage_estimates <- effects_df %>% 
  filter(metric_description != "fluorescence intensity (CFU/ml)") %>% 
  dplyr::select(-RII, -cohensD, -lRR, -var_RII, -sd_RII) %>% 
  group_by(lifestage_general) %>% 
  # standardizing the vital rate estimates
  mutate(avg_mean_symbiotic = mean(mean_symbiotic),
         sd_mean_symbiotic = sd(mean_symbiotic),
         avg_mean_aposymbiotic = mean(mean_aposymbiotic),
         sd_mean_aposymbiotic = sd(mean_aposymbiotic),
         std_mean_symbiotic = (mean_symbiotic - avg_mean_symbiotic)/sd_mean_symbiotic,
         std_mean_aposymbiotic = (mean_aposymbiotic - avg_mean_aposymbiotic)/sd_mean_aposymbiotic) %>% 
pivot_wider(names_from = c(metric_category), values_from = c(mean_symbiotic, mean_aposymbiotic, std_mean_symbiotic, std_mean_aposymbiotic,sd_symbiotic, sd_aposymbiotic, se_symbiotic, se_aposymbiotic, n_symbiotic, n_aposymbiotic))




# trying to visualize the potential conflicting vital rate responses

ggplot(vital_rate_estimates)+
  geom_jitter(aes(x = metric_category, y = std_mean_symbiotic))+
  geom_jitter(aes(x = metric_category, y = std_mean_aposymbiotic), color = "red")


ggplot(filter(vital_rate_estimates, std_mean_symbiotic<2.5))+
  geom_jitter(aes(x = metric_category, y = std_mean_symbiotic))+
  geom_jitter(aes(x = metric_category, y = std_mean_aposymbiotic), color = "red")



ggplot(vital_rate_estimates)+
  geom_jitter(aes(x = metric_category, y = mean_symbiotic))+
  geom_jitter(aes(x = metric_category, y = mean_aposymbiotic), color = "red")




ggplot(stage_estimates)+
  geom_jitter(aes(x = lifestage_general, y = std_mean_symbiotic))+
  geom_jitter(aes(x = lifestage_general, y = std_mean_aposymbiotic), color = "red")

ggplot(filter(stage_estimates, std_mean_symbiotic<2.5))+
  geom_jitter(aes(x = lifestage_general, y = std_mean_symbiotic))+
  geom_jitter(aes(x = lifestage_general, y = std_mean_aposymbiotic), color = "red")


ggplot(stage_estimates)+
  geom_jitter(aes(x = lifestage_general, y = mean_symbiotic))+
  geom_jitter(aes(x = lifestage_general, y = mean_aposymbiotic), color = "red")

bf1 <- bf(std_mean_aposymbiotic_growth| mi() ~ (1|study_number) + (1|experiment_label))
bf2 <- bf(std_mean_aposymbiotic_survival| mi() ~ (1|study_number) + (1|experiment_label))
bf3 <- bf(std_mean_aposymbiotic_reproduction| mi() ~ (1|study_number) + (1|experiment_label))
bf4 <- bf(std_mean_aposymbiotic_recruitment| mi() ~ (1|study_number) + (1|experiment_label))


f2 <- 
  brm(data = stage_estimates, 
      family = gaussian,
      mvbf(bf1,bf2,bf3,bf4, rescor = TRUE),
      prior = c(#prior(gamma(2, .1), class = nu),
                # prior(normal(0, 10), class = Intercept, resp = x),
                # prior(normal(0, 10), class = Intercept, resp = y),
                # prior(normal(0, 10), class = sigma, resp = x),
                # prior(normal(0, 10), class = sigma, resp = y),
                prior(lkj(1), class = rescor)),
      iter = 2000, warmup = 500, chains = 4, cores = 4, 
      seed = 1234)

summary(f2)



# getting and plotting the model prediction
prediction_df <- expand.grid(
                              study_number = NA,
                              experiment_id = NA,
                              host_order =  NA,
                              sd_RII = 0)

# plotting overall mean prediction
preds <- fitted(f2, newdata = prediction_df, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), re_formula = NA)

prediction_df <-  data.frame(t(preds[1,,1]), metric = dimnames(preds)[[3]][1]) #%>% 
prediction_df <- bind_rows(prediction_df, data.frame(t(preds[1,,2]), metric = dimnames(preds)[[3]][2])) #%>% 
prediction_df <- bind_rows(prediction_df, data.frame(t(preds[1,,3]), metric = dimnames(preds)[[3]][3])) #%>% 
prediction_df <- bind_rows(prediction_df, data.frame(t(preds[1,,4]), metric = dimnames(preds)[[3]][4])) #%>% 




ggplot(filter(vital_rate_estimates, std_mean_symbiotic<2.5))+
  # geom_jitter(aes(x = metric_category, y = std_mean_symbiotic))+
  # geom_jitter(aes(x = metric_category, y = std_mean_aposymbiotic), color = "red")+
  geom_point(data = prediction_df, aes(x = metric, y = Estimate))+
  geom_linerange(data = prediction_df, aes(x = metric, ymin = Q2.5, ymax = Q97.5))



summary(f2)
var_corr <- VarCorr(f2)

correlations <- var_corr$residual__$cor

correlations1 <- data.frame(as_tibble(correlations[,,1], rownames = 'from'), "to" = dimnames(correlations)[[1]][[1]])
correlations2 <- data.frame(as_tibble(correlations[,,2], rownames = 'from'), "to" = dimnames(correlations)[[1]][[2]])
correlations3 <- data.frame(as_tibble(correlations[,,3], rownames = 'from'), "to" = dimnames(correlations)[[1]][[3]])
correlations4 <- data.frame(as_tibble(correlations[,,4], rownames = 'from'), "to" = dimnames(correlations)[[1]][[4]])

correlations_df <- bind_rows(correlations1, correlations2, correlations3, correlations4)

ggplot(correlations_df)+
  geom_tile(aes(x = from, y = to, fill = Estimate))+
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, limits = c(-1,1))
  # lims(fill = c(-1,1))


















##### old stuff

#











# Testing out simpler models
fit <- lm(data = effects_df, formula = RII ~ 0 + metric_category)
fit <- lmer(data = effects_df, formula = RII ~ 0 + metric_category + (1|study_number))
fit <- lmer(data = effects_df, formula = RII ~ 0 + metric_category + (1|study_number) + (1|study_number:experiment_id) )
fit <- lmer(data = effects_df, formula = RII ~ 0 + metric_category + (1|study_number) + (1|experiment_label) )
fit <- lmer(data = effects_df, formula = RII ~ 0 + metric_category + (1|study_number) + (1|experiment_label) + (1|host_order/host_family))

# fit <- lmer(data = effects_df, formula = RII ~ 0 + metric_category + (1|study_number) + (1|experiment_label) + (1|treatment_label) + (1|species_label) +(metric_category|study_number))

summary(fit)
# summary(fit)$r.squared



# Getting predictions from the model

prediction_df <- expand.grid( metric_category = unique(effects_df$metric_category),
                              study_number = NA,
                              host_order = c("Poales", "Fabales"),
                              host_family = NA)

preds <- predict(fit, newdat = prediction_df, se.fit = TRUE, re.form =  ~(1|host_order))

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



