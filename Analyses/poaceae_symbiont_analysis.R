# Title: How do effects of Poaceae symbioses differ? #####
# Purpose: Performs analysis on subset of data for Poaceae hosts to evaluate how vital rate effects differ by symbiont taxa
# Authors: Josh Fowler  #####
# Date: Oct 13, 2025 #####

library(tidyverse)

library(lme4)
library(rstan)
library(brms)

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
  mutate(RII_01 = (RII+1)/2) %>% 
  filter(host_family == "Poaceae")

unique(effects_df$symbiont_genus_clean)


poa_df <- effects_df %>% 
  mutate(symbiont_type = case_when(symbiont_category == "Fungi" & !(symbiont_genus_clean %in% c("Epichloe", "Neotyphodium")) ~ "Mycorrhizal Fungi",
                                    symbiont_category == "Fungi" & symbiont_genus_clean %in% c("Epichloe", "Neotyphodium") ~ "Epichloe Fungi",
                                    symbiont_category == "Bacteria" ~ "Bacteria",
                                    TRUE ~ NA)) #%>% 
  # group_by(symbiont_type) %>% 
  # summarize(n())
symbiont_type_summary <- poa_df %>% 
  group_by(symbiont_type, symbiont_genus_clean) %>% 
  summarize(n = length(unique(host_species_clean)))


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

##### Fitting a model for variation in microbial effect across Vital Rate #######

fit <- brm(formula = RII|se(sd_RII, sigma = TRUE) ~ 0 + metric_category*symbiont_type + (1||study_number) + (1||experiment_label) + (1+metric_category||host_genus),
           #(1|study_number) + (1|experiment_label) + (1+metric_category|host_order) + (1+metric_category|host_family),
           data = poa_df, 
           family = "gaussian",
           prior = c(set_prior("normal(0,.2)", class = "b"),
                     set_prior("normal(0,1)", class = "sd"),
                     set_prior("normal(0,1)", class = "sigma")),
           iter = mcmc_pars$iter,
           chains = mcmc_pars$chains,
           warmup = mcmc_pars$warmup, 
           control = list(adapt_delta = 0.99))


# graphical posterior predictive check
VR_ppc_overlay_plot <- pp_check(fit, type = "dens_overlay", ndraws = 100)+
  labs(x = "RII", y = "Density")+
  theme_classic()
VR_ppc_overlay_plot
ggsave(VR_ppc_overlay_plot, filename = "Plots/Poaceae_VR_ppc_overlay_plot.png", width = 4, height = 4)

VR_ppc_hist_plot <- pp_check(fit, type = "stat_grouped", group = "symbiont_type", stat = "mean", ndraws = 500)+
  labs(x = "RII", y = "Density")+
  theme_classic()
VR_ppc_hist_plot
ggsave(VR_ppc_hist_plot, filename = "Plots/Poaceae_VR_ppc_hist_plot.png", width = 5, height = 4)



# getting and plotting the model prediction
prediction_df <- expand.grid( metric_category = unique(poa_df$metric_category),
                              symbiont_type = unique(poa_df$symbiont_type),
                              study_number = NA,
                              experiment_id = NA,
                              host_genus =  NA,
                              sd_RII = 0)

# plotting overall mean prediction
preds <- fitted(fit, newdata = prediction_df, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), re_formula = NA)
prediction_df <- bind_cols(prediction_df, preds) #%>% 



# setting up color schemes
metric_colors <- c("#77AADD", "#EE8866", "#EEDD88", "#FFAABB")
stage_colors <- c("#44BB99", "#BBCC33", "#AAAA00", "#99DDFF")
type_order <- c("Epichloe Fungi", "Mycorrhizal Fungi", "Bacteria")
type_labels <- c(`Epichloe Fungi` = expression("Epichloë Fungi test"),
                `Mycorrhizal Fungi` = "Mycorrhizal Fungi",
                `Bacteria` = "Bacteriatest")

VR_effects_plot <- ggplot(data = prediction_df)+
  geom_hline(aes(yintercept = 0), color = "black", lwd = .1)+
  geom_jitter(data = poa_df, aes( x= metric_category, y = RII, fill = metric_category), shape = 21, color = "white", width = .1, alpha = .8)+
  # geom_linerange(aes(x = metric_category, ymin = Q25, ymax = Q75), lwd = 1.2) + 
  geom_linerange(aes(x = metric_category, ymin = Q2.5, ymax = Q97.5), color = "grey25") + 
  geom_point(aes(x = metric_category, y = Estimate), shape=21, fill = "grey25", color = "white", size = 1.5) +
  scale_fill_manual(values = metric_colors)+
  guides(fill = "none")+
  labs(x = "", color = "Vital Rate")+
  facet_wrap(~factor(symbiont_type, levels = type_order, labels = c("italic(Epichloë)~Fungi", "Mycorrhizal~Fungi", "Bacteria")), labeller = label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1.1)))
VR_effects_plot
ggsave(VR_effects_plot, filename = "Plots/Poaceae_VR_effects_plot.png", width = 7, height = 5)



# some posterior summaries
posts <- tidybayes::epred_draws(fit, newdata = prediction_df, re_formula = NA, ndraws = 500)

posts_summary <- posts %>%
  group_by(symbiont_type, metric_category) %>% 
  summarize(n_draws = n(),
            mean = mean(.epred),
            Q2.5 = quantile(.epred, probs = c(0.025)),
            Q97.5 = quantile(.epred, probs = c(0.975)),
            prop_pos = (sum(.epred>0)/n_draws)*100)




##### Fitting a model for variation in microbial effect across life stage #######

fit_ls <- brm(formula = RII|se(sd_RII, sigma = TRUE) ~ 0 + lifestage_general*symbiont_type + (1||study_number) + (1||experiment_label) + (1+metric_category||host_genus),
           #(1|study_number) + (1|experiment_label) + (1+metric_category|host_order) + (1+metric_category|host_family),
           data = poa_df, 
           family = "gaussian",
           prior = c(set_prior("normal(0,.2)", class = "b"),
                     set_prior("normal(0,1)", class = "sd"),
                     set_prior("normal(0,1)", class = "sigma")),
           iter = mcmc_pars$iter,
           chains = mcmc_pars$chains,
           warmup = mcmc_pars$warmup, 
           control = list(adapt_delta = 0.99))


# graphical posterior predictive check
LS_ppc_overlay_plot <- pp_check(fit_ls, type = "dens_overlay", ndraws = 100)+
  labs(x = "RII", y = "Density")+
  theme_classic()
LS_ppc_overlay_plot
ggsave(LS_ppc_overlay_plot, filename = "Plots/Poaceae_LS_ppc_overlay_plot.png", width = 4, height = 4)

LS_ppc_hist_plot <- pp_check(fit_ls, type = "stat_grouped", group = "symbiont_type", stat = "mean", ndraws = 500)+
  labs(x = "RII", y = "Density")+
  theme_classic()
LS_ppc_hist_plot
ggsave(LS_ppc_hist_plot, filename = "Plots/Poaceae_LS_ppc_hist_plot.png", width = 5, height = 4)



# getting and plotting the model prediction
prediction_df <- expand.grid( lifestage_general = unique(poa_df$lifestage_general),
                              symbiont_type = unique(poa_df$symbiont_type),
                              study_number = NA,
                              experiment_id = NA,
                              host_genus =  NA,
                              sd_RII = 0)

# plotting overall mean prediction
preds <- fitted(fit_ls, newdata = prediction_df, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), re_formula = NA)
prediction_df <- bind_cols(prediction_df, preds) #%>% 



# setting up color schemes
metric_colors <- c("#77AADD", "#EE8866", "#EEDD88", "#FFAABB")
stage_colors <- c("#44BB99", "#BBCC33", "#AAAA00", "#99DDFF")

LS_effects_plot <- ggplot(data = prediction_df)+
  geom_hline(aes(yintercept = 0), color = "black", lwd = .1)+
  geom_jitter(data = poa_df, aes( x= factor(lifestage_general, levels = c("embryo", "juvenile", "adult")), y = RII, fill = lifestage_general), shape = 21, color = "white", width = .1, alpha = .8)+
  # geom_linerange(aes(x = lifestage_general, ymin = Q25, ymax = Q75), lwd = 1.2) + 
  geom_linerange(aes(x = lifestage_general, ymin = Q2.5, ymax = Q97.5), color = "grey25") + 
  geom_point(aes(x = lifestage_general, y = Estimate), shape=21, fill = "grey25", color = "white", size = 1.5) +
  scale_fill_manual(values = c("embryo" = stage_colors[1],
                               "juvenile" = stage_colors[2],
                               "adult" = stage_colors[3])) +
                      guides(fill = "none")+
  labs(x = "", color = "Life Stage")+
  facet_wrap(~factor(symbiont_type, levels = type_order, labels = c("italic(Epichloë)~Fungi", "Mycorrhizal~Fungi", "Bacteria")), labeller = label_parsed)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(1.1)))
LS_effects_plot
ggsave(LS_effects_plot, filename = "Plots/Poaceae_LS_effects_plot.png", width = 7, height = 5)



# some posterior summaries
posts <- tidybayes::epred_draws(fit_ls, newdata = prediction_df, re_formula = NA, ndraws = 500)

posts_summary <- posts %>%
  group_by(symbiont_type, lifestage_general) %>% 
  summarize(n_draws = n(),
            mean = mean(.epred),
            Q2.5 = quantile(.epred, probs = c(0.025)),
            Q97.5 = quantile(.epred, probs = c(0.975)),
            prop_pos = (sum(.epred>0)/n_draws)*100)

sum(posts[posts$lifestage_general == "embryo"& posts$symbiont_type == "Epichloe Fungi",]$.epred>posts[posts$lifestage_general %in% c( "adult") & posts$symbiont_type == "Epichloe Fungi",]$.epred)/500
sum(posts[posts$lifestage_general == "embryo"& posts$symbiont_type == "Epichloe Fungi",]$.epred>posts[posts$lifestage_general %in% c( "juvenile") & posts$symbiont_type == "Epichloe Fungi",]$.epred)/500
sum(posts[posts$lifestage_general == "embryo"& posts$symbiont_type == "Epichloe Fungi",]$.epred>posts[posts$lifestage_general %in% c( "adult","juvenile") & posts$symbiont_type == "Epichloe Fungi",]$.epred)/1000

sum(posts[posts$lifestage_general == "adult"& posts$symbiont_type == "Mycorrhizal Fungi",]$.epred>posts[posts$lifestage_general %in% c( "embryo") & posts$symbiont_type == "Mycorrhizal Fungi",]$.epred)/500
sum(posts[posts$lifestage_general == "adult"& posts$symbiont_type == "Mycorrhizal Fungi",]$.epred>posts[posts$lifestage_general %in% c( "juvenile") & posts$symbiont_type == "Mycorrhizal Fungi",]$.epred)/500
sum(posts[posts$lifestage_general == "adult"& posts$symbiont_type == "Mycorrhizal Fungi",]$.epred>posts[posts$lifestage_general %in% c( "embryo","juvenile") & posts$symbiont_type == "Mycorrhizal Fungi",]$.epred)/1000

sum(posts[posts$lifestage_general == "adult"& posts$symbiont_type == "Bacteria",]$.epred>posts[posts$lifestage_general %in% c( "embryo") & posts$symbiont_type == "Bacteria",]$.epred)/500
sum(posts[posts$lifestage_general == "adult"& posts$symbiont_type == "Bacteria",]$.epred>posts[posts$lifestage_general %in% c( "juvenile") & posts$symbiont_type == "Bacteria",]$.epred)/500
sum(posts[posts$lifestage_general == "adult"& posts$symbiont_type == "Bacteria",]$.epred>posts[posts$lifestage_general %in% c( "embryo","juvenile") & posts$symbiont_type == "Bacteria",]$.epred)/1000



poaceae_symbiont_plot <- VR_effects_plot + LS_effects_plot + plot_annotation(tag_levels = "A") + plot_layout(nrow = 2)
poaceae_symbiont_plot
ggsave(poaceae_symbiont_plot, filename = "Plots/poaceae_symbiont_plot.png", width=6, height = 8)










# Now performing permutation test to compare number of opposing vital rates








####### calculating the proportion of opposing symbiotic effects within a particular study  #######
# first we need to winnow down the data set to only those that have measurements of multiple metrics
# also note that I need to investigate the cases with RII exactly = 0
multiple_metrics_df <- poa_df %>%
  filter(experiment_id != "112-1") %>% 
  filter(!is.na(sd_RII)) %>% 
  filter(RII != 0) %>% 
  mutate(positive = case_when(RII>0 ~ "positive", 
                              RII <0 ~ "negative")) %>% 
  mutate(metricXstage = paste0(metric_description, lifestage_description, sep = "_")) %>% 
  group_by(host_family, host_species_clean,  symbiont_species) %>% 
  dplyr::summarise(n_metrics = length(unique(metric_description)),
                   n_metricsXstage = length(unique(metricXstage)),
                   n_rows = n(),
                   n_pos = sum(positive == "positive"),
                   n_negative = sum(positive == "negative"),
                   prop_pos = n_pos/n_rows,
                   opposing = case_when(prop_pos == 1 ~ "all positive",
                                        prop_pos == 0 ~ "all negative",
                                        prop_pos >0 & prop_pos <1 ~ "opposing")) %>% 
  filter(n_metricsXstage >1) %>% 
  mutate(symbiota_id = paste0(host_species_clean, symbiont_species, sep = "_"))



trt_order <- c("all negative", "opposing", "all positive")


permute_df <- poa_df %>%
  filter(RII != 0) %>% 
  mutate(metricXstage = paste0(metric_description, lifestage_description, sep = "_")) %>% 
  mutate(positive = case_when(RII>0 ~ "positive", 
                              RII <0 ~ "negative")) %>% 
  mutate(symbiota_id = paste0(host_species_clean, symbiont_species, sep = "_")) %>% 
  filter(symbiota_id %in% unique(multiple_metrics_df$symbiota_id))  %>% 
  dplyr::select(host_family, host_species_clean, symbiont_species,metric_category, lifestage_general, metricXstage, positive, symbiont_type, symbiota_id)


# trt_order <- c("all negative", "opposing", "all positive")
# ggplot(multiple_metrics_df)+
#   geom_histogram(aes(x = factor(opposing, level = trt_order), fill = opposing), stat = "count")+
#   theme_classic()



# repeating the permutation, but analyziing by symbiont_type
symbiont_permute_df <-  permute_df %>% 
  mutate(symbiota = paste(host_species_clean, symbiont_species, sep = "_")) %>% 
  group_by(symbiont_type) %>% 
  filter(length(unique(symbiota))>1) %>% 
  rename(real=positive) %>% 
  select(host_family, host_species_clean, symbiont_species, symbiont_type, metric_category, metricXstage, symbiota, real) # re-ordering to put symbiota in front of "real"

symbiont_rank <- c("Bacteria", "Epichloe", "Other Fungi")
permuted_symbiont <- list()
for(i in 1:length(symbiont_rank)){
  filtered_symbiont <- symbiont_permute_df[symbiont_permute_df$symbiont_type == symbiont_rank[i],]
  permuted_symbiont[[i]] <- as_tibble(cbind(filtered_symbiont, replicate(500, sample(filtered_symbiont$real))))  %>% 
    pivot_longer(cols = real:`...508`, names_to = "permutation", values_to = "effect")   %>% 
    group_by(host_family, host_species_clean, symbiont_species,  symbiont_type, permutation) %>% 
    summarize(n_metricsXstage = length(unique(metricXstage)),
              n_rows = n(),
              n_pos = sum(effect == "positive"),
              n_negative = sum(effect == "negative"),
              prop_pos = n_pos/n_rows,
              opposing = case_when(prop_pos == 1 ~ "all positive",
                                   prop_pos == 0 ~ "all negative",
                                   prop_pos >0 & prop_pos <1 ~ "opposing")) %>%  
    filter(n_metricsXstage >1) %>% 
    ungroup() %>% group_by( symbiont_type, permutation) %>% 
    summarize(`All Positive` =  sum(opposing == "all positive"),
              `All Negative` = sum(opposing == "all negative"),
              `Opposing` = sum(opposing == "opposing"))
  
}

permuted_symbionts <- bind_rows(permuted_symbiont)



permuted_symbiont_long <- permuted_symbionts %>% 
  # filter(permutation != "real") %>% 
  pivot_longer(cols = c(`All Positive`, `All Negative`, `Opposing`)) %>% 
  mutate(real_or = case_when(permutation == "real" ~ "Obs. data",
                             permutation != "real" ~ "Permuted data"))

permuted_symbiont_counts <- permuted_symbiont_long %>% filter(permutation != "real")
real_symbiont_counts <- permuted_symbiont_long %>% filter(permutation == "real")



trt_order <- c("All Negative", "Opposing", "All Positive")


# calculating quantiles from permuted data

permuted_symbiont_quant <- permuted_symbiont_long %>% 
  filter(real_or == "Permuted data") %>% 
  group_by(symbiont_type, name, real_or) %>% 
  summarize(mean = mean(value),
            Q02.5 = quantile(value, .025),
            Q97.5 = quantile(value, .975),
            n = n())

permutation_symbiont_plot <- ggplot(permuted_symbiont_counts)+
  geom_jitter(aes(y = factor(name, levels = trt_order), x = value, color= real_or, fill = real_or, alpha = real_or, shape = real_or, size = real_or), height = .3, width = 0.05)+
  geom_linerange(data = permuted_symbiont_quant, aes( y = factor(name, levels = trt_order), xmin = Q02.5, xmax = Q97.5), lwd = .3)+
  geom_point(data = permuted_symbiont_quant, aes(y = factor(name, levels = trt_order), x = mean), shape = 18)+
  geom_point(data = real_symbiont_counts, aes(y = factor(name, levels = trt_order), x = value, color= real_or, fill = real_or, alpha = real_or, shape = real_or, size = real_or))+
  scale_color_manual(values = c("red", "black"))+
  scale_fill_manual(values = c("red", NA))+
  scale_alpha_manual(values = c(1,.1))+
  scale_shape_manual(values = c(16,21))+
  scale_size_manual(values = c(2,1))+
  facet_wrap(~factor(symbiont_type, levels = symbiont_rank), scales = "free", nrow = 4)+
  scale_x_continuous(breaks = scales::breaks_pretty())+
  labs(x = "# of symbiota", y = "", color = "", fill = "", alpha = "", shape = "", size = "")+
  theme_bw()
permutation_symbiont_plot

ggsave(permutation_symbiont_plot, filename = "Plots/permutation_symbiont_plot.png", height = 5, width = 5)


######## repeating this within types of metrics, to see if this is driven by multiple metrics of growth rather than across survival/growth etc. ######
symbiont_permute_df <-  permute_df %>% 
  mutate(symbiota = paste(host_species_clean, symbiont_species, sep = "_")) %>% 
  group_by(symbiont_type) %>% 
  filter(length(unique(symbiota))>1) %>% 
  rename(real=positive) %>% 
  select(host_family, host_species_clean, symbiont_species, symbiont_type, metric_category, metricXstage, symbiota, real) # re-ordering to put symbiota in front of "real"


symbiont_rank <- c("Bacteria", "Epichloe", "Other Fungi")
metric_rank <- unique(poa_df$metric_category)
symbiont_iter <- expand.grid(symbiont_rank = symbiont_rank,metric_rank = metric_rank)
permuted_symbiont <- list()
for(i in 1:nrow(symbiont_iter)){
  filtered_symbiont <- symbiont_permute_df[symbiont_permute_df$symbiont_type == symbiont_iter[i,"symbiont_rank"] & symbiont_permute_df$metric_category == symbiont_iter[i,"metric_rank"],]
  permuted_symbiont[[i]] <- as_tibble(cbind(filtered_symbiont, replicate(500, sample(filtered_symbiont$real))))  %>% 
    pivot_longer(cols = real:`...508`, names_to = "permutation", values_to = "effect")  %>% 
    group_by(host_family, host_species_clean, symbiont_species, metric_category, symbiont_type, permutation) %>% 
    summarize(n_metricsXstage = length(unique(metricXstage)),
              n_rows = n(),
              n_pos = sum(effect == "positive"),
              n_negative = sum(effect == "negative"),
              prop_pos = n_pos/n_rows,
              opposing = case_when(prop_pos == 1 ~ "all positive",
                                   prop_pos == 0 ~ "all negative",
                                   prop_pos >0 & prop_pos <1 ~ "opposing")) %>%  
    filter(n_metricsXstage >1) %>% 
    ungroup() %>% group_by( symbiont_type, metric_category, permutation) %>% 
    summarize(`All Positive` =  sum(opposing == "all positive"),
              `All Negative` = sum(opposing == "all negative"),
              `Opposing` = sum(opposing == "opposing"))

}

permuted_symbionts <- bind_rows(permuted_symbiont)



permuted_symbiont_long <- permuted_symbionts %>% 
  # filter(permutation != "real") %>% 
  pivot_longer(cols = c(`All Positive`, `All Negative`, `Opposing`)) %>% 
  mutate(real_or = case_when(permutation == "real" ~ "Obs. data",
                             permutation != "real" ~ "Permuted data"))

permuted_symbiont_counts <- permuted_symbiont_long %>% filter(permutation != "real")
real_symbiont_counts <- permuted_symbiont_long %>% filter(permutation == "real")



trt_order <- c("All Negative", "Opposing", "All Positive")


# calculating quantiles from permuted data

permuted_symbiont_quant <- permuted_symbiont_long %>% 
  filter(real_or == "Permuted data") %>% 
  group_by(symbiont_type, metric_category, name, real_or) %>% 
  summarize(mean = mean(value),
            Q02.5 = quantile(value, .025),
            Q97.5 = quantile(value, .975),
            n = n())

permutation_metricXsymbiont_plot <- ggplot(permuted_symbiont_counts)+
  geom_jitter(aes(y = factor(name, levels = trt_order), x = value, color= real_or, fill = real_or, alpha = real_or, shape = real_or, size = real_or), height = .3, width = 0.05)+
  geom_linerange(data = permuted_symbiont_quant, aes( y = factor(name, levels = trt_order), xmin = Q02.5, xmax = Q97.5), lwd = .3)+
  geom_point(data = permuted_symbiont_quant, aes(y = factor(name, levels = trt_order), x = mean), shape = 18)+
  geom_point(data = real_symbiont_counts, aes(y = factor(name, levels = trt_order), x = value, color= real_or, fill = real_or, alpha = real_or, shape = real_or, size = real_or))+
  scale_color_manual(values = c("red", "black"))+
  scale_fill_manual(values = c("red", NA))+
  scale_alpha_manual(values = c(1,.1))+
  scale_shape_manual(values = c(16,21))+
  scale_size_manual(values = c(2,1))+
  facet_wrap(~factor(symbiont_type, levels = symbiont_rank) + factor(metric_category, levels = metric_rank), scales = "free")+
  scale_x_continuous(breaks = scales::breaks_pretty())+
  labs(x = "# of symbiota", y = "", color = "", fill = "", alpha = "", shape = "", size = "")+
  theme_bw()
permutation_metricXsymbiont_plot

ggsave(permutation_metricXsymbiont_plot, filename = "Plots/permutation_metricXsymbiont_plot.png", height = 5, width = 5)




######## repeating this within life stage, to see if this is driven by opposing effects across stages rather than within #####
symbiont_permute_df <-  permute_df %>% 
  mutate(symbiota = paste(host_species_clean, symbiont_species, sep = "_")) %>% 
  group_by(symbiont_type) %>% 
  filter(length(unique(symbiota))>1) %>% 
  rename(real=positive) %>% 
  select(host_family, host_species_clean, symbiont_species, symbiont_type, lifestage_general, metricXstage, symbiota, real) # re-ordering to put symbiota in front of "real"


symbiont_rank <- c("Bacteria", "Epichloe", "Other Fungi")
lifestage_rank <- unique(poa_df$lifestage_general)
symbiont_iter <- expand.grid(symbiont_rank = symbiont_rank,lifestage_rank = lifestage_rank)
permuted_symbiont <- list()
for(i in 1:nrow(symbiont_iter)){
  filtered_symbiont <- symbiont_permute_df[symbiont_permute_df$symbiont_type == symbiont_iter[i,"symbiont_rank"] & symbiont_permute_df$lifestage_general == symbiont_iter[i,"lifestage_rank"],]
  permuted_symbiont[[i]] <- as_tibble(cbind(filtered_symbiont, replicate(500, sample(filtered_symbiont$real))))  %>% 
    pivot_longer(cols = real:`...508`, names_to = "permutation", values_to = "effect")  %>% 
    group_by(host_family, host_species_clean, symbiont_species, lifestage_general, symbiont_type, permutation) %>% 
    summarize(n_metricsXstage = length(unique(metricXstage)),
              n_rows = n(),
              n_pos = sum(effect == "positive"),
              n_negative = sum(effect == "negative"),
              prop_pos = n_pos/n_rows,
              opposing = case_when(prop_pos == 1 ~ "all positive",
                                   prop_pos == 0 ~ "all negative",
                                   prop_pos >0 & prop_pos <1 ~ "opposing")) %>%  
    filter(n_metricsXstage >1) %>% 
    ungroup() %>% group_by( symbiont_type, lifestage_general, permutation) %>% 
    summarize(`All Positive` =  sum(opposing == "all positive"),
              `All Negative` = sum(opposing == "all negative"),
              `Opposing` = sum(opposing == "opposing"))
  
}

permuted_symbionts <- bind_rows(permuted_symbiont)



permuted_symbiont_long <- permuted_symbionts %>% 
  # filter(permutation != "real") %>% 
  pivot_longer(cols = c(`All Positive`, `All Negative`, `Opposing`)) %>% 
  mutate(real_or = case_when(permutation == "real" ~ "Obs. data",
                             permutation != "real" ~ "Permuted data"))

permuted_symbiont_counts <- permuted_symbiont_long %>% filter(permutation != "real")
real_symbiont_counts <- permuted_symbiont_long %>% filter(permutation == "real")



trt_order <- c("All Negative", "Opposing", "All Positive")


# calculating quantiles from permuted data

permuted_symbiont_quant <- permuted_symbiont_long %>% 
  filter(real_or == "Permuted data") %>% 
  group_by(symbiont_type, lifestage_general, name, real_or) %>% 
  summarize(mean = mean(value),
            Q02.5 = quantile(value, .025),
            Q97.5 = quantile(value, .975),
            n = n())

permutation_lifestageXsymbiont_plot <- ggplot(permuted_symbiont_counts)+
  geom_jitter(aes(y = factor(name, levels = trt_order), x = value, color= real_or, fill = real_or, alpha = real_or, shape = real_or, size = real_or), height = .3, width = 0.05)+
  geom_linerange(data = permuted_symbiont_quant, aes( y = factor(name, levels = trt_order), xmin = Q02.5, xmax = Q97.5), lwd = .3)+
  geom_point(data = permuted_symbiont_quant, aes(y = factor(name, levels = trt_order), x = mean), shape = 18)+
  geom_point(data = real_symbiont_counts, aes(y = factor(name, levels = trt_order), x = value, color= real_or, fill = real_or, alpha = real_or, shape = real_or, size = real_or))+
  scale_color_manual(values = c("red", "black"))+
  scale_fill_manual(values = c("red", NA))+
  scale_alpha_manual(values = c(1,.1))+
  scale_shape_manual(values = c(16,21))+
  scale_size_manual(values = c(2,1))+
  facet_wrap(~factor(symbiont_type, levels = symbiont_rank) + factor(lifestage_general, levels = lifestage_rank), scales = "free")+
  scale_x_continuous(breaks = scales::breaks_pretty())+
  labs(x = "# of symbiota", y = "", color = "", fill = "", alpha = "", shape = "", size = "")+
  theme_bw()
permutation_lifestageXsymbiont_plot

ggsave(permutation_lifestageXsymbiont_plot, filename = "Plots/permutation_lifestageXsymbiont_plot.png", height = 5, width = 5)

