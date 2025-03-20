# Title: How commonly do microbial effects vary across host life stages? #####
# Purpose: Reads in data from literature search and tests how commonly vital rate effects are of opposing signs
# Authors: Josh Fowler and Gwen Pohlmann #####
# Date: Mar 18, 2025#####

# library(renv)
# renv::init()
# renv::snapshot()
# renv::deactivate()


library(tidyverse)


####### Reading in the data   #######


effects_df <- read_csv("effects_df.csv")



####### calculating the proportion of opposing symbiotic effects within a particular study  #######
# first we need to winnow down the data set to only those that have measurements of multiple metrics
# also note that I need to investigate the cases with RII exactly = 0

multiple_metrics_df <- effects_df %>%
  filter(RII != 0) %>% 
  mutate(positive = case_when(RII>0 ~ "positive", 
                              RII <0 ~ "negative")) %>% 
  mutate(metricXstage = paste0(metric_description, lifestage_description, sep = "_")) %>% 
  group_by(host_species_clean, symbiont_species) %>% 
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
ggplot(multiple_metrics_df)+
  geom_histogram(aes(x = factor(opposing, level = trt_order), fill = opposing), stat = "count")+
  theme_classic()



####### simulating a null distribution by permuting assigned microbial effects #######
permute_df <- effects_df %>% 
  filter(RII != 0) %>% 
  mutate(metricXstage = paste0(metric_description, lifestage_description, sep = "_")) %>% 
  mutate(positive = case_when(RII>0 ~ "positive", 
                              RII <0 ~ "negative")) %>% 
  mutate(symbiota_id = paste0(host_species_clean, symbiont_species, sep = "_")) %>% 
  filter(symbiota_id %in% unique(multiple_metrics_df$symbiota_id)) %>%  
  dplyr::select(host_species_clean, symbiont_species,metricXstage, positive)
  
permuted <- as_tibble(cbind(permute_df, replicate(500, sample(permute_df$positive)))) %>% 
  rename(real=positive) %>% 
  pivot_longer(cols = real:`500`, names_to = "permutation", values_to = "effect")  %>%  
  group_by(host_species_clean, symbiont_species, permutation) %>% 
  summarize(n_metricsXstage = length(unique(metricXstage)),
            n_rows = n(),
            n_pos = sum(effect == "positive"),
            n_negative = sum(effect == "negative"),
            prop_pos = n_pos/n_rows,
            opposing = case_when(prop_pos == 1 ~ "all positive",
                                 prop_pos == 0 ~ "all negative",
                                 prop_pos >0 & prop_pos <1 ~ "opposing")) %>%  
  filter(n_metricsXstage >1) %>% 
  ungroup() %>% group_by(permutation) %>% 
  summarize(`All Positive` =  sum(opposing == "all positive"),
            `All Negative` = sum(opposing == "all negative"),
            `Opposing` = sum(opposing == "opposing"))
  

permuted_long <- permuted %>% 
  # filter(permutation != "real") %>% 
  pivot_longer(cols = c(`All Positive`, `All Negative`, `Opposing`)) %>% 
  mutate(real_or = case_when(permutation == "real" ~ "Obs. data",
                             permutation != "real" ~ "Permuted data"))

permuted_counts <- permuted_long %>% filter(permutation != "real")
real_counts <- permuted_long %>% filter(permutation == "real")



trt_order <- c("All Negative", "Opposing", "All Positive")

# ggplot(permuted_long)+
#   geom_histogram(aes(x = value, fill= name), alpha = .2, bins = 100)+
#   geom_vline(data = real_counts, aes(xintercept = value, color = name))+
#   facet_wrap(~factor(name, levels = trt_order), ncol = 1, strip.position = "top")+
#   labs(x = "# of symbiota")+
#   guides(fill = "none", color = "none")+
#   theme_minimal()+
#   theme(strip.background = element_blank())
  
  

permutation_plot <- ggplot(permuted_counts)+
  geom_jitter(aes(y = factor(name, levels = trt_order), x = value, color= real_or, fill = real_or, alpha = real_or, shape = real_or, size = real_or), height = .3)+
  geom_point(data = real_counts, aes(y = factor(name, levels = trt_order), x = value, color= real_or, fill = real_or, alpha = real_or, shape = real_or, size = real_or))+
  scale_color_manual(values = c("red", "black"))+
  scale_fill_manual(values = c("red", NA))+
  scale_alpha_manual(values = c(1,.2))+
  scale_shape_manual(values = c(16,21))+
  scale_size_manual(values = c(2,1))+
  labs(x = "# of symbiota", y = "", color = "", fill = "", alpha = "", shape = "", size = "")+
  theme_bw()
permutation_plot

ggsave(permutation_plot, filename = "Plots/permutation_plot.png", height = 2, width = 6)






# calculating probabilities


mean(permuted_counts[permuted_counts$name == "All Positive",]$value,)
quantile(permuted_counts[permuted_counts$name == "All Positive",]$value, .025)
quantile(permuted_counts[permuted_counts$name == "All Positive",]$value, .975)
real_counts[real_counts$name == "All Positive",]


mean(permuted_counts[permuted_counts$name == "Opposing",]$value,)
quantile(permuted_counts[permuted_counts$name == "Opposing",]$value, .025)
quantile(permuted_counts[permuted_counts$name == "Opposing",]$value, .975)
real_counts[real_counts$name == "Opposing",]


mean(permuted_counts[permuted_counts$name == "All Negative",]$value,)
quantile(permuted_counts[permuted_counts$name == "All Negative",]$value, .025)
quantile(permuted_counts[permuted_counts$name == "All Negative",]$value, .975)
real_counts[real_counts$name == "All Negative",]
