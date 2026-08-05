#compiling information for Table S2, presenting the list of citations in the meta analysis and information about vital rates
# Author Joshua Fowler
# Date: Aug 3, 2026

library(tidyverse)
library(patchwork)
# setting up color schemes
metric_colors <- c("#77AADD", "#EE8866", "#EEDD88", "#FFAABB")
stage_colors <- c("#44BB99", "#BBCC33", "#AAAA00", "#99DDFF")



effects_df <- read_csv("effects_df.csv") %>% 
  filter(experiment_id != "112-1") %>% 
  filter(!is.na(sd_RII)) %>% 
  mutate(RII = case_when(RII >= (1.0 - 1e-7) ~ 1.0,
                         RII <= (-1.0 + 1e-7) ~ -1.0, TRUE ~ RII),
         RII_01 = (RII+1)/2)


effects_summary <- effects_df %>% 
  group_by(study_number) %>% 
  summarize(`Vital Rates` = paste(unique(metric_category), collapse = ", "),
            `Life Stages` = paste(unique(lifestage_general), collapse = ", "), 
            `Host Family` = paste(unique(host_family), collapse = ", "),
            `Symbiont Family` =  paste(unique(symbiont_family), collapse = ", "),)

web_of_science <- read_csv("WOS_all_records.csv") %>% 
  filter(study_number %in% unique(effects_df$study_number)) %>% 
  select(study_number, Authors, `Publication Year`, `Source Title`, DOI)

table_S2 <- web_of_science %>% left_join(effects_summary) %>% select(-study_number)
write_csv(table_S2, "table_S2.csv")



#### evaluating sample size by calculating saturation of meta-analytic mean
effects_df <- read_csv("effects_df.csv") %>% 
  filter(experiment_id != "112-1") %>% 
  filter(!is.na(sd_RII)) %>% 
  mutate(RII = case_when(RII >= (1.0 - 1e-7) ~ 1.0,
                         RII <= (-1.0 + 1e-7) ~ -1.0, TRUE ~ RII),
         RII_01 = (RII+1)/2)


# full_means <- list()
VR_means <- LS_means <- list()

ids <- unique(effects_df$study_number)

cumulative_ids <- lapply(seq_along(ids), function(i) ids[1:i])

for(i in 1:length(cumulative_ids)){
  df <- effects_df %>% filter(study_number %in% cumulative_ids[[i]])
  # full_means[[i]] <- df %>% 
  #   # group_by(metric_category) %>% 
  #   summarize(mean = mean(RII),
  #             sd = sd(RII),
  #             var = sd(RII)^2) %>% 
  #   mutate(study_count = i)
  
  VR_means[[i]] <- df %>% 
    group_by(metric_category) %>% 
    summarize(mean = mean(RII),
              sd = sd(RII),
              var = sd(RII)^2) %>% 
    mutate(study_count = i)
  
  LS_means[[i]] <- df %>% 
    group_by(lifestage_general) %>% 
    summarize(mean = mean(RII),
              sd = sd(RII),
              var = sd(RII)^2) %>% 
    mutate(study_count = i)
} 

# full_means_df <- bind_rows(full_means)
VR_means_df <- bind_rows(VR_means)
LS_means_df <- bind_rows(LS_means)

# full_means <- ggplot(full_means_df) +
#   geom_line(aes(x = study_count, y = mean)) +
#   geom_point(aes(x = study_count, y = mean), shape = 21, fill = "white") +
#   # facet_wrap(~metric_category, nrow = 1) +
#   # scale_color_manual(values = metric_colors)+
#   labs(y = "Group Means", x = "Cumulative Study Count",  color = "Vital Rate", title = "Vital Rate - Mean in RII")+
#   theme_minimal()


VR_means <- ggplot(VR_means_df) +
  geom_line(aes(x = study_count, y = mean, color = metric_category)) +
  geom_point(aes(x = study_count, y = mean, color = metric_category), shape = 21, fill = "white", size = .5) +
  facet_wrap(~metric_category, nrow = 1) +
  scale_color_manual(values = metric_colors)+
  labs(y = "Group Means", x = "Cumulative Study Count",  color = "Vital Rate", title = "Vital Rate - Mean in RII")+
  theme_minimal()
  
VR_sd <- ggplot(VR_means_df) +
  geom_line(aes(x = study_count, y = var, color = metric_category)) +
  geom_point(aes(x = study_count, y = var, color = metric_category), shape = 21, fill = "white", size = .5) +
  facet_wrap(~metric_category, nrow = 1) +
  scale_color_manual(values = metric_colors)+
  labs(y = "Group Var.", x = "Cumulative Study Count", color = "Vital Rate", title = "Vital Rate - Variance in RII")+
  theme_minimal()

LS_means <- ggplot(LS_means_df) +
  geom_line(aes(x = study_count, y = mean, color = lifestage_general)) +
  geom_point(aes(x = study_count, y = mean, color = lifestage_general), shape = 21, fill = "white", size = .5) +
  facet_wrap(~lifestage_general, nrow = 1) +
  scale_color_manual(values = stage_colors)+
  labs(y = "Group Means", x = "Cumulative Study Count", color = "Life Stage", title = "Life Stage - Mean in RII")+
  theme_minimal()

LS_sd <- ggplot(LS_means_df) +
  geom_line(aes(x = study_count, y = var, color = lifestage_general)) +
  geom_point(aes(x = study_count, y = var, color = lifestage_general), shape = 21, fill = "white", size = .5) +
  facet_wrap(~lifestage_general, nrow = 1) +
  scale_color_manual(values = stage_colors)+
  labs(y = "Group Var.", x = "Cumulative Study Count", color = "Life Stage", title = "Life Stage - Variance in RII")+
  theme_minimal()


cum_means_plot <- VR_means / VR_sd / LS_means / LS_sd + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
ggsave(cum_means_plot, filename = "Plots/cum_means_plot.png", width = 7, height = 10)



# version with randomly selecting individual studies
num_samples <- 100

lengths <- 1:length(ids)

sampled_list <- lapply(lengths, function(len) {
  replicate(num_samples, sample(ids, size = len, replace = FALSE), simplify = FALSE)
})

flat_id_list <- unlist(sampled_list, recursive = FALSE)

sampled_df <- lapply(flat_id_list, function(ids){
  effects_df[effects_df$study_number%in% ids,]
})


VR_sampled_means <- lapply(sampled_df, function(df) {
                        df %>%     
                             group_by(metric_category) %>% 
                             summarize(mean = mean(RII), sd = sd(RII),var = sd(RII)^2) %>% 
                             mutate(study_count = length(unique(df$study_number)))
  })


VR_sampled_means_df <- bind_rows(VR_sampled_means, .id = "id")

VR_sampled_means_summary <- VR_sampled_means_df %>% 
  group_by(metric_category, study_count) %>% 
  summarize(mean_summary = mean(mean),
            sd_summary = mean(sd),
            q2.5 = quantile(mean, probs = c(0.025)),
            q97.5 = quantile(mean, probs = c(0.975)))

VR_means <- ggplot(VR_sampled_means_df) +
  # geom_point(aes(x = study_count, y = mean), color = "black", shape = 21, fill = "white", size = .5) +
  geom_ribbon(data = VR_sampled_means_summary, aes(x = study_count, ymin = q2.5, ymax = q97.5, fill = metric_category), alpha = .5) +
  geom_line(data = VR_sampled_means_summary, aes(x = study_count, y = mean_summary,  color = metric_category)) +
  facet_wrap(~metric_category, nrow = 1) +
  scale_fill_manual(values = metric_colors)+
  scale_color_manual(values = metric_colors)+
  labs(y = "Group Means", x = "Cumulative Study Count",  color = "Vital Rate", fill = "Vital Rate", title = "Vital Rate - Mean RII")+
  theme_minimal()

VR_means






LS_sampled_means <- lapply(sampled_df, function(df) {
  df %>%     
    group_by(lifestage_general) %>% 
    summarize(mean = mean(RII), sd = sd(RII),var = sd(RII)^2) %>% 
    mutate(study_count = length(unique(df$study_number)))
})


LS_sampled_means_df <- bind_rows(LS_sampled_means, .id = "id")

LS_sampled_means_summary <- LS_sampled_means_df %>% 
  group_by(lifestage_general, study_count) %>% 
  summarize(mean_summary = mean(mean),
            sd_summary = mean(sd),
            q2.5 = quantile(mean, probs = c(0.025)),
            q97.5 = quantile(mean, probs = c(0.975)))

LS_means <- ggplot(LS_sampled_means_df) +
  # geom_point(aes(x = study_count, y = mean), color = "black", shape = 21, fill = "white", size = .5) +
  geom_ribbon(data = LS_sampled_means_summary, aes(x = study_count, ymin = q2.5, ymax = q97.5, fill = lifestage_general), alpha = .5) +
  geom_line(data = LS_sampled_means_summary, aes(x = study_count, y = mean_summary,  color = lifestage_general)) +
  facet_wrap(~lifestage_general, nrow = 1) +
  scale_fill_manual(values = stage_colors)+
  scale_color_manual(values = stage_colors)+
  labs(y = "Group Means", x = "Cumulative Study Count",  color = "Life Stage", fill = "Life Stage", title = "Life Stage - Mean RII")+
  theme_minimal()

LS_means



rand_cum_means_plot <- VR_means / LS_means + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
ggsave(rand_cum_means_plot, filename = "Plots/rand_cum_means_plot.png", width = 7, height = 8)

