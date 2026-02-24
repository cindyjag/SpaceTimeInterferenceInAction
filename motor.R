################################################################################
# Title: Motor Task Analysis
# Description: Data processing and analysis pipeline for the Motor Task
# Author: Cindy Jagorska
# Date: 08.10.2025
# Paper: Space-Time Interference in Action
################################################################################

# ==============================#
#         Load Libraries        #
# ==============================#
library(tidyverse)
library(lme4)
library(lmerTest)
library(broom)
library(sjPlot)
library(ggplot2)
library(ggeffects)
library(emmeans)
library(performance)


# ==============================#
#          Data Import          #
# ==============================#

all_data <- read.table("D:/motor_data.txt", header = TRUE, sep = "\t")


# ==============================#
#       Data Preprocessing      #
# ==============================#

all_data <- all_data %>%
  mutate(
    rTime = reaTime - timelastBall,
    timeCon = case_when(
      time1 == 0.3 ~ "short",
      time1 == 0.9 ~ "middle",
      TRUE ~ "long"
    ),
    dist = round(pos3 - pos2, 2),
    distCon = case_when(
      dist == 0.05 ~ "short",
      dist == 0.10 ~ "middle",
      TRUE ~ "long"
    ),
    tcondition = recode(timeCon, short = 0, middle = 1, long = 2),
    scondition = recode(distCon, short = 0, middle = 1, long = 2)
  )

# ==============================#
#           Ratios              #
# ==============================#

all_data <- all_data %>%
  mutate(
    expectedX = pos3 + (pos3 - pos2),
    timeRatio = rTime / time1,
    xRatio = ansx / expectedX
  )

# ==============================#
#         Outlier Removal       #
# ==============================#

all_data_cleaned <- all_data %>%
  group_by(timeCon) %>%
  filter(
    abs(timeRatio - median(timeRatio, na.rm = TRUE)) <=
      5 * mad(timeRatio, constant = 1.4826, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  group_by(distCon) %>%
  filter(
    abs(xRatio - median(xRatio, na.rm = TRUE)) <=
      5 * mad(xRatio, constant = 1.4826, na.rm = TRUE)
  ) %>%
  ungroup()

# ==============================#
#      Linear Mixed Models      #
# ==============================#

# --- Model for timeRatio ---
model_timeRatio <- lmer(timeRatio ~ scondition + tcondition + (tcondition | subj) ,
                        data = all_data_cleaned)

model_timeRatio_qu <- lmer(timeRatio ~ scondition  + I(scondition^2)+ tcondition + (tcondition | subj) ,
                        data = all_data_cleaned)

anova(model_timeRatio, model_timeRatio_qu)


summary(model_timeRatio)
tab_model(model_timeRatio)
res_timeRatio <- simulateResiduals(
  fittedModel = model_timeRatio,
  n = 1000
)
plot(res_timeRatio)

# --- Subject-level betas for timeRatio ---
beta_int <- all_data_cleaned %>%
  group_by(subj) %>%
  group_modify(~ tidy(lm(timeRatio ~ scondition, data = .x))) %>%
  dplyr::select(subj, term, estimate)

# ==============================#
#       Group-level Means       #
# ==============================#

summary_data <- all_data_cleaned %>%
  group_by(scondition) %>%
  summarise(
    mean_timeRatio = mean(timeRatio, na.rm = TRUE),
    se_timeRatio = sd(timeRatio, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# ==============================#
#            Plot 1             #
# ==============================#

large_text_theme <- theme(
  text = element_text(size = 14),
  plot.title = element_text(size = 14),
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12),
  legend.position = "none"
)

ggplot(summary_data, aes(x = scondition, y = mean_timeRatio, color = factor(scondition))) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_timeRatio - se_timeRatio,
                    ymax = mean_timeRatio + se_timeRatio),
                width = 0.2, size = 1) +
  #geom_abline(intercept = 1.153, slope = 0.0143, color = "#00A08A", size = 1) +
  labs(x = "Spatial Interstimulus Interval", y = "Ratio of Temporal Estimate") +
  coord_cartesian(ylim = c(1.125, 1.23)) +
  theme_minimal() +
  large_text_theme +
  scale_color_manual(values = c("#E2D200", "#F98400", "#FF0000")) +
  scale_x_continuous(breaks = 0:2, labels = c("short", "medium", "long"))

# ==============================#
#        Model for xRatio       #
# ==============================#

model_xRatio <- lmer(xRatio ~ tcondition + (1 | subj) + (1 | scondition),
                     data = all_data_cleaned)

model_xRatio_qu <- lmer(xRatio ~ tcondition + + I(tcondition^2) + (1 | subj) + (1 | scondition),
                     data = all_data_cleaned)

anova(model_xRatio, model_xRatio_qu)

r2(model_xRatio)
summary(model_xRatio)
tab_model(model_xRatio)

res_xRatio <- simulateResiduals(
  fittedModel = model_xRatio,
  n = 1000
)
plot(res_xRatio)

# --- Subject-level betas for xRatio ---
betas_int <- all_data_cleaned %>%
  group_by(subj) %>%
  group_modify(~ tidy(lm(xRatio ~ tcondition, data = .x))) %>%
  dplyr::select(subj, term, estimate)

# ==============================#
#        Correlation Tests      #
# ==============================#

# Merge perception and interception betas if available
# (Assumes `beta_perc` and `betas_perc` are loaded from perception analysis)
 all_beta_space <- merge(beta_space, betas_int, by = c("subj", "term"))
 all_beta_space<-all_beta_space %>%
   filter(subj != 25)
 cor.test(all_beta_space[all_beta_space$term == "tcondition",]$estimate.x, all_beta_space[all_beta_space$term == "tcondition",]$estimate.y)
 
 all_beta_time <- merge(beta_time, beta_int, by = c("subj", "term"))
 cor.test(all_beta_time[all_beta_time$term == "scondition",]$estimate.x, all_beta_time[all_beta_time$term == "scondition",]$estimate.y)




 ggplot(all_beta_time[all_beta_time$term == "scondition",], aes(x = estimate.x, y = estimate.y)) +
  geom_point(size = 3, alpha = 0.85, color = "darkgrey") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(title = "space-on-time effect",
    x = "individual betas perceptual task",
    y = "individual betas motor task"
  ) +
  theme_minimal(base_size = 14)+
  large_text_theme
 
 
 ggplot(all_beta_space[all_beta_space$term == "tcondition",], aes(x = estimate.x, y = estimate.y)) +
   geom_point(size = 3, alpha = 0.85, color = "darkgrey") +
   geom_smooth(method = "lm", se = FALSE, color = "black") +
   labs(title = "time-on-space effect",
     x = "individual betas perceptual task",
     y = "individual betas motor task"
   ) +
   theme_minimal(base_size = 14)+
   large_text_theme
 

 
 
 plot_data <- bind_rows(
   all_beta_time %>%
     filter(term == "scondition") %>%
     mutate(effect = "space-on-time effect"),
   all_beta_space %>%
     filter(term == "tcondition") %>%
     mutate(effect = "time-on-space effect")
 )
 
 ggplot(plot_data, aes(x = estimate.x, y = estimate.y)) +
   geom_point(
     size = 3,
     alpha = 0.8,
     color = "grey40"
   ) +
   geom_smooth(
     method = "lm",
     se = FALSE,
     color = "black",
     linewidth = 0.8
   ) +
   facet_wrap(~ effect) +
   labs(
     x = "linear beta estimates (perceptual task)",
     y = "linear beta estimates (motor task)"
   ) +
   theme_minimal(base_size = 14) +
   large_text_theme +
   theme(
     strip.text = element_text(face = "bold"),
     panel.grid.minor = element_blank()
   )
 
 
 
 
# ==============================#
#       Group Means for xRatio  #
# ==============================#

summary_data2 <- all_data_cleaned %>%
  group_by(tcondition) %>%
  summarise(
    mean_spaceRatio = mean(xRatio, na.rm = TRUE),
    se_spaceRatio = sd(xRatio, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# ==============================#
#            Plot 2             #
# ==============================#

ggplot(summary_data2, aes(x = tcondition, y = mean_spaceRatio, color = factor(tcondition))) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_spaceRatio - se_spaceRatio,
                    ymax = mean_spaceRatio + se_spaceRatio),
                width = 0.2, size = 1) +
  #geom_abline(intercept = 1.208, slope = -0.0174, color = "#00A08A", size = 1) +
  labs(x = "Temporal Interstimulus Interval", y = "Ratio of Spatial Estimate") +
  coord_cartesian(ylim = c(1.125, 1.23)) +
  theme_minimal() +
  large_text_theme +
  scale_color_manual(values = c("#E2D200", "#F98400", "#FF0000")) +
  scale_x_continuous(breaks = 0:2, labels = c("short", "medium", "long"))



################################################################################
# End of Script
################################################################################
