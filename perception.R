################################################################################
# Title: Perception Task Analysis
# Description: Data processing and analysis pipeline for the Perception Task
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
library(DHARMa)


# ==============================#
#        Data Import            #
# ==============================#





all_data <- read.table("D:/perceptual_data.txt", header = TRUE, sep = "\t")

# ==============================#
#      Data Preprocessing       #
# ==============================#

# Response time
all_data <- all_data %>%
  mutate(
    rTime = reaTime - timeText,
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
    ans = if_else(rX < 0, 0, 1),
    tcondition = recode(timeCon, short = 0, middle = 1, long = 2),
    scondition = recode(distCon, short = 0, middle = 1, long = 2),
    task = case_when(
      task == "viz.VizText(49)" ~ "time",
      task == "viz.VizText(51)" ~ "time",
      TRUE ~ "space"
    )
  )

# Remove outliers (rTime > 5)
all_data <- all_data %>% filter(rTime <= 5)


# ==============================#
#    Descriptive Statistics     #
# ==============================#

## Mean per participant and condition
agg_data <- all_data %>%
  group_by(subj, timeCon, distCon, task) %>%
  summarise(
    mean_ans = mean(ans),
    mean_rTime = mean(rTime),
    sd_ans = sd(ans),
    sd_rTime = sd(rTime),
    .groups = "drop"
  )

## Mean across participants
desc_summary <- all_data %>%
  group_by(task, timeCon, distCon) %>%
  summarise(
    mean_ans = mean(ans),
    mean_rTime = mean(rTime),
    se_ans = sd(ans) / sqrt(n()),
    se_rTime = sd(rTime) / sqrt(n()),
    .groups = "drop"
  )

# ==============================#
#     Linear Mixed Models       #
# ==============================#


# ---- TIME task ----
##standardize for quadratic function
all_data_t <- all_data %>%
  mutate(
    scondition_c = scondition - mean(scondition, na.rm = TRUE),
    scondition_c2 = scondition_c^2
  )


model_time_qu <- glmer(
  ans ~ scondition_c + scondition_c2 +
    (1 | subj)+ (1 | tcondition),
  data = filter(all_data_t, task == "time"),
  family = binomial
)

model_time_li <- glmer(
  ans ~ scondition +
    (1 | subj)+ (1 | tcondition),
  data = filter(all_data_t, task == "time"),
  family = binomial
)
anova(model_time_li, model_time_qu)

summary(model_time_qu)
tab_model(model_time_qu)




# ---- SPACE task ----
all_data_s <- all_data %>%
  mutate(
    tcondition_c = tcondition - mean(tcondition, na.rm = TRUE),
    tcondition_c2 = tcondition_c^2
  )


model_space_qu <- glmer(
  ans ~ tcondition_c + tcondition_c2 +
    (1 | subj)+ (1 | scondition),
  data = filter(all_data_s, task == "space"),
  family = binomial
)

model_space_li <- glmer(
  ans ~ tcondition +
    (1 | subj)+ (1 | scondition),
  data = filter(all_data_s, task == "space"),
  family = binomial
)


anova(model_space_li, model_space_qu)

summary(model_space_qu)
tab_model(model_space_qu)



res_space <- simulateResiduals(
  fittedModel = model_space_qu,
  n = 1000
)
plot(res_space)

testDispersion(res_space)

# ==============================#
#     Extract Beta Estimates    #
# ==============================#




# --- TIME betas ---

beta_time <- all_data %>%
  filter(task == "time") %>%
  group_by(subj) %>%
  group_modify(~ 
                 tidy(glm(ans ~ scondition,
                          data = .x,
                          family = binomial))
  ) %>%
  select(subj, term, estimate)



# --- SPACE betas ---

beta_space <- all_data %>%
  filter(task == "space") %>%
  group_by(subj) %>%
  group_modify(~ 
                 tidy(glm(ans ~ tcondition,
                          data = .x,
                          family = binomial))
  ) %>%
  select(subj, term, estimate)


# ==============================#
#            Plots              #
# ==============================#

large_text_theme <- theme(
  text = element_text(size = 12),
  plot.title = element_text(size = 12),
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 10),
  axis.title.x = element_text(margin = margin(t = 20)),
  legend.position = "none"
)

# ---- TIME plot ----
plot_time <- all_data %>%
  filter(task == "time") %>%
  group_by(scondition) %>%
  summarise(
    mean_ans = mean(ans),
    se_ans = sd(ans) / sqrt(n()),
    .groups = "drop"
  )

ggplot(plot_time, aes(x = scondition, y = mean_ans, color = factor(scondition))) +
  geom_point(size = 3) +
  #geom_abline(intercept = 0.5639, slope = -0.0159, color = "#00A08A", size = 1) +
  geom_errorbar(aes(ymin = mean_ans - se_ans, ymax = mean_ans + se_ans), width = 0.2, size = 1) +
  labs(x = "Spatial Interstimulus Interval", y = "% too late") +
#  geom_smooth( aes(x = scondition, y = ans), color = "#00A08A", method="lm", formula = y ~ x + I(x^2),se = FALSE, data = filter(all_data, task == "time"))+
  theme_minimal() +
  large_text_theme +
  coord_cartesian(ylim = c(0.45, 0.65)) +
  scale_color_manual(values = c("#E2D200", "#F98400", "#FF0000")) +
  scale_x_continuous(breaks = 0:2, labels = c("short", "medium", "long"))

# ---- SPACE plot ----
plot_space <- all_data %>%
  filter(task == "space") %>%
  group_by(tcondition) %>%
  summarise(
    mean_ans = mean(ans),
    se_ans = sd(ans) / sqrt(n()),
    .groups = "drop"
  )

library(ggeffects)

pred_space <- ggpredict(
  model_space,
  terms = "tcondition"
)

ggplot(plot_space, aes(x = tcondition, y = mean_ans, color = factor(tcondition))) +
  geom_point(size = 3) +
 # geom_line(
  #  data = pred_space,
   # aes(x = x, y = predicted),
    #color = "#00A08A",
#    size = 1.2
 # ) +
  geom_errorbar(aes(ymin = mean_ans - se_ans, ymax = mean_ans + se_ans), width = 0.2, size = 1, data = plot_space) +
  labs(x = "Temporal Interstimulus Interval", y = "% too far right") +
  theme_minimal() +
  large_text_theme +
  coord_cartesian(ylim = c(0.45, 0.65)) +
  scale_color_manual(values = c("#E2D200", "#F98400", "#FF0000")) +
  scale_x_continuous(breaks = 0:2, labels = c("short", "medium", "long"))
# ==============================#
#             Tests             #
# ==============================#

## One-sample t-tests against 0.5
t.test(all_data %>% filter(task == "space") %>% group_by(subj) %>%
         summarise(mean_ans = mean(ans)) %>% pull(mean_ans),
       mu = 0.5)

t.test(all_data %>% filter(task == "time") %>% group_by(subj) %>%
         summarise(mean_ans = mean(ans)) %>% pull(mean_ans),
       mu = 0.5)






################################################################################
# End of Script
################################################################################
