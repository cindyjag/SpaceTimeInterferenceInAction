# SpaceTimeInterferenceInAction

# Temporal and Spatial Estimation Tasks – Data and Analysis

This repository contains the data and R analysis scripts for the **Perception** and **motor** tasks described in our study on the relationship between temporal and spatial estimation in action.

├── data/
│ ├── perceptual_data.txt # Combined data file for Perception Task
│ ├── motor_data.txt # Combined data file for Interception Task
│
├── scripts/
│ ├── perception.R # Full analysis pipeline for Perception Task
│ ├── motor.R # Full analysis pipeline for Motor Task
│
├── README.md # This file

## ⚙️ Requirements

These analyses were run in **R (≥ 4.3)** using the following packages:

```r
tidyverse
lme4
lmerTest
broom
sjPlot
ggeffects
emmeans

You can install them all at once with:

install.packages(c(
  "tidyverse", "lme4", "lmerTest", "broom",
  "sjPlot", "ggeffects", "emmeans"
))
