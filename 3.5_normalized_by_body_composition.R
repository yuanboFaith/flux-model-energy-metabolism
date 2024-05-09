rm(list=ls())

library(ComplexHeatmap)
library(scales) # for color scale
library(mgcv)
library(ggpmisc)
library(ggbeeswarm)
library(readxl)
library(RColorBrewer)
library(viridis)
library(ggsci)
library(rebus)
library(matrixStats)
library(limSolve)
library(xlsx)
library(cowplot)
library(rstatix)
library(ggpubr)
library(tidyverse)

# load the theme
load("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/1_core_13CO2_import_data.RData")

# dataset of body composition by MRI
d.bodyComposition <- read_excel(
  "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/body composition.xlsx") %>% 
  mutate(phenotype = factor(phenotype, levels = ordered.phenotype))
  

# mean and SD of mass of fat and lean tissues
d.bodyComposition.tidy <- d.bodyComposition %>%  
  group_by(phenotype) %>% 
  summarise(across(c(fat, lean),
                   .fn = list(mean = ~mean(.x), sd = ~sd(.x)))
  ) %>% 
  pivot_longer(-phenotype, names_sep = "_",
               names_to = c("part", ".value")) %>% 
  # error bar position
  group_by(phenotype) %>% 
  arrange(desc(part)) %>% 
  mutate(y.error = cumsum(mean))

# plot
plt.fat.lean <- d.bodyComposition.tidy %>% 
  ggplot(aes(x = phenotype, y = mean, fill = part)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = y.error - sd, ymax = y.error),
                width = .2) +
  theme.myClassic +
  scale_y_continuous(expand = expansion(mult = 0),
                     breaks = seq(0, 70, 10)) +
  scale_fill_manual(values = c("fat" = "snow2", "lean" = "snow4")) +
  labs(y = "mass (g)")

plt.fat.lean

# calculate fraction of fat and lean - ------------------------
d.bodyComposition.summary <- d.bodyComposition %>% 
  mutate(sum = fat + lean,
         fat.frac = fat/sum,
         lean.frac = lean / sum) %>% 
  group_by(phenotype) %>% 
  summarise(across(contains("frac"), 
                   .fn = list(
                     mean = ~ mean(.x, na.rm = T), # calcualte the mean
                     sd = ~ sd(.x, na.rm = T))))  # calculate standard deviation
# tidy up
d.bodyComposition.summary.tidy <- d.bodyComposition.summary %>% 
  pivot_longer(-phenotype, names_to = c("part", ".value"), names_sep = "_") %>% 
  group_by(phenotype) %>% 
  arrange(desc(part)) %>% 
  mutate(y.error = cumsum(mean)) %>% 
  mutate(part = str_remove(part, ".frac"))

# plot of distribution fraction
plt.fat.lean.frac <- d.bodyComposition.summary.tidy %>% 
  ggplot(aes(x = phenotype, y = mean, fill = part)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = y.error - sd, ymax = y.error),
                width = .2) +
  theme.myClassic +
  scale_y_continuous(expand = expansion(mult = 0),
                     breaks = seq(0, 1, .2)) +
  scale_fill_manual(values = c("fat" = "snow2", "lean" = "snow4")) +
  labs( y = "fraction")

plt.fat.lean.frac

plot_grid(plt.fat.lean, plt.fat.lean.frac)



# save project
save(d.bodyComposition.summary.tidy, 
     d.bodyComposition.tidy, 
     file = "3.5_core_13CO2_infusion_physiology_basics.RData")



