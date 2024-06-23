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

load("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/7_core_production_flux.RData")

# hyperinsulinemia clamp 
d.hyperInsulinemia <- d.normalized.tidy.23 %>% 
  mutate(mouse_ID = str_sub(sample, -4)) %>% 
  separate(mouse_ID, into = c("time", "mouse_ID")) %>% 
  group_by(Compound) %>% 
  mutate(enrich.weighted = C_Label / max(C_Label) * enrichment) %>% 
  group_by(mouse_ID, time, .add = T) %>% 
  summarise(enrich = sum(enrich.weighted))

d.hyperInsulinemia

# see additional processing in Excel 

d.clamp <- read_excel("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/hyperinsulinemia.xlsx",
                      range = "A15:F45")  
d.clamp  

# extract endogenous production of glucose and palmitate in WT
d.clamp.ALL <- d.Fcirc.standard.atom_art.or.tail %>% 
  filter(phenotype == "WT" & Compound %in% c("Glucose", "C16:0")) %>% 
  select(Compound, Fcirc_atom.g.BW, mouse_ID) %>% 
  mutate(treatment = "basal", timepoint = "T0", Index = "Endo_Ra") %>% 
  relocate(treatment, mouse_ID, timepoint, Fcirc_atom.g.BW, Compound, Index) %>% 
  # combine with clamp data
  bind_rows(d.clamp)


func.plt.clamp <- function(data, myTitle){
  data %>% 
    ggplot(aes(x = timepoint, y = Fcirc_atom.g.BW, fill = timepoint)) +
    stat_summary(geom = "bar", fun = mean, color = "black", alpha = .7) +
    stat_summary(geom = "errorbar", fun.data = mean_sdl, fun.args = list(mult = 1),
                 width = .3) +
    geom_quasirandom(width = .2, size = 2, color = "grey20") +
    facet_grid(.~Compound, scales = "free", space = "free") +
    theme.myClassic +
    theme(axis.title.y = element_text(margin = margin(r = 15)),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, .1)), name = myTitle,
                       # limits = c(0, 2100)
    ) +
    scale_x_discrete(labels = c("T0" = "basal", "T1" = "Clamp T1", "T2" = "Clamp T2"), name = NULL,
                     expand = expansion(add = .8)) +
    scale_fill_manual(values = c("steelblue2", "tomato", "tomato"))
}

# C16:0 endogenous production
p1 <- d.clamp.ALL %>% filter(Index == "Endo_Ra" & Compound == "C16:0") %>% 
  func.plt.clamp(myTitle = "Endogenous Ra = Rd (nmol C/min/g)")
p1

# Glucose endogenous production
p2 <- d.clamp.ALL %>% filter(Index == "Endo_Ra" & Compound == "Glucose") %>% 
  func.plt.clamp(myTitle = "Endogenous Ra (nmol C/min/g)")
p2

# Glucose total consumption
p3 <- d.clamp.ALL %>% 
  filter(Index == "Total Ra" & Compound == "Glucose") %>% 
  func.plt.clamp(myTitle = "Rd (nmol C/min/g)")
p3

plot_grid(p1, p2, p3, rel_widths = c(2, 2.8, 2), nrow = 1)

ggsave(filename = "Clamp.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 4, width = 8)

# t test of lipolysis suppression
d.clamp.ALL %>% filter(Index == "Endo_Ra" & Compound == "C16:0") %>% 
  t.test(Fcirc_atom.g.BW ~ treatment, data = . )

# t test of glucose endogenous production suppression
d.clamp.ALL %>% filter(Index == "Endo_Ra" & Compound == "Glucose") %>% 
  t.test(Fcirc_atom.g.BW ~ treatment, data = . )


# Plot glycemia
d.gly.ins <- read_excel(
  path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/hyperinsulinemia.xlsx",
  sheet = "insulin_glycemia")

d.gly.ins %>% 
  ggplot(aes(x = time, y = conc., fill = status)) +
  stat_summary(geom = "bar", fun = mean, color = "black", alpha = .7) +
  stat_summary(geom = "errorbar", fun.data = mean_sdl, fun.args = list(mult = 1),
               width = .3) +
  geom_quasirandom(width = .2, size = 2, color = "grey20") +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  scale_x_discrete(labels = c("T0" = "basal", "T1" = "Clamp T1", "T2" = "Clamp T2"), name = NULL,
                   expand = expansion(add = .8)) +
  scale_fill_manual(values = c("steelblue2", "tomato")) +
  facet_wrap(~analyte, scales = "free") +
  
  theme.myClassic +
  theme(axis.title.y = element_text(margin = margin(r = 15)),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.spacing = unit(60, "pt")) 

ggsave(filename = "Clamp_glycemia_insulin.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 3.5, width = 6.5)

