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
    ggplot(aes(x = timepoint, y = Fcirc_atom.g.BW)) +
    stat_summary(geom = "bar", fun = mean, fill = "snow4", color = "black") +
    stat_summary(geom = "errorbar", fun.data = mean_sdl, fun.args = list(mult = 1),
                 width = .3) +
    geom_quasirandom(width = .1) +
    facet_grid(.~Compound, scales = "free", space = "free") +
    theme.myClassic +
    theme(axis.title.y = element_text(margin = margin(r = 15))) +
    scale_y_continuous(expand = expansion(mult = c(0, .1)), name = myTitle,
                       # limits = c(0, 2100)
    ) +
    scale_x_discrete(labels = c("T0" = "basal", "T1" = "Clamp T1", "T2" = "Clamp T2"), name = NULL,
                     expand = expansion(add = .8))
}

# endogenous production
p1 <- d.clamp.ALL %>% filter(Index == "Endo_Ra") %>% 
  func.plt.clamp(myTitle = "Endogenous Ra (nmol C/min/g)")
p1

# Total consumption
p2 <- d.clamp.ALL %>% 
  filter(Index == "Total Ra" & Compound == "Glucose") %>% 
  func.plt.clamp(myTitle = "Total Consumption (nmol C/min/g)")
p2

plot_grid(p1, p2, rel_widths = c(4, 2.2))
