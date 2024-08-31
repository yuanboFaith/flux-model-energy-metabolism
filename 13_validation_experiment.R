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

d.corrected.tidy.23 %>% filter(C_Label == 0 & Compound == "C16:0") %>% view()

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


func.plt.clamp <- function(data, myTitle, pt.size = 2){
  data %>% 
    ggplot(aes(x = timepoint, y = Fcirc_atom.g.BW, fill = timepoint)) +
    stat_summary(geom = "bar", fun = mean, color = "black", alpha = .7) +
    stat_summary(geom = "errorbar", fun.data = mean_sdl, fun.args = list(mult = 1),
                 width = .3) +
    geom_quasirandom(width = .2, size = pt.size, color = "grey20") +
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
  filter(Index %in% "Total Ra" & Compound == "Glucose" |
           Index == "Endo_Ra" & Compound == "Glucose" & treatment == "basal") %>% 
  func.plt.clamp(myTitle = "Rd (nmol C/min/g)",
                  pt.size = c(rep(1, 17), rep(2, 10))) +
  scale_fill_manual(values = c("steelblue2", "tomato", "tomato"))
p3

plot_grid(p1, p2, p3, rel_widths = c(2, 2.8, 2.8), nrow = 1)

ggsave(filename = "Clamp.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 4, width = 9)

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
       height = 3.5, width = 9)





# -<>-# -<>-# -<>-# -<>-# -<>-# -<>-# -<>-# -<>-# -<>-# -<>-# -<>-# -<>-# -<>-
path.iso <- "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/isoproterenol.xlsx"
d.iso.norm <- read_excel(path.iso, sheet = "Normalized")  

# Fcirc
d.iso.norm.summary <- d.iso.norm %>% select(C_Label, contains("D505")) %>% 
  pivot_longer(-C_Label, names_to = "sample", values_to = "enrich") %>% 
  # average labeling
  mutate(enrich.weighted = enrich * C_Label / 18) %>% 
  group_by(sample) %>% 
  summarise(enrich = sum(enrich.weighted)) %>% 
  # calculate Fcirc: nmol C/min/g, 10 uL/min/animal infusion, 11.8 mM concentration of 
  mutate(Fcirc.anim = (1 - enrich) / enrich * 11.8 * 10 * 18,
         Fcirc.anim.gBW = Fcirc.anim / 29) %>% 
  # extract sample info
  separate(sample, into = c("Diet", "Mouse", "treatment", "run")) %>% 
  mutate(ID = str_c(Diet, ".", Mouse, ".", run))



# TIC
d.iso.TIC <- read_excel(path.iso, sheet = "Corrected")  
d.iso.TIC.summary <-  d.iso.TIC %>% 
  filter(C_Label == 0) %>% 
  select(contains("D505")) %>% 
  pivot_longer(everything(), names_to = "sample", values_to = "TIC") %>% 
  
  # extract sample info
  separate(sample, into = c("Diet", "Mouse", "treatment", "run")) %>% 
  mutate(ID = str_c(Diet, ".", Mouse, ".", run))

# Combine the dataset
d.iso.all <- d.iso.norm.summary %>% left_join(d.iso.TIC.summary) %>% 
  # remove outlier, the changing trend is opposite in both TIC and labeling
  filter(! ID == "D5053.M1.run2") %>% 
  arrange(ID) %>% 
  # calculate delta change
  group_by(ID) %>% 
  mutate(delta.Fcirc = (max(Fcirc.anim.gBW) - min(Fcirc.anim.gBW))/min(Fcirc.anim.gBW) * 100,
         delta.TIC = (max(TIC) - min(TIC))/min(TIC) * 100) %>% 
  # remove 3 observations where there is only mild change in Fcirc
  filter(delta.Fcirc > 10)

# plot
f.plt.iso <- function(myY){
  d.iso.all %>% 
    ggplot(aes(x = treatment, y = {{myY}})) +
    geom_point(size = 3) +
    geom_line(aes(group = ID)) +
    coord_cartesian(ylim = c(0, NA)) +
    theme.myClassic 
  # theme(legend.position = "none")
}

p.iso.1 <- f.plt.iso(myY = TIC) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)),
                     labels = function(x)(x/10^8),
                     name = "TIC ( x 10^8)\n") 

p.iso.2 <- f.plt.iso(myY = Fcirc.anim.gBW) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)),
                     name = "nmol C / min / g \n") 


# calculate relative change 
p.iso.3 <- d.iso.all %>% 
  group_by(ID) %>% 
  summarise(delta.Fcirc = (max(Fcirc.anim.gBW) - min(Fcirc.anim.gBW))/min(Fcirc.anim.gBW) * 100,
            delta.TIC = (max(TIC) - min(TIC))/min(TIC) * 100,
            basalTIC = min(TIC),
            basalFcirc = min(Fcirc.anim.gBW)) %>% 
  
  ggplot(aes(x = delta.Fcirc, y = delta.TIC)) +
  geom_point(color = "black", size = 3) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
  geom_abline(slope = 1, linetype = "dashed", color = "grey") +
  theme.myClassic

plot_grid(p.iso.1, p.iso.2, p.iso.3, nrow = 1)

ggsave(filename = "isoproterenol.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 3, width = 9)

d.Ra <- d.iso.all %>% select(treatment, ID, Fcirc.anim.gBW) %>% spread(treatment,  Fcirc.anim.gBW)
t.test(x = d.Ra$basal, y = d.Ra$iso, paired = T)  

d.TIC <- d.iso.all %>% select(treatment, ID, TIC) %>% spread(treatment,  TIC)
t.test(x = d.TIC$basal, y = d.TIC$iso, paired = T)  
