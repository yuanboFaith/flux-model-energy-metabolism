# This script includes data analysis for peer review:
# 1. glucose and C16:0 flux under hyperinsulinemia condition
# 2. C18:1 Ra change under isoproterenol treatment
# 3. ITT, GTT
# 4. TAG Ra
# 5. flux comparison with the literature data
# 6. comparison with fluxes without non-negativity constraint

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
  func.plt.clamp(myTitle = "Rd (nmol C/min/g)") +
  scale_fill_manual(values = c("steelblue2", "tomato", "tomato"))
p3

plot_grid(p1, ggplot() + theme_void(),
          p2, ggplot() + theme_void(),
          p3, 
          rel_widths = c(2, .7, 2.8, .7, 2.8), nrow = 1)

ggsave(filename = "Clamp.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 3, width = 9)

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
       height = 2.5, width = 9)





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



# - <>-# - <>-# - <>-# - <>-# - <>-# - <>-# - <>-# - <>-# - <>-# - <>-# - <>-


d.ITT.GTT <- read_excel(
  "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/hyperinsulinemia.xlsx",
  sheet = "ITT_GTT")  

func.ITT.GTT <- function(whichTest = "ITT"){
  dg <- 5
  
  d.ITT.GTT.i <- d.ITT.GTT %>% filter(experiment == whichTest)
  
  plt.curve.ITT.GTT <- d.ITT.GTT.i %>%
    select(-AUC) %>% 
    pivot_longer(-c(animals, phenotype, BW, experiment), names_to = "time", values_to = "glycemia") %>% 
    mutate(time = as.numeric(time)) %>% 
    ggplot(aes(x = time, y = glycemia, color = phenotype)) +
    
    stat_summary(
      geom = "line", fun = mean, aes(group = phenotype), 
      position = position_dodge(dg), linewidth = .8) +
    stat_summary( 
      geom = "errorbar", fun.data = mean_sdl,
      fun.args = list(mult = 1), width = 13, 
      position = position_dodge(dg), linewidth = .8) +
    stat_summary(
      geom = "point", fun = mean, position = position_dodge(dg), 
      size = 3, shape = 21, fill = "white", stroke = .8) +
    
    # geom_point(shape = 21) +
    
    scale_y_continuous(limits = c(0, NA), breaks = seq(0, 600, 100)) +
    scale_x_continuous(breaks = seq(0, 120, 30)) +
    scale_color_manual(values = color.phenotype) +
    labs(y = "Glycemia (mg/dL)\n", x = "time (min)") +
    theme.myClassic +
    theme(legend.position = "none")
  #plt.curve.ITT.GTT
  
  
  plt.AUC <- d.ITT.GTT.i %>% 
    select(animals, phenotype, experiment, AUC) %>% 
    pivot_longer(-c(animals, phenotype, experiment), names_to = "time", values_to = "AUC") %>% 
    mutate(phenotype = factor(phenotype, levels = ordered.phenotype)) %>% 
    mutate(time = as.numeric(time)) %>% 
    ggplot(aes(x = phenotype, y = AUC, fill = phenotype)) + 
    stat_summary(geom = "bar", fun = mean, width = .8, color = "black", alpha = .7) +
    stat_summary(geom = "errorbar", fun.data = mean_sdl, 
                 fun.args = list(mult = 1),  width = .3) +
    geom_quasirandom() +
    scale_y_continuous(limits = c(0, NA), 
                       # breaks = seq(0, 50000, 10000),
                       n.breaks = 6,
                       expand = expansion(mult = c(0, .1))) +
    scale_x_discrete(expand = expansion(add = .8)) +
    scale_fill_manual(values = color.phenotype) +
    labs(y = "AUC (mg / dL × 120 min)\n", x = NULL) +
    theme.myClassic +
    theme(legend.position = "none",
          axis.text.x = element_blank())
  
  plot_grid(plt.curve.ITT.GTT, plt.AUC, align = "h")
  
}

GTT <- func.ITT.GTT(whichTest = "GTT" )
ITT <- func.ITT.GTT(whichTest = "ITT" )


plot_grid(GTT, ITT, align = "h")


ggsave(filename = "GTT_ITT.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 3, width = 12)


# t test
func.test.GTT.ITT <- function(whichTest){
  
  m <- d.ITT.GTT %>% 
    select(animals, phenotype, experiment, AUC) %>% 
    filter(experiment == whichTest) %>% aov(formula = AUC ~ phenotype)
  tukey_hsd(m)
}

func.test.GTT.ITT(whichTest = "ITT")
func.test.GTT.ITT(whichTest = "GTT")




# -<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-

setwd("/Users/boyuan/Desktop/Harvard/Research/db db mice/TAG decay kinetics/lipidomcis method/formal experiment")
path <- "TAG kinetics data_all_data.csv"
mydata <- read.csv(path) 

mydata <- mydata  %>% filter(isotopeLabel %in% c(
  "C12 PARENT", "C13-label-1", "C13-label-2", "C13-label-3", 
  "C13-label-16", "C13-label-17", "C13-label-18"))

d.corrected <- mydata %>% 
  group_by(isotopeLabel, compound, formula) %>% 
  summarise(across(M1.0:last_col(), sum)) %>% 
  ungroup() %>% 
  accucor::natural_abundance_correction(resolution = 120000)


d.tidy <- d.corrected$Normalized %>% 
  # tidy up
  pivot_longer(-c(Compound, C_Label), names_to = "sample", values_to = "fraction") %>% 
  separate(sample, into = c("mouse", "time"))  %>% 
  mutate(time = as.numeric(time)) %>% 
  
  # filter  
  # filter(! Compound %>% str_detect("20:1")) %>% 
  # filter(! Compound %>% str_detect("16:0_16:0"))  %>% 
  # filter(! Compound %>% str_detect("18:0")) %>% 
  filter(mouse %in% c("M1", "M6", "M8", "M9")) %>% 
  filter(time > 0) %>% 
  filter(time < 100) 

# first time point as normalization basis
d.normalizationBasis <- d.tidy %>% 
  group_by(Compound, mouse) %>% 
  filter(time == min(time) & C_Label == 16) %>% 
  rename(normBasis = fraction) %>% 
  arrange(mouse, time) %>% 
  select(-c(C_Label, time))

d.normalizationBasis


# combine tidy data with normalization basis
d.tidy.normalized <- d.tidy %>% 
  left_join(d.normalizationBasis) %>% 
  # add normalized fraction
  group_by(mouse, Compound) %>%
  mutate(fraction.norm = fraction / normBasis) 


d.kinetics.16 <- d.tidy.normalized %>% 
  filter(C_Label == 16) %>% 
  filter(time <= 60)

p <- d.kinetics.16 %>% 
  ggplot(aes(x = time, y = fraction.norm, color = Compound, shape = mouse))+
  ggbeeswarm::geom_quasirandom(width = 3, size = 3) 
# geom_line() +
# facet_wrap(~mouse, scales = "free") +
# scale_y_continuous(transform = "log2") + annotation_logticks(sides = "l") +
# geom_smooth(method = "lm", aes(group = 1), se = F, size = .5, color = "black") 
# facet_wrap(~mouse)

p

# Fit the non-linear regression model for the entire dataset
nls_model <- nls(fraction.norm ~ A * exp(-k * time), 
                 data = d.kinetics.16, 
                 start = list(A = 1, k = 0.05))

summary(nls_model)
# Extract the fitted parameters
params <- coef(nls_model)
print(params)

# Define the equation function using the fitted parameters
eqn <- function(x) params["A"] * exp(-params["k"] * x)

# Create the equation text to display on the plot
eqn_text <- paste0("y = ", round(params["A"], 2), " * e^{-", round(params["k"], 4), " * x}")


# fit regression for each TAG species

# Define a function to fit nls model for each group
fit_nls <- function(data) {
  nls(fraction.norm ~ A * exp(-k * time),
      data = data, 
      start = list(A = 1, k = 0.05))
}

nls_fits <- d.kinetics.16 %>%
  group_by(Compound) %>%
  do(model = fit_nls(.))

# Predict values based on the model
predictions <- nls_fits %>%
  rowwise() %>%
  do(data.frame(time = seq(15, 60, 1), 
                predicted = predict(.$model, newdata = data.frame(time = seq(15, 60, 1))),
                Compound = .$Compound))



# Plot the data and the fitted regression line

p + 
  # fitted line for each TAG species
  geom_line(data = predictions, 
            aes(x = time, y = predicted, color = Compound), inherit.aes = F,
            size = .5, show.legend = F) +
  
  stat_function(fun = eqn, color = "black", linewidth = 1.2) +
  # geom_smooth(method = "lm", se = F, linewidth = 1, color = "black", aes(group = 1)) +
  # theme_minimal(base_size = 19) +
  theme.myClassic +
  # Add the equation text on the plot
  # annotate("text", 
  #          x = 30, y = .4,
  #          # x = max(d.kinetics.16$time) * 0.7, 
  #          # y = max(d.kinetics.16$fraction) * 0.9, 
  #          label = eqn_text, size = 5, color = "blue") +
  
  # log scale with linear regression
  # scale_y_continuous(transform = "log2")  + annotation_logticks(sides = "l") +
  
  scale_x_continuous(breaks = c(15, 30, 45, 60)) +
  coord_cartesian(ylim = c(.2, 1)) +
  labs(x = "time (min)", y = "normalized labeling\n") +
  scale_color_brewer(palette = "Set3") +
  guides(shape = "none") +
  theme(
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 15), 
    legend.title = element_blank()
  )

ggsave(
  filename = "TAG decay curve.pdf",
  path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures", 
  width = 7, height = 5)


halfLife <- log(2)/params[["k"]]
halfLife

# the physiological meaning of K is fractional turnover rate
# Ra = k * C * V

# 1.5 mM TAG, or 4.5 mM equivalent TAG-FA; mM = nmol/uL
# volumne
V <-  30 * 49 # 49 uL plasma/g BW

params[["k"]] * V * (.35 * 3) / 25 # nmol/g/min

# x <- 1:80
# y = 1.02 * exp(-.0124 * x)
# plot(x, y)





# 5. flux comparison with the literature data

d <- read_excel(
  "/Users/boyuan/Desktop/Harvard/Research/db db mice/Infusion data/flux comparision.xlsx")

d %>% 
  mutate(Pathway = fct_reorder(Pathway, `nmol/min/g`)) %>% 
  ggplot(aes(x = Pathway, y = `nmol/min/g`, color = reference, fill = reference)) +
  geom_point(size = 5, shape = 21) +
  coord_flip() +
  scale_y_continuous(trans = "log10") +
  annotation_logticks(sides = "b") +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c(rep("white", 4), "black")) +
  theme(axis.text = element_text(color = "black")) +
  labs(x = NULL) +
  scale_color_brewer(palette = "Set1")




# 6. Compare with fluxes computed without non-negativity constraint


func.calculateInterConvertFluxes <- function(constraint = 0){
  
  d.directFlux_interConversion = tibble()
  
  matrixDimension = d.enrich.atom.summary$Compound %>% n_distinct()
  
  for (targetCompound in d.enrich.atom.summary$Compound %>% unique()  ){
    for (phenotypeSelected in c("WT", "ob/ob", "HFD")){
      
      # targetCompound = "C16:0"
      # phenotypeSelected = "HFD"
      
      # subset of selected phenotype
      d.enrich.atom.summary_Pheno.i = d.enrich.atom.summary %>%
        filter(phenotype == phenotypeSelected) %>% ungroup()
      
      
      # labeling of target metabolite (TM) under TM as tracer
      L.TM.TM = d.enrich.atom.summary_Pheno.i %>% 
        filter(infused.tracer == targetCompound & Compound == targetCompound)
      # one row matrix: mean of TM tracer, TM labeling
      mat_TM.TM_mean =  L.TM.TM$enrich.original.mean %>% 
        rep(time = matrixDimension) %>% matrix(nrow = 1)
      # one row matrix: standard error of the mean of TM tracer, TM labeling
      mat_TM.TM_sem =  L.TM.TM$enrich.original.sem %>% 
        rep(time = matrixDimension) %>% matrix(nrow = 1)
      
      
      # labeling of TM under Source Metabolites (S) as tracer (matrix naming as L . tracer . labeled_metabolite)
      L.S.TM = d.enrich.atom.summary_Pheno.i %>% 
        filter(infused.tracer != targetCompound & Compound == targetCompound)
      # mean matrix, TM labeling under Sources infusion
      mat_L.S.TM_mean =  L.S.TM$enrich.original.mean %>% 
        rep(time = matrixDimension) %>%  matrix(nrow = matrixDimension-1)
      # sem matrix, TM labeling under Sources infusion
      mat_L.S.TM_sem =  L.S.TM$enrich.original.sem %>% rep(time = matrixDimension) %>% 
        matrix(nrow = matrixDimension-1) 
      
      # matrix L1: labeling of TM with infusion of TM and S as tracer, combining above matrices
      mat_L1.mean = rbind(mat_TM.TM_mean, mat_L.S.TM_mean)
      rownames(mat_L1.mean) = c(targetCompound, L.S.TM$infused.tracer %>% as.character())
      
      mat_L1.sem = rbind(mat_TM.TM_sem, mat_L.S.TM_sem)
      rownames(mat_L1.sem) = rownames(mat_L1.mean)
      
      
      # --- Create matrix L2
      
      # labeling of source metabolites under TM as tracer
      L.TM.S = d.enrich.atom.summary_Pheno.i %>% 
        filter(infused.tracer == targetCompound & Compound != targetCompound) 
      # a) mean matrix
      mat_L.TM.S_mean = L.TM.S %>% select(infused.tracer, Compound, enrich.original.mean) %>% 
        spread(key = Compound, value = enrich.original.mean)
      # b) SEM matrix
      mat_L.TM.S_sem = L.TM.S %>% select(infused.tracer, Compound, enrich.original.sem) %>% 
        spread(key = Compound, value = enrich.original.sem)
      
      # labeling of source metabolites under S as tracer
      L.S.S = d.enrich.atom.summary_Pheno.i %>% 
        filter(! (infused.tracer == targetCompound | Compound == targetCompound))
      # a) mean matrix
      mat_L.S.S_mean = L.S.S %>% select(c(infused.tracer, Compound, enrich.original.mean)) %>% 
        spread(key = Compound, value = enrich.original.mean)
      # b) SEM matrix
      mat_L.S.S_sem = L.S.S %>% select(c(infused.tracer, Compound, enrich.original.sem)) %>% 
        spread(key = Compound, value = enrich.original.sem)
      
      # Combine: 
      # matrix of mean
      mat_L2_mean =  rbind(mat_L.TM.S_mean, mat_L.S.S_mean) %>% 
        # add zero column as the last column
        cbind(storage = rep(0, time = matrixDimension))
      # add row name; remove tracer column
      rownames(mat_L2_mean) = mat_L2_mean$infused.tracer
      mat_L2_mean = mat_L2_mean %>% select(-1)
      
      # matrix of SEM
      mat_L2_sem =  rbind(mat_L.TM.S_sem, mat_L.S.S_sem) %>% 
        # add zero column as the last column
        cbind(storage = rep(0, time = matrixDimension))
      # add row name; remove tracer column
      rownames(mat_L2_sem) = mat_L2_sem$infused.tracer
      mat_L2_sem = mat_L2_sem %>% select(-1)
      
      # Infusion rate 
      R.TM.nmol.g.min = (d.normalized.tidy %>% filter(infused.tracer == targetCompound & phenotype == phenotypeSelected))$infusion_nmol13C.atoms_perMin.gBW %>% unique()
      
      # labeling of TM under infusion of TM
      L.TM.TM.mean = L.TM.TM$enrich.original.mean 
      L.TM.TM.sem = L.TM.TM$enrich.original.sem
      
      
      # Monte Carlo simulation
      d.directFluxes.monteCarloRecorder = c()
      for (i in 1:100){
        # standard Gaussian distribution
        mat.standardGaussain = rnorm(n = matrixDimension * matrixDimension, mean = 0, sd = 1) %>% 
          matrix(nrow = matrixDimension)
        
        # L1 and L2 matrix on the left side
        mat_L1.mean.MonteCarlo = mat_L1.mean + mat_L1.sem * mat.standardGaussain
        mat_L2_mean.MonteCarlo = mat_L2_mean + mat_L2_sem * mat.standardGaussain
        
        # The right side r vector
        L.TM.TM.monteCarlo = L.TM.TM.mean + rnorm(n = 1, mean = 0, sd = 1) * L.TM.TM.sem 
        v.R.mean.monteCarlo = c(R.TM.nmol.g.min * (1 -  L.TM.TM.monteCarlo),  # the first element
                                rep(0, matrixDimension-1)) # the remaining zeros 
        
        
        # create matrix of 1 
        mat_1 <- matrix(1, nrow = matrixDimension, ncol = matrixDimension)
        
        # solve matrix operation 
        # WITH non-negative constraints
        list.directFlx.i = lsei(
          A = (mat_L1.mean.MonteCarlo - mat_L2_mean.MonteCarlo) / (mat_1 - mat_L2_mean.MonteCarlo) ,
          B = v.R.mean.monteCarlo, # this solves Ax = B
          G = diag(1, nrow = matrixDimension), # G is the identity matrix, dimension = B-vector length
          H = rep(constraint, matrixDimension), # vector of zero, length = B-vector length
          verbose = FALSE)
        v.F.nmol.min.g = list.directFlx.i$X
        
        # record output for each iteration
        d.directFluxes.monteCarloRecorder = d.directFluxes.monteCarloRecorder %>% cbind(v.F.nmol.min.g)
      }
      
      # plot simulated each iteration
      # d.directFluxes.monteCarloRecorder %>% Heatmap(cluster_columns = F, cluster_rows = F)
      
      # summarize
      d.directFlux_pheno.i = tibble(nmol.min.g.mean = d.directFluxes.monteCarloRecorder %>% apply(MARGIN = 1, FUN = mean),
                                    # calculate the STD, but derived from simulated SEM monte carlo matrix, thus name as SEM
                                    nmol.min.g.SEM = d.directFluxes.monteCarloRecorder %>% apply(MARGIN = 1, FUN = sd)) %>% 
        mutate(phenotype = phenotypeSelected) %>% 
        mutate(targetCompound = targetCompound) %>% 
        mutate(sources = c(L.S.TM$infused.tracer %>% as.character(), "storage")) %>% 
        select(phenotype, targetCompound, sources, contains("nmol.min"))
      
      # this is nmol carbon atoms / min / g BW
      d.directFlux_interConversion = rbind(d.directFlux_interConversion, d.directFlux_pheno.i)  
    }
  }
  
  d.directFlux_interConversion = d.directFlux_interConversion %>% 
    mutate(targetCompound = factor(targetCompound, levels = c(ordered.Compound), ordered = F),
           sources = factor(sources, levels = ordered.Compound %>% rev(), ordered = F)) %>% 
    arrange(sources) %>% 
    group_by(phenotype, targetCompound) %>% 
    mutate(nmol.min.g.SEM.Y.axis.position = (func.cumulatedSum(abs(nmol.min.g.mean)) + func.cumulatedSum(nmol.min.g.mean))/2 ) 
  
  # Calculate per animal fluxes
  d.directFlux_interConversion = d.directFlux_interConversion %>%  # good
    # combine with body weight data
    left_join(d.BW %>% rename(targetCompound = Compounds), 
              by = c("targetCompound", "phenotype")) %>% 
    # calculate animal-based flux
    mutate(nmol.min.animal.mean = nmol.min.g.mean * BW.mean,
           nmol.min.animal.SEM =  nmol.min.g.SEM * BW.mean,
           nmol.min.animal.SEM.Y.axis.position = func.cumulatedSum(abs(nmol.min.animal.mean))) %>% 
    mutate(targetCompound = factor(targetCompound, levels = ordered.Compound),
           phenotype = factor(phenotype, levels = ordered.phenotype))
  
  d.directFlux_interConversion <- d.directFlux_interConversion %>% 
    filter(phenotype == "WT") %>% 
    ungroup() %>% 
    select(-phenotype, targetCompound, sources, contains("animal")) %>% 
    mutate(constraint = as.character(constraint), .before = 1)
  
  d.directFlux_interConversion
  
}

d.withConstraint <- func.calculateInterConvertFluxes(constraint = 0) 
d.withConstraint

d.noConstraint <- func.calculateInterConvertFluxes(constraint = -100000)
d.noConstraint

dwnC <- bind_rows(d.withConstraint, d.noConstraint)

# per g BW
dwnC %>% 
  # filter(targetCompound %in% c("Glucose", "Lactate", "Glycerol", "Glutamine", "3-HB")) %>% 
  # filter(! targetCompound %in% c("Alanine")) %>% 
  ggplot(aes(x = constraint, y = nmol.min.g.mean, fill = sources)) +
  geom_bar(stat = "identity", color = "black", size = .4, alpha = .9) +
  geom_errorbar(aes(ymin = nmol.min.g.SEM.Y.axis.position - nmol.min.g.SEM,
                    ymax = nmol.min.g.SEM.Y.axis.position),
                width = .3, size = .5) +
  facet_wrap(~targetCompound, scales = "free", nrow = 1) +
  # scale_color_manual(values = color.Compounds) +
  scale_fill_manual(values = color.Compounds, guide = guide_legend(reverse=T)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), n.breaks = 7) +
  expand_limits(y = 0) +
  scale_x_discrete(expand = expansion(mult = c(.8, .8)),
                   labels = c("w/o", "w/")) +
  geom_hline(yintercept = 0, linewidth = 1) +
  theme.myClassic +
  theme(axis.title.x = element_blank(),
        strip.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 17),
        axis.line.x.bottom = element_blank(),
        # axis.text.x = element_text(size = 10),
        panel.spacing = unit(20, "pt"),
        legend.title = element_text(color = "black", face = "bold", size = 13),
        # legend.position = c(.93, .23)
  ) +
  guides(fill = guide_legend(reverse = F)) +
  labs(y = "nmol C / min / g\n", fill = "sources")+
  coord_cartesian(clip = "off")

ggsave(filename = "flux comparison constraints.pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 3.5, width = 20)



save.image(file = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/13_validation_experiment.RData")
