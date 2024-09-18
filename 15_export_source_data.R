rm(list = ls())

library(plyr)
library(rebus)
library(viridis)
library(lubridate)
library(readxl)
library(purrr)
library(broom)
library(RColorBrewer)
library(splines)
library(pryr)
library(cowplot)
library(gridExtra)
library(ggrepel)
library(ggbeeswarm)
library(scales)
library(ggsci)
library(tidyverse)
library(xlsx)

# export data to a single spreadsheet, with specified sheet name

ff <- function(data, sheet, append = T){ 
  data %>% write.xlsx(
    file = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/Data S1 - Source Data/Data S1.xlsx", 
    sheetName = sheet, 
    row.names = F, 
    append = append)} # append as the default, except for the first sheet


# Figure 1
load("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/2_core_13CO2_bolus_injection.RData")


# Fig. 1B (middle), S1B (up) 
d.all.treated.bolus.injection %>% 
  filter(phenotype == "WT" ) %>%
  ungroup() %>% 
  select(tracer, time.h, round, cage, P1.umol.min) %>% 
  mutate(mouse_ID = paste0(round, "-", cage), .keep = "unused") %>% 
  mutate(`time (h)` = time.h,
         `13CO2 exhaled (nmol/min/animal)` = P1.umol.min * 1000, 
         .keep = "unused") %>% 
  arrange(tracer, mouse_ID, `time (h)`) %>% 
  as.data.frame() %>% 
  ff(sheet = "1B (middle), S1B (up)", append = F) # first sheet

# Fig. 1B (right), S1B (bottom)
d.all.treated.recovery.timeSections.clean %>% 
  filter(phenotype == "WT") %>% 
  ungroup() %>% 
  select(tracer, time.h, round, cage, recovery.timeSection) %>% 
  mutate(mouse_ID = paste0(round, "-", cage), .keep = "unused") %>% 
  mutate(`time (h)` = time.h,
         `recovery` = recovery.timeSection, 
         .keep = "unused") %>% 
  arrange(tracer, mouse_ID, `time (h)`) %>% 
  as.data.frame() %>% 
  ff(sheet = "1B (right), S1B (bottom)")



load("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/3_core_13CO2_infusion_physiology_basics.RData")


# 1C
d.14C[1:3, ] %>% 
  select(mouse, `DPM measured`, inj.DPM, carc.recovery) %>% 
  rename(`DPM measured in carcass` = `DPM measured`,
         `DPM administered` = inj.DPM,
         `recovery (storage fraction)` = carc.recovery) %>% 
  as.data.frame() %>% 
  ff(sheet = "1C")



# Fig. 1D, S1C
d.infusion.expCurve %>% 
  filter(phenotype == "WT") %>% 
  ungroup() %>% 
  select(tracer, time.h, round, cage, P1.umol.min.normalized) %>% 
  mutate(mouse_ID = paste0(round, "-", cage), .keep = "unused") %>% 
  mutate(`time (h)` = time.h,
         `13CO2 exhaled µmol/min/animal` = P1.umol.min.normalized, 
         .keep = "unused") %>% 
  arrange(tracer, mouse_ID, `time (h)`) %>% 
  as.data.frame() %>% 
  ff(sheet = "1D, S1C")


# 1E, S1E
f.clean %>% select(-phenotype) %>% arrange(tracer, inj) %>% 
  select(tracer, inj, recovery, fox) %>% 
  rename(`administration method` = inj) %>% 
  as.data.frame() %>% ff("1E, S1E")


# S1D
A %>% 
  filter(phenotype == "WT") %>% 
  rename(`bolus recovery` = bolus,
         `bolus recovery standard deviation` = bolus.SD,
         `infusion recovery` = infusion,
         `infusion recovery standard deviation` = infusion.SD) %>%
  ungroup() %>% select(-phenotype) %>% 
  as.data.frame() %>% ff("S1D")  


# S1F
load("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/13_validation_experiment.RData")

d.kinetics.16 %>% select(mouse, Compound, time, fraction.norm) %>% 
  arrange(mouse, Compound, time) %>% 
  rename("M+16 normalized fraction" = "fraction.norm",
         `time (min)` = time,
         mouse_ID = mouse) %>% 
  as.data.frame() %>% ff("S1F")  


# 2C, S2C
load("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/6_core_consumption_flux.RData")

d.ox.sink.overal.2 %>% 
  filter(phenotype == "WT") %>% 
  ungroup() %>% select(-c(phenotype, yAxis.error, phenotype)) %>% 
  mutate(destiny = ifelse(destiny == "sink.overal", "storage", "oxidation")) %>% 
  rename("nmol C/min/animal"     = nmol.min.animal,
         "nmol C/min/animal SEM" = nmol.min.animal.SEM,
         "overall destiny"       = destiny) %>% 
  arrange(Compound, `overall destiny`) %>% 
  
  # per g BW
  mutate(`nmol C/min/g` = `nmol C/min/animal` / 29,
         `nmol C/min/g SEM` = `nmol C/min/animal SEM` / 29) %>% 
  as.data.frame() %>% ff("2C, S2C")  



# S2(A2)
x1 <- d.clamp.ALL %>% filter(Compound == "C16:0") %>% 
  filter(Index == "Endo_Ra")

x2 <- d.clamp.ALL %>% filter(Compound == "Glucose") 

bind_rows(x1, x2) %>% 
  mutate(Index = str_replace(Index, "Endo_", "Endogenous ") %>% 
           str_replace("Total Ra", "Total Ra (Endogenous+infusate; Rd)")) %>% 
  select(Compound, treatment, mouse_ID, timepoint, Index, Fcirc_atom.g.BW) %>% 
  rename(flux = Index, 
         `nmol C/min/g` = Fcirc_atom.g.BW) %>% 
  as.data.frame() %>% ff("S2(A2)")  



# S2(A3)
d.gly.ins %>% rename(`concentration` = conc.) %>% 
  select(analyte, time, Mouse, concentration, unit) %>% 
  arrange(analyte, time, Mouse) %>% 
  mutate(unit = ifelse(unit == "-", "relative MS intensity", unit)) %>% 
  as.data.frame() %>% ff("S2(A3)")  


# S2(B2)
d.iso.all %>% ungroup() %>%
  mutate(mouse_ID = rep(1:5, each = 2), .before = 1) %>% 
  select(mouse_ID, treatment, Fcirc.anim, Fcirc.anim.gBW, TIC) %>% 
  mutate(treatment = ifelse(treatment == "iso", "isoproterenol", "basal")) %>% 
  rename(`C nmol/min/animal` = Fcirc.anim,
         `C nmol/min/g` = Fcirc.anim.gBW,
         `relative intensity` = TIC) %>% 
  as.data.frame() %>% ff("S2(B2)")  



# 3B
d.enrich.atom.summary %>% 
  filter(phenotype == "WT") %>% 
  ungroup() %>% 
  select(infused.tracer, Compound, contains("enrich"), n.replicate) %>% 
  rename(`infused tracer` = infused.tracer,
         `labeled metabolite` = Compound,
         `averaged C atom labeling` =  enrich.original.mean,
         `averaged C atom labeling standard deviation` =  enrich.original.sd,
         `averaged C atom labeling SEM` =  enrich.original.sem,
         `replicate N` = n.replicate) %>% 
  as.data.frame() %>% ff("3B")  



# 3C(left), 5C(left), S5D
deE <- d.energyExpenditure %>% 
  ungroup() %>% 
  select(phenotype, round.cage, contains("CO2"), contains("O2."), RER, contains("cal"), BW) %>% 
  arrange(phenotype) %>% 
  rename(mouse_ID = round.cage) 
colnames(deE) <- c("phenotype", "mouse ID", "CO2 umol/min", "CO2 L/min", "O2 umol/min", "O2 L/min", "RER", "kcal/h", "cal/min", "body weight")

deE %>% as.data.frame() %>% ff("3C(left), 5C(left), S5D")  


# 3C (right), 5C (right)
fCn <- d.flx.CO2.nonOx.sink %>% 
  select(phenotype, Compounds, destiny, contains(".g"), contains("animal"), -y.error.animal)
colnames(fCn) <- c(
  "phenotype", "nutrient", "destiny (direct)", "nmol C/min/g", "nmol/min/g SEM", "nmol/min/animal", "nmol/min/animal SEM")
fCn %>% as.data.frame() %>% ff("3C (right), 5C (right)")  



# 3D, 5D
load("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/8_core_flux extrapolation.RData")

d.flx.CO2.nonOx.sink.categorized.summary %>% 
  mutate(destiny = ifelse(destiny == "CO2", "oxidation", "storage")) %>% 
  select(phenotype, category, destiny, contains("nmol")) %>% 
  arrange(phenotype, category, destiny) %>% 
  rename(`nmol C/min/animal` = categoryTotal.nmol.min.animal.mean,
         `nmol C/min/animal SEM` = categoryTotal.nmol.min.animal.SEM) %>% 
  as.data.frame() %>% ff("3D, 5D")  



# 4B, S4A
dFiC <- d.directFlux_interConversion %>% 
  filter(phenotype == "WT") %>% ungroup() %>% 
  select(targetCompound, sources, 
         nmol.min.g.mean, nmol.min.g.SEM, nmol.min.animal.mean, nmol.min.animal.SEM)

colnames(dFiC) <- c("nutrient", "sources", "nmol/min/g", "nmol/min/g SEM", "nmol/min/animal", "nmol/min/animal SEM")

dFiC %>% as.data.frame() %>% ff("4B, S4A")  


# 4C, S4B
dtP <- d.destiny.totalProduction %>% 
  filter(phenotype == "WT") %>% ungroup() %>% 
  select(Compounds, destiny, 
         nmol.min.g.mean, nmol.min.g.SEM, nmol.min.animal.mean, nmol.min.animal.SEM) %>% 
  arrange(Compounds, destiny)

colnames(dtP) <-   c("nutrient", "destiny", "nmol/min/g", "nmol/min/g SEM", "nmol/min/animal", "nmol/min/animal SEM")

dtP %>% as.data.frame() %>% ff("4C, S4B")


# S4C
d.ultimate.S.frac %>% 
  select(phenotype, compound, store, Ult.frac.mean, Ult.frac.SEM) %>% 
  arrange(phenotype, compound, store) %>% 
  rename(`ultimate contribution fraction` = Ult.frac.mean,
         `ultimate contribution fraction SEM` = Ult.frac.SEM) %>% 
  as.data.frame() %>% ff("S4C")


# S4D, S8B
load("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/9_core_ATP_estimate.RData")

d.ATP.tidy4 %>% 
  select(phenotype, cycle, wasted, wasted.SEM) %>% 
  rename(`ATP wasted (umol/min/animal)` = wasted,
         `ATP wasted SEM (umol/min/animal)` = wasted.SEM) %>% 
  arrange(phenotype, cycle) %>% 
  as.data.frame() %>% ff("S4D, S8B")



# 5A
load("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/3_core_13CO2_infusion_physiology_basics.RData")

d.fox.infuse.bolus.subset %>% 
  relocate(phenotype, tracer, inj) %>% 
  arrange(phenotype, tracer, inj) %>% 
  rename(nutrient = tracer,
         `adminitration` = inj,
         `fraction of oxidation` = fox) %>% 
  mutate(adminitration = ifelse(adminitration == "bolus inj.", "bolus injection", adminitration)) %>% 
  as.data.frame() %>% ff("5A")



# 5B, S6B
pcmn <-  d.production.consumption.mirror.normalized %>% 
  relocate(phenotype, about, index, from, to) %>% 
  # rename(`nutrients (subplots)` = about,
  #        `panel (production or sinks)` = index) %>% 
  # check data structure
  # filter(`nutrients (subplots)` == "Glucose" & from %in% c("Glucose", "Lactate") & to %in% c("Glucose", "Lactate")  & phenotype == "WT")
  
  # convert lower half part of y-axis to positive values
  mutate(across(where(is.numeric), ~abs(.))) %>% 
  ungroup() %>% 
  select(-c(x.axis, myfill, contains("err.Y"), BW, fat, lean)) %>% 
  mutate(phenotype = factor(phenotype, levels = ordered.phenotype)) %>% 
  arrange(phenotype, about, index, from, to)

colnames(pcmn) <- c(
  "phenotype", "nutrients (subplots)", "panel (production or sinks)", "from", "to",
  "nmol/min/animal", "nmol/min/animal SEM", 
  "nmol/min/g body mass", "nmol/min/g body mass SEM", 
  "nmol/min g fat mass", "nmol/min/g fat mass SEM", 
  "nmol/min/g lean mass", "nmol/min/g lean mass SEM")

pcmn %>% as.data.frame() %>% ff("4D, 5B, S6B, 6")




# S5A
d.BW2 %>% group_by(mouse_ID, phenotype) %>% summarise(BW = mean(BW)) %>% 
  relocate(phenotype) %>% 
  rename(`body weight (g)` = BW) %>% 
  as.data.frame() %>% ff("S5A")


# S5B
d.glycemia %>% 
  select(phenotype, round, metCage, glycemia) %>% 
  filter(! is.na(glycemia)) %>% 
  mutate(`ID` = str_c(round, "-", metCage), 
         .keep = "unused", .after = phenotype) %>% 
  as.data.frame() %>% ff("S5B")


# S5C
d.insulin %>% select(phenotype, samples, conc_corrected) %>% 
  filter(! is.na(conc_corrected)) %>%
  arrange(phenotype) %>%
  mutate(ID = samples, .keep = "unused", .after = 1) %>% 
  as.data.frame() %>% ff("S5C")


# S5E
d.intensity.norm.eachMouse %>% 
  arrange(phenotype, Compound) %>% 
  rename(`relative abundance` = intens.norm.mean) %>% 
  relocate(phenotype, Compound) %>% 
  as.data.frame() %>% ff("S5E")  

# S5(F-G)
d.ITT.GTT %>% select(-BW) %>% 
  relocate(experiment, phenotype) %>% 
  rename(`0 min`   = `0`,
         `15 min`  = `15`,
         `30 min`  = `30`,
         `60 min`  = `60`,
         `120 min` = `120`,
         `area under curve` = AUC) %>% 
  as.data.frame() %>% ff("S5(F-G)")  



# S6A
d.enrich.atom.selected %>% filter(phenotype != "db/db") %>% 
  select(phenotype, infused.tracer, Compound, infusion_mouseID, enrich.atom.corrected) %>% 
  arrange(phenotype, infused.tracer, Compound) %>% 
  rename(`infused tracer` = infused.tracer,
         `labeled metabolite` = Compound,
         `mouse ID` = infusion_mouseID,
         `averaged C atom labeling` =  enrich.atom.corrected) %>% 
  as.data.frame() %>% ff("S3A,S6A")  


# S7A
load("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/8_core_flux_extrapolation.RData")

d.volcano %>% 
  relocate(compare, from, to, flux, ratio, pvalue, p.adj, contains("neg")) %>% 
  arrange(compare, from, to) %>% 
  mutate(pathway = paste(from, "→", to), .after = 1, .keep = "unused") %>% 
  rename(comparison = compare, 
         `average flux umol/min/animal` = flux,
         `fold change` = ratio,
         p = pvalue,
         `adjusted p` = p.adj,
         `-log(p)` = neg.log.p.adj,
         `-log(adjusted p)` = neg.log.pvalue) %>% 
  mutate(comparison = ifelse(comparison == "HFD_WT", "HFD vs. WT", "ob/ob vs. WT")) %>% 
  as.data.frame() %>% ff("S7A")  


# S7B
d.Fcirc.BW %>% 
  select(phenotype, infusion_mouseID,  Compound, Fcirc_animal, BW) %>% 
  rename(`mouse ID` = infusion_mouseID,
         `Ra (µmol molecules/min/animal)` = Fcirc_animal, 
         `body weight g` = BW) %>% 
  arrange(phenotype, Compound) %>% 
  as.data.frame() %>% ff("S7B")  

# S7C
d.Fcirc.BW.matched %>% 
  select(phenotype, Compound, infusion_mouseID, BW.range, Fcirc_animal) %>% 
  rename(`mouse ID` = infusion_mouseID,
         `Ra (µmol molecules/min/animal)` = Fcirc_animal,
         `body weight range (g)` = BW.range) %>% 
  arrange(phenotype, Compound) %>% 
  as.data.frame() %>% ff("S7C")  


# S8D
dwnC %>% 
  select(constraint, targetCompound, sources, nmol.min.g.mean, nmol.min.g.SEM) %>% 
  rename(nutrient = targetCompound,
         `nmol/min/g` = nmol.min.g.mean,
         `nmol/min/g SEM` = nmol.min.g.SEM) %>% 
  mutate(constraint = ifelse(constraint == "0", "with constraint", "without constraint")) %>% 
  as.data.frame() %>% ff("S8D")  


