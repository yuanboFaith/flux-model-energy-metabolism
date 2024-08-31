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

load("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/8_core_flux extrapolation.RData")


ordered.Compound.numbered <- 1:length(ordered.Compound)
names(ordered.Compound.numbered) <- ordered.Compound


d.paired <- d.from.to %>% 
  # rename non-Ox sink as storage
  mutate(from = str_replace(from, "non-Ox sink", "storage"),
         to = str_replace(to, "non-Ox sink", "storage")) %>% 
  # recode nutrients with numbers
  mutate(from2 = ordered.Compound.numbered[as.character(from)],
         to2 = ordered.Compound.numbered[as.character(to)]) %>% 
  # reorder the `from` and `to` nutrients 
  mutate(index = 1:nrow(d.from.to)) %>% 
  group_by(index) %>% 
  mutate(from3 = min(from2, to2),
         to3 = max(from2, to2)) %>% 
  # cycle name
  mutate(A = ordered.Compound[from3],
         B = ordered.Compound[to3]) %>% 
  mutate(cycle = str_c(A, "_", B)) %>% 
  # clean up
  ungroup() %>% 
  select(-c(from2, to2, from3, to3, index)) %>% 
  # select(-c(contains("from"), contains("to"))) %>% 
  arrange(cycle) 

# this dataset by this step contains both directions of each cycle
d.paired

# select the minimal flux value of the two opposing reactions in each cycle
d.cycle.CO2 <- d.paired %>% 
  # filter(phenotype == "WT") %>% 
  group_by(phenotype, cycle) %>% 
  filter(nmol.min.animal.mean == min(nmol.min.animal.mean)) %>% 
  # remove `from` and `to`, as the directionality is no longer important
  select(-c(from, to)) %>% 
  arrange(phenotype) %>% 
  # remove futile cycles with minor fluxes 
  filter(nmol.min.animal.mean > 10) %>% 
  arrange(phenotype, desc(nmol.min.animal.mean))

d.cycle.CO2

# remove CO2 entries, keep only the cycles
d.cycle.CO2 %>% 
  filter(B != "CO2") %>% 
  filter(B != "storage")




# - <>- -TAG-NEFA cycle - <>- - <>- - <>- - <>- - <>- - <>- - <>- - <>- 
# rearrange reesterification data, expressed in molecules nmol/min/animal
x <- d.reesterified.tidy %>% 
  filter(destiny != "CO2") %>% select(-c(error.y, fraction))
# mean values
x1 <- x %>% select(-nmol.min.animal.SEM) %>% 
  # convert to molecule molarity
  mutate(nmol.min.animal.mean = nmol.min.animal.mean / 16) %>% 
  spread(destiny, nmol.min.animal.mean)
# SEM
x2 <- x %>% select(-nmol.min.animal.mean) %>% 
  # convert to molecule molarity
  mutate(nmol.min.animal.SEM = nmol.min.animal.SEM / 16) %>% 
  spread(destiny, nmol.min.animal.SEM) %>% 
  rename(ER.SEM = ER, IR.SEM = IR)
# carbon atom flux of esterification
d.ester <- left_join(x1, x2)

# combine with glycerol to sink flux
d.ester2 <- d.flx.CO2.nonOx.sink %>% 
  filter(destiny == "non-Ox sink" & Compounds == "Glycerol") %>% 
  ungroup() %>% select(phenotype, nmol.min.animal.mean, nmol.min.animal.SEM) %>% 
  # convert to molecules basis
  mutate(nmol.min.animal.mean = nmol.min.animal.mean / 3,
         nmol.min.animal.SEM  = nmol.min.animal.SEM / 3) %>% 
  # S.glyL: carbon atoms sink flux of glycerol
  rename(S.gly = nmol.min.animal.mean, S.gly.SEM = nmol.min.animal.SEM) %>% 
  right_join(d.ester)

# combine with glucose to sink flux
d.ester3 <- d.flx.CO2.nonOx.sink %>% 
  filter(destiny == "non-Ox sink" & Compounds == "Glucose") %>% 
  ungroup() %>% select(phenotype, nmol.min.animal.mean, nmol.min.animal.SEM) %>% 
  # convert to molecules basis
  mutate(nmol.min.animal.mean = nmol.min.animal.mean / 6,
         nmol.min.animal.SEM  = nmol.min.animal.SEM / 6) %>% 
  rename(S.glu = nmol.min.animal.mean, S.glu.SEM = nmol.min.animal.SEM) %>% 
  left_join(d.ester2) %>% 
  # reorder columns
  select(phenotype, !contains("SEM"), contains("SEM"))

# counting ATP cost 
d.ester4 <- d.ester3 %>% 
  # needed glycerol molecules for ER 
  mutate(ER.gly.need = ER / 3,
         ER.gly.needed.SEM = ER.SEM / 3) %>% 
  mutate(
    # Sink flux of glycerol to TAG via lipoprotein pathway, cost 11 ATP / TAG
    ATP.ER.lipop = S.gly * 11,
    ATP.ER.lipop.SEM = S.gly.SEM * 11,
    
    # needed glycerol molecules for ER only from lipolysis G3P pathway
    ATP.ER.G3P = (ER.gly.need - S.gly) * 8,
    ATP.ER.G3P.SEM = sqrt(ER.gly.needed.SEM^2 + S.gly.SEM^2) * 8,
    
    ATP.ER.total = ATP.ER.lipop + ATP.ER.G3P,
    ATP.ER.total.SEM = sqrt(ATP.ER.lipop.SEM^2 +  ATP.ER.G3P.SEM^2)
  ) %>% 
  
  # part of ER, and all IR gets its glycerol from glucose via glycolysis G3P pathway
  # mutate(
  #   # glycerol molecule flux to TAG, with glycerol from glycolysis and G3P pathway
  #   gly.from.glu = (ER + IR)/3 - S.gly, 
  #   # glucose molecules that supports to generate glycerol backbone in TAG
  #   glu.to.TAG = gly.from.glu/2,
  #   # such TAG sink flux of glucose relative to total glucose sink flux
  #   frac.glu.TAG =  glu.to.TAG / S.glu) %>% 
  # # ATP cost by the glycolysis G3P pathway
  # mutate(ATP.G3Ppath = gly.from.glu * 8) %>% 

# Total ATP cost for ER reesterification
mutate(ATP.TAG.resynthesis = ATP.ER.total,
       ATP.TAG.resynthesis.SEM = ATP.ER.total.SEM)

d.ester4
# 74~92% of glucose sink flux is directed to TAG glycerol...

d.ATP.TAG <- d.ester4 %>% 
  select(phenotype, ATP.TAG.resynthesis, ATP.TAG.resynthesis.SEM)



# nmol ATP consumed / min / animal
d.ATP.TAG
# ATP fold change
22952/13547

# ER fold change
99.2/43.2

# ER + IR fold change
(28.4+99.2)/(26.4+43.2)

#- --- Protein turnover ---<>----<>----<>----<>----<>----<>-
# average protein oxidized or resynthesized, in mg
MW.aminoAcids = 100

d.ATP.totalAA <- d.flx.CO2.nonOx.sink %>% 
  filter(Compounds == "Valine" & destiny == "non-Ox sink") %>% 
  ungroup() %>% select(phenotype, nmol.min.animal.mean, nmol.min.animal.SEM) %>% 
  # convert to molecules basis, and extrapolate to total amino acids
  mutate(nmol.min.animal.mean = nmol.min.animal.mean / 5 * correct.factor.aminoAcids.molecule,
         nmol.min.animal.SEM  = nmol.min.animal.SEM / 5 * correct.factor.aminoAcids.molecule) %>% 
  # S.glyL: carbon atoms sink flux of glycerol
  rename(totalAA = nmol.min.animal.mean, totalAA.SEM = nmol.min.animal.SEM) %>% 
  # calculate ATP consumed: 4 ATP consumed per peptide bond
  mutate(ATP.totalAA = 4 * totalAA,
         ATP.totalAA.SEM = 4 * totalAA.SEM) %>% 
  select(phenotype, contains("ATP"))

d.ATP.totalAA


# 12.2kcal/mole

#- --- Cori cycle ---<>----<>----<>----<>----<>----<>-
d.ATP.cori <- d.production.consumption %>% ungroup() %>% 
  select(phenotype, from, to, nmol.min.animal.mean, nmol.min.animal.SEM) %>% distinct()  %>% 
  filter(from  == "Lactate" & to == "Glucose") %>% 
  # convert to molecules basis of glucose
  mutate(nmol.min.animal.mean = nmol.min.animal.mean / 6,  
         nmol.min.animal.SEM = nmol.min.animal.SEM / 6 * .6) %>% 
  # 4 ATP cost per cycle of 1 glucose molecule
  mutate(ATP.Cori = nmol.min.animal.mean * 4,
         ATP.Cori.SEM = nmol.min.animal.SEM * 4) %>% 
  select(phenotype, contains("ATP"))


# ATP consumed: nmol ATP / min / animal
d.ATP <- left_join(d.ATP.totalAA, d.ATP.TAG) %>% left_join(d.ATP.cori) %>% 
  select(phenotype, !contains("SEM"), contains("SEM")) 
# simplify column names
colnames(d.ATP)<- colnames(d.ATP)  %>% str_remove("ATP.")

x1 <- d.ATP %>% select(phenotype, !contains("SEM")) %>% 
  pivot_longer(-phenotype, names_to = "cycle", values_to = "wasted") 

x2 <- d.ATP %>% select(phenotype, contains("SEM")) %>% 
  pivot_longer(-phenotype, names_to = "cycle", values_to = "wasted.SEM") %>% 
  mutate(cycle = str_remove(cycle, ".SEM"))

# ATP in nmol/min/animal
d.ATP.tidy <-  left_join(x1, x2, by = c("phenotype", "cycle")) 



# read total energy expenditure (TEE) dataset from 13CO2 measurement experiments
load("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/3_core_13CO2_infusion_physiology_basics.RData")

d.EE <- d.energyExpenditure %>% ungroup() %>% 
  group_by(phenotype) %>% 
  summarise(cal.min.SEM = sd(cal.min) / length(phenotype),
            cal.min = mean(cal.min)) %>% 
  # equivalent ATP of total energy expendigure, nmol/min/animal
  # 12.2 kcal / mole ATP (hydrolysis of one phosphate bond) , or 12.2 cal / mmol
  mutate(TEE.equi = cal.min / (12.2) * 10^6,
         TEE.equi.SEM = cal.min.SEM / 12.2 * 10^6 ) %>% 
  select(phenotype, contains("TEE"))

d.EE # TEE-equivalent ATP value, nmol / min / animal


# combine total ATP consumed with futile cycle ATP
d.ATP.tidy2 <- d.ATP.tidy %>% 
  right_join(d.EE, by = "phenotype") %>% 
  # calculate percentage of waste relative to TEE-equivalent ATP 
  mutate(frac.EE = wasted / TEE.equi,
         frac.EE.SEM = sqrt((wasted.SEM / wasted)^2 + (TEE.equi.SEM / TEE.equi)^2) * frac.EE )

d.ATP.tidy2
# plot


# calculate ATP generated based on direct oxidation flux
d.OX.cateogry <-  d.flx.CO2.nonOx.sink.categorized.summary %>% 
  filter(destiny == "CO2") %>% ungroup() %>% 
  select(phenotype, category, categoryTotal.nmol.min.animal.mean, categoryTotal.nmol.min.animal.SEM)

# Glucose 6 carbon - 32 ATP, or 5 ATP / carbon
# fat: oleic acid gives 118.5 ATP, or ~6.5 ATP / carbon
# protein: assume 6 ATP / carbon

# add # of ATP per carbon oxidation
ATP.perCarbon        <- c( 32/6,   118.5/18,   5,    22/4)
names(ATP.perCarbon) <- c("carbs", "fat", "protein", "KB")

d.ATP.produced <- d.OX.cateogry %>% mutate(ATP.perC = ATP.perCarbon[category]) %>% 
  # ATP total produced: nmol/min/animal
  mutate(ATP.produced = categoryTotal.nmol.min.animal.mean * ATP.perC,
         ATP.produced.SEM = categoryTotal.nmol.min.animal.SEM * ATP.perC)

d.ATP.produced.summary <- d.ATP.produced %>% 
  group_by(phenotype) %>% 
  summarise(prod.byCO2 = sum(ATP.produced),
            prod.byCO2.SEM = sum(ATP.produced.SEM^2) %>% sqrt())
d.ATP.produced.summary


# calculate O2 consumption and P-O ratio based ATP production
PO = 2
d.ATP.produced.summary.byOxy <- d.energyExpenditure %>% ungroup() %>% 
  group_by(phenotype) %>% 
  summarise(O2.L.min.SEM = sd(O2.L.min) / length(phenotype),
            O2.L.min = mean(O2.L.min)) %>% 
  # calculate ATP production based on PO ratio
  mutate(prod.byOxy = O2.L.min / 22.4 * 2 * 10^9 * PO,
         prod.byOxy.SEM = O2.L.min.SEM / 22.4 * 2 * 10^9 * PO) %>% 
  select(phenotype, contains("prod"))


# combine with the dataset of ATP wasted, TEE-equivalent ATP, and total ATP production estimated by CO2 and Oxygen
d.ATP.tidy3 <- d.ATP.tidy2 %>% left_join(d.ATP.produced.summary, by = "phenotype") %>% 
  # calculate percentage of waste relative to actually produced ATP calculated based on CO2 data
  mutate(frac.prod.CO2 = wasted / prod.byCO2,
         frac.prod.CO2.SEM = sqrt((wasted.SEM / wasted)^2 + (prod.byCO2.SEM / prod.byCO2)^2) * frac.prod.CO2 ) %>% 
  
  # combine ATP produced estimated by oxygen usage
  left_join(d.ATP.produced.summary.byOxy, by = "phenotype") %>% 
  
  # how different is the ATP production estimate by O2 and CO2 measurement
  mutate(prod.CO2.vs.Oxy = prod.byCO2 / prod.byOxy) %>% 
  # ATP production efficiency
  mutate(efficiency = prod.byCO2 / TEE.equi) 

d.ATP.tidy3


# compare the total ATP produced calculated by different means

d.ATP.tidy3 %>% 
  select(phenotype, TEE.equi, contains("prod."), -contains("frac"), -contains("SEM"))




# check total ATP produced by 2 methods: CO2 based vs. O2 based methods
d.ATP.tidy3 %>% select(phenotype, contains("prod"), -contains("frac")) %>% 
  distinct()


d.ATP.tidy4 <- d.ATP.tidy3 %>% 
  mutate(cycle = factor(cycle, levels = c("TAG.resynthesis", "totalAA", "Cori") %>% rev())) %>% 
  arrange(desc(cycle)) %>% 
  group_by(phenotype) %>% 
  mutate( err.y.wasted = cumsum(wasted),
          err.y.EE = cumsum(frac.EE),
          err.y.prod.CO2 = cumsum(frac.prod.CO2)) 

d.ATP.tidy4

# visualization
d.ATP.tidy4.WT <- d.ATP.tidy4 %>% filter(phenotype == "WT")

frac.max.prod.CO2 <- d.ATP.tidy4.WT$frac.prod.CO2 %>% sum()

d.ATP.tidy4.WT$frac.EE %>% sum()

d.ATP.tidy4.WT %>% 
  ggplot(aes(x = 1.7, y = wasted, fill = cycle)) +
  
  # anchoring point
  # geom_point(aes(x=  .5, y = 20000), inherit.aes = F, color = "grey", alpha = 0) +
  # geom_point(aes(x=  3.5, y = 20000), inherit.aes = F, color = "grey", alpha = 0) +
  
  geom_col(color = "black", alpha = .5, width = .7) +
  geom_errorbar(aes(ymin = err.y.wasted - wasted.SEM,
                    ymax = err.y.wasted ),
                width = .2) +
  # secondary axis on the primary axis
  annotate(geom = "text", 
           x = (s <- .12),
           label = ((u <- seq(.02, frac.max.prod.CO2 + .02, .02)) * 100) %>% paste("%"),
           y = (j <- u * unique(d.ATP.tidy4.WT$prod.byCO2)), # number of ATP wasted
           hjust = 0, size = 4.5) +
  # secondary ticks on y axis on the left
  annotate(geom = "segment", 
           x = 0, xend = .06, 
           y = j, yend = j) +
  scale_x_continuous(expand = expansion(mult = c(0, .2))) +
  
  scale_y_continuous(
    # limits = c(0, sum(d.ATP.tidy4.WT$wasted)),
    expand = expansion(mult = c(0, .1)),
    labels = function(x){x/1000},
    breaks = seq(0, 35 * 1000, 5 * 1000),
    name = "ATP wasted  mmole / min",
    sec.axis = sec_axis(trans = ~. / unique(d.ATP.tidy4.WT$TEE.equi) * 100,
                        # breaks = seq(.05, 5, 1),
                        labels = function(x) {paste(x, "%")} )
  ) +
  # scale_x_discrete(expand = expansion(add = .8)) +
  scale_fill_brewer(palette = "Accent", 
                    label = c("Cori" = "Cori cycle", 
                              "TAG.resynthesis" = "TAG synthesis", 
                              "totalAA" = "protein synthesis")) +
  theme.myClassic +
  theme(axis.title.y.left = element_text(margin = margin(r = 10)),
        axis.title.y.right = element_text(margin = margin(l = 15)),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y.left = element_text(margin = margin(r = 5, unit = "pt")),
        axis.text.y.right = element_text(margin = margin(r = 5, l = 5, unit = "pt")),
        axis.line.y.right = element_line(linetype = "dashed"),
        legend.position = "bottom") 




ggsave(filename = "ATP cost_WT.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 3.5, width = 3.4)


# plot all 3 phenotypes

# fraction of total whole body ATP consumption (energy expenditure)

func.polish <- function(plt) {
  plt + theme.myClassic +
    scale_x_discrete(expand = expansion(add = .8),
                     labels = function(x){str_replace(x, "WT", "control")}) +
    scale_fill_brewer(palette = "Accent", 
                      label = c("Cori" = "Cori cycle", 
                                "TAG.resynthesis" = "TAG resynthesis", 
                                "totalAA" = "protein resynthesis")) +
    theme(axis.title.y = element_text(margin = margin(r = 10)),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) 
}

plt.ATP.wasted <-  (
  d.ATP.tidy4 %>% 
    ggplot(aes(x = phenotype, y = wasted, fill = cycle)) +
    geom_col(color = "black", alpha = .5) +
    geom_errorbar(aes(ymin = err.y.wasted - wasted.SEM,
                      ymax = err.y.wasted ),
                  width = .2) +
    scale_y_continuous(expand = expansion(mult = c(0, .1)),
                       labels = function(x){x / 1000},
                       name = "ATP wasted  µmol / min")) %>% 
  func.polish() 

# % of total ATP produced
plt.frac.prod.CO2 <- (
  d.ATP.tidy4 %>% 
    ggplot(aes(x = phenotype, y = frac.prod.CO2, fill = cycle)) +
    geom_col(color = "black", alpha = .5) +
    geom_errorbar(aes(ymin = err.y.prod.CO2 - frac.prod.CO2.SEM,
                      ymax = err.y.prod.CO2 ),
                  width = .2) +
    scale_y_continuous(expand = expansion(mult = c(0, .1)),
                       labels = function(x){x * 100},
                       breaks = seq(0, .16, .02),
                       name = "% of total ATP production"))  %>% 
  func.polish()
plt.frac.prod.CO2

# % of total EE
plt.frac.TEE <- (
  d.ATP.tidy4 %>% 
    ggplot(aes(x = phenotype, y = frac.EE, fill = cycle)) +
    geom_col(color = "black", alpha = .5) +
    geom_errorbar(aes(ymin = err.y.EE - frac.EE.SEM,
                      ymax = err.y.EE ),
                  width = .2) +
    scale_y_continuous(expand = expansion(mult = c(0, .1)),
                       labels = function(x){x * 100},
                       breaks = seq(0, .08, .005),
                       name = "% of total energy expenditure"))  %>% 
  func.polish()

plt.frac.TEE

# plot all
plot_grid(plt.ATP.wasted + theme(legend.position = "none"), 
          ggplot() + theme_void(),
          
          plt.frac.prod.CO2 + theme(legend.position = "none"),
          ggplot() + theme_void(),
          
          plt.frac.TEE,
          
          nrow = 1, rel_widths = c(2.4, .1, 
                                   2.4, .1, 
                                   5))

ggsave(filename = "ATP cost_3 pheno.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 3.5, width = 10)


#---<>-#---<>-#---<>-#---<>-#---<>-#---<>-#---<>-
# glycolysis and TCA flux of circulating glucose

# total glycolysis rate of circulating glucose
# glucose to direct flux to CO2
d.glc.to.CO2 <- d.production.consumption %>% 
  filter(from == "Glucose" & to == "CO2") %>% ungroup() %>% 
  select(phenotype,  nmol.min.animal.mean, nmol.min.animal.SEM) %>% 
  rename(glc.to.CO2 = nmol.min.animal.mean,
         glc.to.CO2.SEM = nmol.min.animal.SEM)

# lactate direct flux to CO2
d.lac.to.CO2 <- d.production.consumption %>% 
  filter(from == "Lactate" & to == "CO2") %>% ungroup() %>% 
  select(phenotype, from, to, nmol.min.animal.mean, nmol.min.animal.SEM)

# glucose contribution fraction to lactate
d.glc.to.lac <- d.production.consumption %>% ungroup() %>% 
  select(phenotype, from, to, contains("animal")) %>% distinct() %>% 
  group_by(phenotype, to) %>% 
  mutate(contri.frac = nmol.min.animal.mean / sum(nmol.min.animal.mean)) %>% 
  select(phenotype, from, to, contri.frac) %>% 
  filter(from == "Glucose" & to == "Lactate") %>% distinct() %>% 
  rename(frac.glu.to.lac = contri.frac) %>% ungroup() %>% 
  select(phenotype, frac.glu.to.lac)

# carbon flux of glucose to CO2 via lactate
d.glc.CO2.via.lac <- d.lac.to.CO2 %>% left_join(d.glc.to.lac, by = "phenotype") %>% 
  mutate(glc.CO2.via.lac = nmol.min.animal.mean * frac.glu.to.lac,
         glc.CO2.via.lac.SEM = nmol.min.animal.SEM * frac.glu.to.lac) %>% 
  select(phenotype, contains("glc.CO2.via.lac"))

# calculate total TCA flux of circulating glucose (direct Ox or via lactate)
d.circ.glc.TCA <- d.glc.to.CO2 %>% left_join(d.glc.CO2.via.lac, by = "phenotype") %>% 
  group_by(phenotype) %>% 
  mutate(glc.TCA = glc.to.CO2 + glc.CO2.via.lac,
         glc.TCA.SEM = sqrt(glc.to.CO2.SEM^2 + glc.CO2.via.lac.SEM^2),
         .keep = "none")
d.circ.glc.TCA

# an alternative way to calculate TCA flux, via overall oxidation flux
d.circ.glc.TCA <- d.ox.sink.overal.2 %>% 
  filter(Compound == "Glucose" & destiny == "ox.overal") %>% 
  filter(phenotype != "db/db") %>% 
  ungroup() %>% 
  select(phenotype, contains("anim") ) %>% 
  rename(glc.TCA = nmol.min.animal,
         glc.TCA.SEM = nmol.min.animal.SEM)



# glc to lactate flux
d.glc.to.lac <- d.production.consumption %>% 
  filter(from == "Glucose" & to == "Lactate") %>% ungroup() %>% 
  select(phenotype, from, to, nmol.min.animal.mean, nmol.min.animal.SEM) %>% distinct() %>% 
  rename(glc.to.lac = nmol.min.animal.mean,
         glc.to.lac.SEM = nmol.min.animal.SEM) %>% 
  select(phenotype, contains("glc.to.lac"))

# glycolysis flux of circulating glucose
d.glycolysis <- d.glc.to.CO2 %>% left_join(d.glc.to.lac) %>% 
  group_by(phenotype) %>% 
  mutate(glycolysis = glc.to.CO2 + glc.to.lac,
         glycolysis.SEM = sqrt(glc.to.CO2.SEM^2 + glc.to.lac.SEM^2),
         .keep = "none")

# combine glycolysis and TCA flux
x1 <- d.circ.glc.TCA %>% rename(C.mean = glc.TCA, C.SEM = glc.TCA.SEM) %>% 
  mutate(path = "TCA",
         ATP.mean = C.mean / 6 * 28, 
         ATP.SEM =  C.SEM / 6 * 28)

x2 <- d.glycolysis %>% rename(C.mean = glycolysis, C.SEM = glycolysis.SEM) %>%
  mutate(path = "Glycolysis",
         ATP.mean = C.mean / 6 * 2, 
         ATP.SEM =  C.SEM / 6 * 2)

d.glc.glycolysis.TCA <- rbind(x1, x2) %>% 
  mutate(path = factor(path, levels = c("TCA", "Glycolysis"))) %>% 
  group_by(phenotype) %>% 
  arrange(desc(path)) %>% 
  mutate(y.err.C = cumsum(C.mean),
         y.err.ATP = cumsum(ATP.mean)) %>% 
  arrange(phenotype)
d.glc.glycolysis.TCA

8640 / 6710

x1 <- rnorm(3, mean = 8640, sd = 939)
x2 <- rnorm(3, mean = 6710, sd = 365)
t.test(x1, x2)


# visualize
func.plt.glycolysis.TCA <-  function(mydata){
  
  # mydata <- d.glc.glycolysis.TCA
  plt1 <- mydata %>% 
    ggplot(aes(x = phenotype, y = C.mean, fill = path)) +
    geom_col(color = "black", position = "dodge") +
    geom_errorbar(aes(ymax = C.mean + C.SEM,
                      ymin = C.mean - C.SEM),
                  width = .2,
                  position = position_dodge(.9)) +
    # text
    # geom_text(aes(label = round(C.mean/1000, 1) %>% paste("\nµmol/min")),
    #           position = position_stack(vjust = .5)) +
    
    scale_x_discrete(expand = expansion(add = .8),
                     labels = function(x){str_replace(x, "WT", "control")}) +
    scale_y_continuous(expand = expansion(mult = c(0, .1)),
                       # breaks = seq(0, 15, 3)* 1000, 
                       n.breaks = 7,
                       labels = function(x){x/1000}, 
                       name = "C µmol / min\n") + 
    theme.myClassic +
    theme(legend.position = "bottom",
          text = element_text(size = 21)) +
    scale_fill_manual(values = c("grey60", "snow")) 
  
  plt2 <- mydata %>% 
    ggplot(aes(x = phenotype, y = ATP.mean, fill = path)) +
    geom_col(color = "black", position = position_dodge()) +
    geom_errorbar(aes(ymax = ATP.mean + ATP.SEM,
                      ymin = ATP.mean - ATP.SEM),
                  width = .2, 
                  position = position_dodge(.9)) +
    scale_x_discrete(expand = expansion(add = .8),
                     labels = function(x){str_replace(x, "WT", "control")}) +
    scale_y_continuous(expand = expansion(mult = c(0, .1)),
                       n.breaks = 7,
                       labels = function(x){x/1000}, 
                       name = "ATP µmol / min\n") +
    theme.myClassic +
    theme(legend.position = "bottom",
          text = element_text(size = 21)) +
    scale_fill_manual(values = c("grey60", "snow"))
  
  plot_grid(plt1, plt2, nrow = 1)
}

d.glc.glycolysis.TCA %>% 
  filter(phenotype == "WT") %>% func.plt.glycolysis.TCA()

ggsave(filename = "glycolysis TCA_WT.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 4.5, width = 5)

d.glc.glycolysis.TCA %>% func.plt.glycolysis.TCA()

ggsave(filename = "glycolysis TCA_3pheno.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 5, width = 7)



#-,>-#-,>-#-,>-#-,>-#-,>-#-,>-#-,>-#-,>-#-,>-#-,>-#-,>-
# Calculate ATP Gibbs free energy 

conc.ATP <- seq(from = 1, to = 10, 1) 
conc.Pi <-  seq(from = 1, to = 10, 1) 
conc.ADP <- seq(from = .1, to = 1, .1) 

d.ATP.DG <- tibble(conc.ATP, conc.ADP, conc.Pi)

DG.standard <- -31.5
R = 8.314 / 1000 # KJ/mol/K
K = 310 # Kelvin

d.DG <- d.ATP.DG %>% 
  # create all possible combinations
  expand(conc.ATP, conc.ADP, conc.Pi) %>% 
  mutate(DG.physio = DG.standard + R * K * log( conc.ADP/1000  *  conc.Pi/1000  / (conc.ATP /1000) ),
         DG.physio = -DG.physio)

d.DG


d.DG %>% 
  ggplot(aes(x = conc.ATP, y = DG.physio, color = conc.ADP)) +
  geom_line(aes(group = conc.ADP)) + 
  scale_color_distiller(palette = "Spectral", breaks = seq(0.1, 1, .2)) +
  facet_wrap(~conc.Pi, nrow = 1) +
  scale_x_continuous(breaks = seq(1, 11, 4)) +
  
  guides(color = guide_colorbar(barwidth = unit(300, "pt"),
                                barheight = unit(5, "pt"))) +
  
  geom_point(data = tibble(
    conc.ATP = 5, DG.physio = 51, conc.ADP = 5, conc.Pi = 5),
    color = "black", shape = "diamond", size = 3) +
  geom_hline(yintercept = 51, linetype = "dashed", linewidth = .5) +
  scale_y_continuous(
    name = "Energy (KJ/mol)\n",
    sec.axis = sec_axis(trans = ~(. - 51)/(51)*100,
                        labels = function(x){paste(x, "%")},
                        name = "relative error"
    )
  ) +
  labs(x = "ATP (mM)", 
       color = "ADP (mM)    ") +
  theme.myClassic +
  theme(legend.position = "bottom",
        panel.spacing.x = unit(10, "pt"),
        legend.title = element_text(color = "black", face = "bold", size = 13),
        axis.text.x = element_text(size = 12)) 


ggsave(filename = "ATP Gibbs Free Energy.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 4, width = 10)
