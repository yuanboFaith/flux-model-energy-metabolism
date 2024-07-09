rm(list = ls())

library(writexl)
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
library(broom)
library(tidyverse)


load(file = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/5_core_labeling_analysis.RData")



# combine with body weight
# add body weight
d.BW = d.normalized.tidy %>% 
  select(mouse_ID, phenotype, BW, infused.tracer) %>% 
  distinct() %>% 
  group_by(phenotype, infused.tracer) %>% 
  summarise(BW.mean = mean(BW),
            BW.sd = sd(BW)) %>% 
  rename(Compounds = infused.tracer)
d.BW




# Integrate with 13CO2 data: consumption fate of circulating nutrients on a whole body level
d.CO2.fox = read_excel(
  "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/nutrients 13CO2 recovery and fox.xlsx", 
  sheet = "fox_pooled")
d.CO2.fox = d.CO2.fox %>% select(-1) %>% rename(infused.tracer = tracer)

colnames(d.CO2.fox) <- 
  colnames(d.CO2.fox) %>% str_replace("fox", replacement = "fox") 

# for lactate, force fox.mean to be a max of 1
d.CO2.fox <- d.CO2.fox %>% 
  mutate(fox.mean = ifelse(fox.mean > 1, 1, fox.mean))


# Calculate overall oxidation and sink fluxes

d.ox.sink.overal.0 <-  d.Fcirc_standard.atom.summary %>% 
  # select the classic Fcirc-atom quantification
  select(Compound, phenotype, n.rep, contains("atom"), -contains("position")) %>% 
  
  # select whole animal data
  select(-contains("g.BW")) %>% 
  
  mutate(Fcirc_animal.atom.SEM = Fcirc_animal.atom.sd / sqrt(n.rep)) %>% 
  select(-c(Fcirc_animal.atom.sd, n.rep)) %>% 
  left_join(d.CO2.fox %>% rename(Compound = infused.tracer) %>% select(-c(n.rep, method, fox.SD)),
            by = c("phenotype", "Compound")) 

d.ox.sink.overal.1 <- d.ox.sink.overal.0 %>% 
  mutate(
    ox.overal.nmol.min.animal = Fcirc_animal.atom.mean * fox.mean,
    sink.overal.nmol.min.animal = Fcirc_animal.atom.mean * (1-fox.mean),
    
    ox.overal.nmol.min.animal.SEM = 
      sqrt((Fcirc_animal.atom.SEM/Fcirc_animal.atom.mean)^2 + (fox.SEM/fox.mean)^2) * ox.overal.nmol.min.animal,
    
    sink.overal.nmol.min.animal.SEM = ox.overal.nmol.min.animal.SEM) %>% 
  select(-c(contains("Fcirc"), contains("fox"))) %>% 
  mutate(Compound = factor(Compound, levels = ordered.Compound),
         phenotype = factor(phenotype, levels = ordered.phenotype))

x1 <- d.ox.sink.overal.1 %>% 
  select(-contains("SEM")) %>% 
  pivot_longer(-c(Compound, phenotype), names_to = "destiny", values_to = "nmol.min.animal") %>% 
  mutate(destiny = str_remove(destiny, ".nmol.min.animal"))

x2 <- d.ox.sink.overal.1 %>% 
  select(Compound, phenotype, contains("SEM")) %>% 
  pivot_longer(-c(Compound, phenotype), names_to = "destiny", values_to = "nmol.min.animal.SEM") %>% 
  mutate(destiny = str_remove(destiny, ".nmol.min.animal.SEM"))

d.ox.sink.overal.2 <- x1 %>% 
  left_join(x2, by = c("Compound", "phenotype", "destiny")) %>% 
  mutate(destiny = factor(destiny, levels = c("ox.overal", "sink.overal"))) %>% 
  arrange(desc(destiny)) %>% 
  group_by(Compound, phenotype) %>% 
  mutate(yAxis.error = cumsum(nmol.min.animal))


# plot 1 : WT dodged position

u <- d.ox.sink.overal.2 %>% ungroup() %>% 
  filter(phenotype == "WT") %>%
  filter(destiny == "ox.overal") %>% 
  mutate(Compound = fct_reorder(Compound, nmol.min.animal, sum, .desc = T)) %>% 
  arrange(Compound)

# WT flux per animal
plt.overal.ox.sink.WT.dodge <-   
  d.ox.sink.overal.2 %>% ungroup() %>% 
  filter(phenotype == "WT") %>% 
  mutate(Compound = factor(Compound, levels = u$Compound)) %>% 
  ggplot(aes(x = Compound, y = nmol.min.animal, fill = destiny)) +
  geom_col(color = "black", width = .7, alpha = .7, position = position_dodge(.7)) +
  geom_errorbar(aes(ymax = nmol.min.animal + nmol.min.animal.SEM,
                    ymin = nmol.min.animal - nmol.min.animal.SEM),
                width = .2,
                position = position_dodge(.7)) +
  scale_y_continuous(breaks = seq(0, 16000, 4000),
                     position = "left",
                     labels = function(x){x/1000},
                     name = "overal fluxes (µmol / min )",
                     expand = expansion(mult = c(0, .1))) +
  scale_fill_manual(values = c("orange", "black"),
                    labels = c("OX", "NOS")) + # NOX for non-oxidative sink
  theme.myClassic +
  theme(
    legend.text = element_text(size = 15),
    legend.key.size = unit(20, "pt"),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(.84, .8),
    axis.title.y = element_text(margin = margin(r = 10, "pt")),
    plot.margin = margin(l = 20)) +
  labs(x = NULL) +
  coord_cartesian(ylim = c(0, NA))

plt.overal.ox.sink.WT.dodge

ggsave(
  filename = "ox.sink.overal.WT_dodge.pdf",
  path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures", 
  width = 7/1.1, height = 3.5/1.1)



# WT flux per g BW

d.BW %>% group_by(phenotype) %>% summarise(BW = mean(BW.mean))

plt.overal.ox.sink.WT.dodge +
  scale_y_continuous(
    breaks = seq(0, 16000, 4000),
    position = "right",
    labels = function(x){x/1000},
    name = "overal fluxes (µmol / min )",
    expand = expansion(mult = c(0, .1)),
    sec.axis = sec_axis(trans =  ~ . / 29, name = "nmol / min / g")) 

ggsave(
  filename = "ox.sink.overal.WT_dodge.gBW.pdf",
  path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures", 
  width = 7, height = 3.5/1.1)


# plot 2: WT, bars stacked
plt.overal.ox.sink.WT <- d.ox.sink.overal.2 %>% ungroup() %>% 
  filter(phenotype == "WT") %>% 
  mutate(Compound = fct_reorder(Compound, nmol.min.animal, sum, .desc = T)) %>% 
  ggplot(aes(x = Compound, y = nmol.min.animal, fill = destiny)) +
  geom_col(color = "black", width = .7, alpha = .7) +
  geom_errorbar(aes(ymax = yAxis.error,
                    ymin = yAxis.error - nmol.min.animal.SEM),
                width = .2) +
  scale_y_continuous(n.breaks = 8,
                     position = "left",
                     labels = function(x){x/1000},
                     name = "overal fluxes (µmol / min )",
                     expand = expansion(mult = c(0, .1))) +
  scale_fill_manual(values = c("orange", "black"),
                    labels = c("ox", "non-ox sink")) +
  theme.myClassic +
  theme(
    legend.text = element_text(size = 15),
    legend.key.size = unit(20, "pt"),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(.75, .8),
    axis.title.y = element_text(margin = margin(r = 10, "pt")),
    plot.margin = margin(l = 20)) +
  labs(x = NULL) +
  coord_cartesian(ylim = c(0, NA))

plt.overal.ox.sink.WT

ggsave(
  filename = "ox.sink.overal.WT.pdf",
  path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures", 
  width = 6, height = 3.5)



# plot 3: all three phenotypes
plt.overal.ox.sink.3pheno <-  d.ox.sink.overal.2 %>% 
  filter(phenotype != "db/db") %>% 
  ggplot(aes(x = phenotype, y = nmol.min.animal, fill = destiny)) +
  geom_col(color = "black", alpha = .7) +
  geom_errorbar(aes(ymax = yAxis.error,
                    ymin = yAxis.error - nmol.min.animal.SEM),
                width = .5) +
  facet_wrap(~Compound, nrow = 2, scales = "free") +
  scale_y_continuous(labels = function(x){x/1000},
                     name = "overal fluxes (µmol / min / animal)",
                     expand = expansion(mult = c(0, .1))) +
  scale_x_discrete(expand = expansion(add = 1),
                   labels = function(x){str_replace(x, "WT", "control")}) +
  theme.myClassic +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(r = 10, "pt"))) +
  scale_fill_manual(values = c("orange", "black"),
                    labels = c("ox.overal" = "oxidation",
                               "sink.overal" = "storing")) +
  coord_cartesian(ylim = c(0, NA))

plt.overal.ox.sink.3pheno

ggsave(
  filename = "ox.sink.overal.3pheno.pdf",
  path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures", 
  width = 10, height = 6)


plt.overal.ox.sink.3pheno +
  geom_text(aes(label = round(nmol.min.animal)),
            position = position_stack(vjust = .4),
            color = "red") +
  # add percentage
  geom_text(data = d.ox.sink.overal.2 %>% mutate(frac = nmol.min.animal /sum(nmol.min.animal) ) %>% 
              filter(phenotype != "db/db"),
            aes(label = round(frac*100, 1) %>% paste("%")),
            position = position_stack(vjust = .6), color = "cyan",
            size = 3)

ggsave(
  filename = "ox.sink.overal.3pheno_numbered.pdf",
  path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures", 
  width = 15, height = 8)


# <>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-# <>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-<>-

# Calculate direct fluxes of oxidation / sink 

func.flx.CO2.nonOx.sink = function(whichDestiny = "CO2"){
  
  # absolute fraction of oxidation 
  d.absolute.fox = tibble() 
  d.fittedError = tibble()
  
  for (myPhenotype in c("WT", "ob/ob", "HFD")){ 
    
    # myPhenotype = "WT"
    
    d.enrich.atom.summary_Pheno.i = d.enrich.atom.summary %>% 
      filter(phenotype == myPhenotype)
    
    # labeling of all-set metabolites under tracer of all-set metabolites
    d_L.mean = d.enrich.atom.summary_Pheno.i %>% ungroup() %>% 
      select(infused.tracer, Compound, enrich.original.mean) %>% 
      spread(Compound, enrich.original.mean) 
    # convert to matrix with tracer being row names
    mat_L.mean = d_L.mean %>% select(-1) %>% as.matrix()
    rownames(mat_L.mean) = d_L.mean$infused.tracer
    
    # matrix dimension or vector length
    mat.dim <- nrow(mat_L.mean)
    
    # SEM of labeling
    d_L.SEM = d.enrich.atom.summary_Pheno.i %>% ungroup() %>% 
      select(infused.tracer, Compound, enrich.original.sem) %>% 
      spread(Compound, enrich.original.sem) 
    mat_L.SEM = d_L.SEM %>% select(-1) %>% as.matrix()
    
    
    # infusion rate of all-set tracers
    d.infusion_pheno.i = d.normalized.tidy %>% filter(phenotype == myPhenotype) %>% 
      select(infused.tracer, infusion_nmol13C.atoms_perMin.gBW) %>% distinct()
    
    # arrange compounds in ordered factor
    d.infusion_pheno.i = d.infusion_pheno.i %>% 
      # tracer as arranged ordered factor
      mutate(infused.tracer = factor(infused.tracer, levels = ordered.Compound, ordered = F)) %>%
      arrange(infused.tracer) %>% 
      filter(infused.tracer %in% (d.enrich.atom.summary_Pheno.i$Compound %>% unique() )) %>% 
      # add 13CO2 fox column
      left_join(d.CO2.fox %>% filter(phenotype == myPhenotype)) 
    
    # Infusion rate R 
    R_pheno.i <- d.infusion_pheno.i$infusion_nmol13C.atoms_perMin.gBW
    mat.R.pheno.i.reciprocal <- 1/matrix(rep(R_pheno.i, each = mat.dim ), nrow = mat.dim, byrow = T)
    
    
    # the right-side vector
    if (whichDestiny == "CO2"){
      # for oxidative CO2 destiny, use fox
      v_R.fox.mean = d.infusion_pheno.i$fox.mean 
    } else if ( whichDestiny == "non-Ox sink") { 
      # use (1 - fox) for non-oxidative sink; using the same object name as CO2 case
      # then what is solved is the absolute fraction of non-Ox sink
      v_R.fox.mean = (1 - d.infusion_pheno.i$fox.mean)
    }
    
    names(v_R.fox.mean) = d.infusion_pheno.i$infused.tracer
    v_R.fox.SEM = d.infusion_pheno.i$fox.SEM
    
    
    # Check L matrix compounds order match with R.rec vector
    if(! (!rownames(mat_L.mean) == names(v_R.fox.mean)) %>% sum() == 0){stop("Compounds names not matching.")}
    
    
    # constraint matrix
    G = diag(1, nrow = mat.dim) # absolute fraction of oxidation >= 0
    H = rep(0, mat.dim) # flux to CO2 >= 0
    
    
    # Simulate with Monte Carlo for the selected phenotype
    d.absolute.fox_pheno.i = c()
    d.fittedError_pheno.i = c()
    
    for (i in 1: (n.run <- 100) ){
      
      # standard normal matrix and vector
      mat.standardGaussain <- rnorm(n = mat.dim^2, mean = 0, sd = 1) %>% matrix(mat.dim)
      v.standardGaussain <- rnorm(n = mat.dim, mean = 0, sd = 1)
      
      # Simulating terms of the matrix 
      mat_L.MonteCarlo = mat_L.mean + mat_L.SEM * mat.standardGaussain
      M_MonteCarlo = mat.R.pheno.i.reciprocal * (mat_L.MonteCarlo / (1 - mat_L.MonteCarlo))
      
      # Simulating terms on the right side of the equation sign
      v_R.fox.MonteCarlo = v_R.fox.mean + v_R.fox.SEM * v.standardGaussain
      
      # Solve with non-negative constraints
      list.absolute.fox = lsei(A = M_MonteCarlo,
                               B = v_R.fox.MonteCarlo, # this solves Ax = B
                               # non-negative direct flux to CO2
                               G = G, # G is the identity matrix, dimension = B-vector length
                               H = H, # vector of zero, length = B-vector length
                               verbose = FALSE)
      
      # keep recording the absolute fox for each iteration
      d.absolute.fox_pheno.i = rbind(d.absolute.fox_pheno.i, list.absolute.fox$X)
      
      # calculate the fitting values and fitted error
      fitted.value <- (M_MonteCarlo %*% (list.absolute.fox$X)) %>% as.data.frame() # absolute fitted error
      # residual: actual - fitted
      fitted.error.rel.frac <- (fitted.value$V1 - v_R.fox.MonteCarlo )/v_R.fox.MonteCarlo
      fitted.error.rel.frac
      d.fittedError_pheno.i = rbind(d.fittedError_pheno.i, fitted.error.rel.frac)
    }
    
    # solved fluxes
    d.absolute.fox_pheno.i = d.absolute.fox_pheno.i %>% 
      as_tibble() %>% mutate(phenotype = myPhenotype) 
    # combine all phenotypes solved vector
    d.absolute.fox = rbind(d.absolute.fox, d.absolute.fox_pheno.i)
    
    # fitted error
    d.fittedError_pheno.i <- d.fittedError_pheno.i %>% as_tibble() %>% 
      mutate(phenotype = myPhenotype, run = 1:n.run)
    # combine all phenotype errors
    d.fittedError <- rbind(d.fittedError, d.fittedError_pheno.i)
    
  }
  return(list(fluxes = d.absolute.fox, fitted.rel.error = d.fittedError))
}


# First check the fitted relative error fraction
d.fitted.relErr.CO2 <- func.flx.CO2.nonOx.sink(whichDestiny = "CO2")$fitted.rel.error %>% mutate(destiny = "CO2") 
d.fitted.relErr.sink <- func.flx.CO2.nonOx.sink(whichDestiny = "non-Ox sink")$fitted.rel.error %>% mutate(destiny = "non-Ox sink")
d.fitted.relErr.tidy <- rbind(d.fitted.relErr.CO2, d.fitted.relErr.sink) %>% 
  pivot_longer(-c(phenotype, destiny, run), names_to = "compound", values_to = "fitted.relErr") %>% 
  mutate(compound = factor(compound, levels = ordered.Compound),
         phenotype = factor(phenotype, levels = ordered.phenotype))

# plot relative fitted error fraction
d.fitted.relErr.tidy %>% 
  ggplot(aes(x = compound, y = fitted.relErr, color = phenotype)) +
  # geom_quasirandom(alpha = .8, size = 1) +
  # geom_hline(yintercept = 0) +
  # geom_hline(yintercept = -.1, linetype = "dashed", color = "snow4") +
  # geom_hline(yintercept = .1, linetype = "dashed", color = "snow4") +
  geom_boxplot(position = position_dodge(.8)) +
  facet_wrap(~destiny, nrow = 2, scales = "free") +
  # scale_fill_manual(values = color.phenotype) +
  scale_color_manual(values = color.phenotype,
                     labels = function(x){str_replace(x, "WT", "control")}) +
  # (fitted value - actual) / actual value
  scale_y_continuous(breaks = seq(-.8, .8, .2), name = "relative fitted error fraction\n") +
  coord_cartesian(ylim = c(-.5, .5)) +
  theme.myClassic +
  theme(axis.text.x = element_text(hjust = 1, angle = 60))

ggsave(filename = "direct sink ox fitted error rel frac.pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 6, width = 7)

d.fitted.relErr.tidy %>% 
  filter(destiny != "CO2") %>% 
  filter(compound == "Lactate")


# absolute flux to oxidation 
d.OX = func.flx.CO2.nonOx.sink(whichDestiny = "CO2")$fluxes %>% mutate(destiny = "CO2")
# output flux to non-oxidative sink
d.NOS = func.flx.CO2.nonOx.sink(whichDestiny = "non-Ox sink")$fluxes %>% mutate(destiny = "non-Ox sink")
# combine and calculate mean and SD
d.flx.CO2.nonOx.sink <- rbind(d.OX, d.NOS) %>% 
  pivot_longer(-c(phenotype, destiny), names_to = "Compounds", values_to = "flx") %>% 
  group_by(phenotype, destiny, Compounds) %>% 
  summarise(nmol.min.g.mean = mean(flx),
            nmol.min.g.SEM = sd(flx)) %>% 
  mutate(Compounds = factor(Compounds, levels = ordered.Compound),
         phenotype = factor(phenotype, levels = ordered.phenotype)) 




# note that the direct OX and NOS fluxes are calculated based on all infusion data, 
# but the body weight is stratified based on each tracer; this may lead to small differences vs. Tony's Fcirc values
d.flx.CO2.nonOx.sink <- d.flx.CO2.nonOx.sink %>% 
  left_join(d.BW, by = c("phenotype", "Compounds")) %>% 
  mutate(nmol.min.animal.mean = nmol.min.g.mean * BW.mean,
         nmol.min.animal.SEM = nmol.min.g.SEM * BW.mean)

d.flx.CO2.nonOx.sink

# Visualize
# 1) plot WT direct OX and NOS flux per animal, dodged bars
d.flx.CO2.nonOx.sink %>% 
  filter(phenotype == "WT") %>%
  mutate(Compounds = fct_reorder(Compounds, nmol.min.animal.mean, .desc = T)) %>% 
  ggplot(aes(x = Compounds, y = nmol.min.animal.mean, fill = destiny)) +
  geom_col(position = position_dodge(.8),
           alpha = .7, color = "black", linewidth = .5,
           width = .8) +
  geom_errorbar(aes(ymin = nmol.min.animal.mean - nmol.min.animal.SEM,
                    ymax = nmol.min.animal.mean + nmol.min.animal.SEM),
                width = .5, 
                position = position_dodge(.8)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)),
                     labels = function(x)x/1000) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme.myClassic +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        legend.position = c(.8, .8)) +
  scale_fill_manual(values = c("orange", "black"),
                    labels = c("CO2" = "direct oxidation", 
                               "non-Ox sink" = "direct storing")) + # NOX for non-oxidative sink
  labs(x = NULL, y = "µmol C/min/animal\n") +
  coord_cartesian(ylim = c(0, NA))

ggsave(filename = "direct Ox NOS WT.pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       width = 6, height = 4)


# 2) plot stackec bars of OX and NOS

d.flx.CO2.nonOx.sink <- d.flx.CO2.nonOx.sink %>% 
  mutate(Compounds = factor(Compounds, levels = ordered.Compound %>% rev()),
         phenotype = factor(phenotype, levels = ordered.phenotype)) %>% 
  # add position of the error bar
  arrange(phenotype, destiny, desc(Compounds)) %>% 
  mutate(y.error.animal = cumsum(nmol.min.animal.mean)) 

# per animal, all phenotypes
func.plt.stacked.OX.NOS <- function(dataset = d.flx.CO2.nonOx.sink){
  dataset %>% 
    ggplot(aes(x = phenotype, y = nmol.min.animal.mean, fill = Compounds)) +
    geom_col(color = "black") +
    facet_wrap(~destiny, nrow = 1) +
    scale_fill_manual(values = color.Compounds, guide = guide_legend(reverse=T)) +
    geom_errorbar(aes(ymin = y.error.animal - nmol.min.animal.SEM,
                      ymax = y.error.animal), width = .3) +
    scale_y_continuous(expand = expansion(mult = c(0, .1)),
                       breaks = seq(0, 160000, 20000),
                       labels = function(x){x/1000} ) + # to umol /min/animal)
    scale_x_discrete(expand = expansion(add = 1),
                     labels = function(x){str_replace(x, "WT", "control")}) +
    theme.myClassic +
    theme(axis.title.x = element_blank(),
          panel.spacing = unit(0.1, "pt"))  +
    labs(y = "µmol C-atoms / min / animal\n", 
         title = "Nutrients direct flux to CO2 / non-Ox sink") +
    guides(fill = guide_legend(reverse = F))  
}

# - WT only
plt.stacked.CO2.nonOx.sink.WT <- 
  d.flx.CO2.nonOx.sink %>% 
  filter(phenotype == "WT") %>% 
  func.plt.stacked.OX.NOS

plt.stacked.CO2.nonOx.sink.WT

ggsave(filename = "CO2 sink stacked bars WT.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 5, width = 5)

# 3 phenotype
plt.flx.CO2.nonOx.sink.3pheno <- 
  d.flx.CO2.nonOx.sink %>% 
  func.plt.stacked.OX.NOS

plt.flx.CO2.nonOx.sink.3pheno

ggsave(filename = "CO2 sink stacked bars 3 phenotypes.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 5, width = 7)



# mark with numbers

# compare with the total CO2 production (umol/min/animal) measured directly with met cage
load("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/3_core_13CO2_infusion_physiology_basics.RData")


# calcualte mean total CO2 production by met cage direct measurement
d.totalCO2 <- d.energyExpenditure %>% 
  group_by(phenotype) %>% 
  summarise(totalCO2.umol.min = mean(CO2.umol.min))


# calculate pct contribution from each nutrient to total CO2 output
d.pct <- d.flx.CO2.nonOx.sink %>% 
  filter(destiny == "CO2") %>% 
  select(phenotype, destiny, Compounds, nmol.min.animal.mean) %>% 
  left_join(d.totalCO2, by = "phenotype") %>% 
  # for destiny to CO2, calculate % relative to total CO2 output; 
  # for NOS, input NA
  mutate(pct.total = nmol.min.animal.mean/1000/totalCO2.umol.min * 100) 


# label bars with numbers

NNN <- d.pct$Compounds %>% n_distinct()

plt.flx.CO2.nonOx.sink.3pheno +
  # absolute fluxes
  geom_text(aes(label = round(nmol.min.animal.mean/1000, 3)),
            position = position_stack(vjust = .5)) +
  # for CO2, label contribution % of total CO2 measured directly by met cage
  geom_text(data = d.pct,
            aes(x = rep((1:3)+.2, each = NNN), y = nmol.min.animal.mean, 
                label = round(pct.total, 1) %>% paste("%")),
            inherit.aes = F, color = "blue",
            position = position_stack(vjust = .5),
            hjust  = 0, size = 3) +
  
  # sum of all nutrients
  geom_text(data = d.flx.CO2.nonOx.sink %>% group_by(phenotype, destiny) %>% 
              summarise(nmol.min.animal.mean.sum = sum(nmol.min.animal.mean)),
            aes(x = phenotype,
                y = nmol.min.animal.mean.sum + 4000, # move upward by 5000 nmol/min/anim
                label = paste("sum\n", round(nmol.min.animal.mean.sum/1000, 1))),
            inherit.aes = F,
            color = "green2", fontface = "bold") +
  
  # total CO2 output rate directly measured
  # bar
  geom_col(data = SS <- d.pct %>% select(phenotype, totalCO2.umol.min) %>% distinct(),
           aes(x = phenotype, 
               y = totalCO2.umol.min*1000), # convert the total CO2 from umol to nmol
           inherit.aes = F, linetype = "dashed", color = "black", alpha = .05) +
  # text
  geom_text(data = SS,
            aes(x = phenotype, 
                y = totalCO2.umol.min*1000, # convert the total CO2 from umol to nmol
                label = paste("total CO2\n", round(totalCO2.umol.min, 1))),
            fontface = "bold", color = "cyan2",
            inherit.aes = F) +
  
  # adjust y scale
  coord_cartesian(ylim = c(0, 90 * 1000), clip = "off") +
  theme(strip.text = element_text(size = 18))

ggsave(filename = "CO2 sink stacked bars 3 pheno numbered.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 8, width = 12)



d.energyExpenditure %>% 
  ggplot(aes(x = CO2.umol.min, fill = phenotype)) +
  geom_density()

d.energyExpenditure %>% 
  ggplot(aes(x = CO2.umol.min, fill = phenotype)) +
  geom_histogram(alpha = .4, position = "dodge")+
  facet_wrap(~phenotype, nrow = 3)


# WT only
f.plt.WT.CO2 <- function(whichY){
  d.energyExpenditure %>% 
    filter(phenotype == "WT") %>% 
    filter(CO2.umol.min > 35) %>% # remove one outliner in WT
    ggplot(aes(x = phenotype, y = {{whichY}}, fill = phenotype)) + 
    # mean
    stat_summary(fun = mean, geom = "bar", alpha = .7, color = "black") +
    # standard deviation
    stat_summary(fun.data = mean_se, geom = "errorbar",
                 fun.args = list(mult = 1),
                 width = .3) +
    # points
    geom_quasirandom(show.legend = F, alpha = .5, shape = 21) +
    scale_fill_manual(values = color.phenotype) +
    facet_wrap(~phenotype) +
    theme.myClassic +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none") +
    coord_cartesian(clip = "off") +
    scale_x_discrete(expand = expansion(add = .7)) +
    scale_y_continuous(breaks = seq(0, 70, 10),
                       limits = c(0, 70),
                       expand = expansion(mult = c(0, .01))) 
}

# WT per animal
plt.CO2.total.WT <- f.plt.WT.CO2(whichY = CO2.umol.min)
plt.CO2.total.WT 

# WT per g BW
plt.CO2.total.WT.BW <- f.plt.WT.CO2(whichY = CO2.umol.min/ 29 ) +
  scale_y_continuous(
    breaks = seq(0, 3, .3),
    expand = expansion(mult = c(0, .01)),
    name = "CO2 umol / min /g BW") 
plt.CO2.total.WT.BW 



# WT only: combine direct total CO2 with stacked nutrients contribution to CO2

# update the format of stacked bars to be aligned with the total direct CO2 measurement 
plt.flx.CO2.nonOx.sink.animal.WT.2 <- 
  plt.stacked.CO2.nonOx.sink.WT +
  theme.myClassic +
  theme(panel.spacing = unit(0 , "pt"),
        axis.title.x = element_blank(),
        plot.title = element_blank()) +
  coord_cartesian(clip = "off") +
  scale_x_discrete(expand = expansion(add = .7)) +
  scale_y_continuous(breaks = seq(0, 70000, 10000),
                     expand = expansion(mult = c(0, .01)),
                     limits = c(0, 70000),
                     labels = function(x){x/1000}) # convert nmol to µmol

plt.flx.CO2.nonOx.sink.animal.WT.2

plt.CO2.NOS.WT <-  
  plot_grid(
    plt.CO2.total.WT + theme(axis.text.x = element_blank()), 
    # ggplot() + theme_void(),
    plt.flx.CO2.nonOx.sink.animal.WT.2 + theme(axis.title = element_blank(), 
                                               axis.title.x = element_blank(),
                                               axis.text.x = element_blank(),
                                               axis.text.y = element_blank(),
                                               axis.line.y = element_blank(),
                                               axis.ticks.y = element_blank(),
                                               axis.ticks.x = element_blank(),
                                               panel.spacing = unit(4, "pt")),
    rel_widths = c(5, 12), nrow = 1)

plt.CO2.NOS.WT

ggsave(filename = "CO2 sink stacked nutrients WT_gBW.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 4, width = 5)




# total CO2, CO2 nutrient-wise, and non-ox sink
plt.CO2.total.3pheno <- d.energyExpenditure %>% 
  filter(CO2.umol.min < 90) %>% # remove one outliner of ob/ob 
  filter(CO2.umol.min > 35) %>% # remove one outliner in WT
  
  ggplot(aes(x = phenotype, y = CO2.umol.min, fill = phenotype)) + 
  # mean
  stat_summary(fun = mean, geom = "bar", alpha = .7, color = "black") +
  # standard deviation
  stat_summary(fun.data = mean_se, geom = "errorbar",
               fun.args = list(mult = 1),
               width = .3) +
  # points
  geom_quasirandom(show.legend = F, alpha = .5, shape = 21) +
  scale_fill_manual(values = color.phenotype) +
  theme.myClassic +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  scale_x_discrete(expand = expansion(add = 1),
                   labels = function(x){str_replace(x, "WT", "control")}) +
  scale_y_continuous(breaks = seq(0, 120, 25),
                     limits = c(0, 100),
                     expand = expansion(mult = c(0, .1))) 

plt.CO2.total.3pheno

# all 3 phenotypes: combine direct total CO2 with stacked nutrients contribution to CO2
plt.CO2.sink.3pheno <-  plot_grid(
  plt.CO2.total.3pheno, 
  
  plt.flx.CO2.nonOx.sink.3pheno + 
    theme(strip.text = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1),
          plot.title = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_y_continuous(limits = c(0, 100000),
                       breaks = seq(0, 120000, 25000),
                       labels = function(x){x/1000},
                       expand = expansion(mult = c(0, .1))) ,
  
  rel_widths = c(5.5, 12))

plt.CO2.sink.3pheno

ggsave(filename = "CO2 sink stacked nutrients 3 pheno.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 4, width = 6)


# output ordinary infusion parameters
d.tracer.conc <-  d.normalized.tidy %>% 
  group_by(infused.tracer) %>% 
  summarise(tracer.con.mM = mean(tracer_conc_mM)) 

d.ordinary.infusion.info <-  d.normalized.tidy %>% 
  # filter(phenotype != "db/db") %>% 
  group_by(infused.tracer) %>% 
  summarise(uL.g.min.min = (infusion_uL_perMin / BW) %>% min(),
            uL.g.min.max = (infusion_uL_perMin / BW) %>% max(),
            tracer_conc_mM = mean(tracer_conc_mM)) %>% 
  mutate(uL.g.min = paste(uL.g.min.min, "-", uL.g.min.max)) %>% 
  select(infused.tracer, tracer_conc_mM, uL.g.min) %>% 
  rename(tracer = infused.tracer)

# combine with the 13CO2 infusion and bolus information
d.tracer.info.all <-  d.13C.info %>% left_join(d.ordinary.infusion.info)



# Output dataset
file.supp.table <- "/Users/boyuan/Desktop/Harvard/Research/db db mice/Infusion data/A_V_ratio_publish_formatted.xlsx"

write.xlsx(x = d.tracer.info.all, 
           file = file.supp.table, sheetName = "tracer info")


# output the A-V ratio dataset
# mean and SEM
write.xlsx(x = d.art.vs.tail.ratio.atom.summary.paper %>% 
             mutate(phenotype = factor(phenotype, levels = ordered.phenotype)) %>% 
             arrange(infused.tracer, phenotype) %>% as.data.frame(),
           file = file.supp.table,
           sheetName = "A_V_ratio", append = T)




# Compare the overall and direct oxidation flux
x1 <- d.ox.sink.overal.2 %>% 
  filter(phenotype == "WT" & destiny == "ox.overal") %>% 
  ungroup() %>% 
  select(Compound, nmol.min.animal, nmol.min.animal.SEM) %>% 
  mutate(type = "overall")

x2 <- d.flx.CO2.nonOx.sink %>% 
  filter(phenotype == "WT" & destiny == "CO2") %>% 
  ungroup() %>% 
  select(Compounds, nmol.min.animal.mean, nmol.min.animal.SEM) %>% 
  mutate(type = "direct") %>% 
  rename(nmol.min.animal = nmol.min.animal.mean, Compound = Compounds)

bind_rows(x1, x2) %>% 
  mutate(Compound = factor(Compound, levels = u$Compound)) %>% 
  ggplot(aes(x = Compound, y = nmol.min.animal, fill = type)) +
  geom_col(color = "black", width = .7, alpha = .5, position = position_dodge(.7)) +
  geom_errorbar(aes(ymax = nmol.min.animal + nmol.min.animal.SEM,
                    ymin = nmol.min.animal - nmol.min.animal.SEM),
                width = .2,
                position = position_dodge(.7)) +
  scale_y_continuous(breaks = seq(0, 16000, 4000),
                     position = "left",
                     labels = function(x){x/1000},
                     name = "oxidation (µmol C / min)",
                     expand = expansion(mult = c(0, .1))) +
  scale_fill_manual(values = c("orange", "black")) + # NOX for non-oxidative sink
  theme.myClassic +
  theme(
    legend.text = element_text(size = 15),
    legend.key.size = unit(20, "pt"),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(.84, .8),
    axis.title.y = element_text(margin = margin(r = 10, "pt")),
    plot.margin = margin(l = 20)) +
  labs(x = NULL) +
  coord_cartesian(ylim = c(0, NA))


ggsave(
  filename = "direct and overall oxidation comparison.pdf",
  path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures", 
  width = 6, height = 3.5)

# save the project
save.image(file = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/6_core_consumption_flux.RData")
