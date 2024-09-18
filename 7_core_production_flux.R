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

load("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/6_core_consumption_flux.RData")


# Interconversions and production fluxes ----------
# Direct contributions fractions: Using Tony's method (Cell metabolism 2020)

d.contribute.fraction = tibble()
d.monteCarloRecord = tibble()

for (destinationCompound in myCompounds ){
  for (destinationPhenotype in ordered.phenotype[1:3]){
    
    # destinationCompound = "Glucose"
    # destinationPhenotype = "WT"
    
    d.i = d.enrich.atom.normalized.summary %>% ungroup() %>%
      # select phenotype
      filter(phenotype == destinationPhenotype) %>%
      # remove the Compound whose contributing source is calculated
      filter(infused.tracer != destinationCompound & Compound != destinationCompound) %>%
      select(-phenotype)
    
    # mean matrix
    d.i.mean = d.i %>%
      # create matrix for mean or sd
      select(infused.tracer, Compound, enrich.normalized.mean) %>%
      spread(key = Compound, value = enrich.normalized.mean)
    # change to matrix and add row names
    mat.mean = d.i.mean %>% select(-infused.tracer) %>% as.matrix()
    rownames(mat.mean) = d.i.mean$infused.tracer
    
    # SEM matrix
    d.i.SEM = d.i %>%
      # create matrix for mean or sd
      select(infused.tracer, Compound, enrich.normalized.sem) %>%
      spread(key = Compound, value = enrich.normalized.sem)
    # change to matrix and add row names
    mat.SEM = d.i.SEM %>% select(-infused.tracer) %>% as.matrix()
    rownames(mat.SEM) = d.i.SEM$infused.tracer
    
    # check row names match
    if (((! (rownames(mat.mean) == rownames(mat.SEM))) %>% sum()) > 1) stop("Matrices Compound order does not match")
    
    
    # vector : destinationCompound normalized enrichment
    d.destinationCompound = d.enrich.atom.normalized.summary %>% ungroup() %>%
      # select phenotype
      filter(phenotype == destinationPhenotype) %>%
      # remove the Compound whose contributing source is calculated
      filter((infused.tracer != destinationCompound) & (Compound == destinationCompound)) %>%
      select(-c(phenotype, Compound))
    
    # vector of mean
    v.mean = d.destinationCompound$enrich.normalized.mean
    names(v.mean) = d.destinationCompound$infused.tracer
    # vector of standard deviation
    v.SEM = d.destinationCompound$enrich.normalized.sem
    names(v.SEM) = d.destinationCompound$infused.tracer
    
    # check matrix Compound order match vector Compound order
    if (((! (rownames(mat.mean) == names(v.mean))) %>% sum()) > 1)stop(
      "Matrices Compound order does not match vector Compound order.")
    
    
    # Monte Carlo simulation
    mat.contributionFraction.repeat = vector()
    for (i in 1:100){
      mat.standardGaussain = rnorm(n = nrow(mat.mean)^2, mean = 0, sd = 1) %>% matrix(nrow = nrow(mat.mean))
      mat.MonteCarlo =  mat.mean + mat.SEM * mat.standardGaussain
      
      v.standardGaussain = rnorm(n = nrow(mat.mean), mean = 0, sd = 1)
      v.MonteCarlo = v.mean + v.SEM * v.standardGaussain
      
      # Apply non-negative constraint; use lsei package; this is the applied method in this publication!!
      # citation: Package limSolve , solving linear inverse models in R, by Karline Soetaert et al.
      list.contributionFraction.i = lsei(A = mat.MonteCarlo, B = v.MonteCarlo, # this solves Ax = b
                                         G = diag(1, nrow = length(v.MonteCarlo)), # G is the identity matrix, dimension = B-vector length
                                         H = rep(0, length(v.MonteCarlo)), # vector of zero, length = B-vector length
                                         verbose = FALSE)
      v.contributionFraction.i = list.contributionFraction.i$X
      
      mat.contributionFraction.repeat = cbind(mat.contributionFraction.repeat, v.contributionFraction.i)
    }
    
    # a record of Monte Carlo each round of simulation
    d.monteCarloRecord.i = mat.contributionFraction.repeat %>% as_tibble() %>%
      mutate(phenotype = destinationPhenotype,
             destinationCompound = destinationCompound,
             originCompound = rownames(mat.contributionFraction.repeat)) %>%
      select(phenotype, destinationCompound, originCompound,  1:ncol(mat.contributionFraction.repeat))
    d.monteCarloRecord = rbind(d.monteCarloRecord, d.monteCarloRecord.i)
    
    
    # summary result from Monte Carlo simulation
    d.contribute.fraction.i =
      tibble(OriginCompound = rownames(mat.mean),
             contribute.fraction.mean = apply(mat.contributionFraction.repeat,  MARGIN = 1, mean),
             contribute.fraction.SEM = apply(mat.contributionFraction.repeat,  MARGIN = 1, sd)) %>%
      mutate(phenotype = destinationPhenotype,
             destinationCompound = destinationCompound) %>%
      select(phenotype, destinationCompound, OriginCompound, contains("fraction")) %>%
      mutate(totalCounted = sum(contribute.fraction.mean),
             totalCounted.SEM = sum(contribute.fraction.SEM^2) %>% sqrt() )
    d.contribute.fraction = rbind(d.contribute.fraction, d.contribute.fraction.i)
  }
}



# Compound as ordered factor
d.contribute.fraction = d.contribute.fraction %>%
  mutate(destinationCompound = factor(destinationCompound, levels = ordered.Compound, ordered = F),
         OriginCompound = factor(OriginCompound, levels = ordered.Compound, ordered = F))

# mark y-axis position for error bars
d.contribute.fraction = d.contribute.fraction %>%
  group_by(phenotype, destinationCompound) %>%
  mutate(yaxis.errorBar = func.cumulatedSum(contribute.fraction.mean))


# Plot
plt.labeling.directContribution = d.contribute.fraction %>%
  mutate(phenotype = factor(phenotype, levels = ordered.phenotype)) %>% 
  ggplot(aes(x = phenotype, y = contribute.fraction.mean,
             fill = OriginCompound, color = OriginCompound)) +
  geom_bar(stat = "identity",  color = "black") +
  geom_errorbar(aes(ymin = yaxis.errorBar - contribute.fraction.SEM,
                    ymax = yaxis.errorBar ),
                width = .2, show.legend = F, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 1, .1)) +
  scale_x_discrete(expand = expansion(add = 0),
                   labels = function(x){str_replace(x, "WT", "control")}) +
  labs(y = "Direct contribution fraction\n",
       title = "Direct contribution fraction, mean ± SEM\nUsing Tony's method (Cell Metabolism 2020)\n ") +
  theme.myClassic +
  theme(axis.title.x = element_blank()) +
  # scale_fill_nejm() +scale_color_nejm() +
  facet_wrap(~destinationCompound, nrow = 2) +
  scale_color_manual(values = color.Compounds) +
  scale_fill_manual(values = color.Compounds) +
  # add y = 1 bar border
  geom_col(data = NULL, aes(y = 1), fill = NA, color = "black", 
           position = "identity", linetype = "longdash", linewidth = .5) 

plt.labeling.directContribution




# Calculate fluxes directly with refined method  ---------



# Calculate circulating nutrients inter-conversion fluxes
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
        H = rep(0, matrixDimension), # vector of zero, length = B-vector length
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
  mutate(nmol.min.g.SEM.Y.axis.position = func.cumulatedSum(nmol.min.g.mean)) 


# Calculate per animal fluxes
d.directFlux_interConversion = d.directFlux_interConversion %>%  # good
  # combine with body weight data
  left_join(d.BW %>% rename(targetCompound = Compounds), 
            by = c("targetCompound", "phenotype")) %>% 
  # calculate animal-based flux
  mutate(nmol.min.animal.mean = nmol.min.g.mean * BW.mean,
         nmol.min.animal.SEM =  nmol.min.g.SEM * BW.mean,
         nmol.min.animal.SEM.Y.axis.position = func.cumulatedSum(nmol.min.animal.mean)) %>% 
  mutate(targetCompound = factor(targetCompound, levels = ordered.Compound),
         phenotype = factor(phenotype, levels = ordered.phenotype))


# function to apply theme to inter-convertion flux
func.plot.theme = function(plot){
  plot + 
    facet_wrap(~targetCompound, scales = "free", nrow = 2) +
    # scale_color_manual(values = color.Compounds) +
    scale_fill_manual(values = color.Compounds, guide = guide_legend(reverse=T)) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    expand_limits(y = 0) +
    scale_x_discrete(expand = expansion(mult = c(.5, .5)),
                     labels = function(x){str_replace(x, "WT", "control")}) +
    theme.myClassic +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          panel.spacing = unit(20, "pt") ) +
    guides(fill = guide_legend(reverse = F))
}

# per g BW
plt.flx.interConversion.g.BW = 
  d.directFlux_interConversion %>% 
  ggplot(aes(x = phenotype, y = nmol.min.g.mean, fill = sources)) +
  geom_bar(stat = "identity", color = "black", size = .4, alpha = .9) +
  geom_errorbar(aes(ymin = nmol.min.g.SEM.Y.axis.position - nmol.min.g.SEM,
                    ymax = nmol.min.g.SEM.Y.axis.position),
                width = .3, size = .5) 
plt.flx.interConversion.g.BW = func.plot.theme(plt.flx.interConversion.g.BW) +
  labs(title = "direct sources of circulating nutrients", y = " nmol C-atoms / g BW / min")
plt.flx.interConversion.g.BW



# per animal - plot1: WT only
plt.source.WT <- d.directFlux_interConversion %>% 
  filter(phenotype == "WT") %>% 
  ungroup() %>%  # Critical!!!! the reordered column cannot be a grouping variable!!! 
  mutate(targetCompound = as.character(targetCompound)) %>% 
  mutate(targetCompound = fct_reorder(targetCompound, nmol.min.animal.mean, 
                                      .fun = sum, .desc = T)) %>% 
  ggplot(aes(x = targetCompound, y = nmol.min.animal.mean, fill = sources)) +
  geom_col(color = "black", size = .4, alpha = .9) +
  geom_errorbar(aes(ymin = nmol.min.animal.SEM.Y.axis.position - nmol.min.animal.SEM,
                    ymax = nmol.min.animal.SEM.Y.axis.position),
                width = .3, size = .5) +
  scale_fill_manual(values = color.Compounds, guide = guide_legend(reverse=T)) +
  scale_x_discrete(expand = expansion(add = .8)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)),
                     n.breaks = 8,
                     labels = function(x){x / 1000}) +
  expand_limits(y = 0) +
  theme.myClassic +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(r =10, "pt"), size = 18),
        panel.spacing = unit(20, "pt"),
        plot.title = element_blank(),
        legend.position = c(.86, .65)) +
  guides(fill = guide_legend(reverse = F)) +
  labs(title = "direct sources of circulating nutrients\n", 
       y = " µmol C-atoms / animal / min)") 
plt.source.WT

ggsave(
  filename = "production source flux.WT.pdf",
  path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures", 
  width = 6.3, height = 5.5)


# normalize by body weight --------
plt.source.WT +
  scale_y_continuous(
    expand = expansion(mult = c(0, .1)),
    n.breaks = 8,
    name = "C umol / min /animal",
    position = "right",
    labels = function(x){x / 1000},
    sec.axis = sec_axis(trans = ~ . / 29, breaks = seq(0, 1000, 100),
                        name = "C nmol / min / g BW")) 

ggsave(
  filename = "production source flux.WT_gBW.pdf",
  path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures", 
  width = 6.6, height = 5.5)





# per animal - plot 2: all phenotypes
plt.flx.interConversion.animal = d.directFlux_interConversion %>% 
  #filter(phenotype == "WT") %>% 
  ggplot(aes(x = phenotype, y = nmol.min.animal.mean, fill = sources)) +
  geom_col(color = "black", size = .4, alpha = .9) +
  geom_errorbar(aes(ymin = nmol.min.animal.SEM.Y.axis.position - nmol.min.animal.SEM,
                    ymax = nmol.min.animal.SEM.Y.axis.position),
                width = .3, size = .5)

plt.flx.interConversion.animal = func.plot.theme(plt.flx.interConversion.animal) +
  labs(title = "direct sources of circulating nutrients\n", 
       y = " nmol C-atoms / animal / min)")

plt.flx.interConversion.animal 




# accumulated circulatory atom 
d.totalProduction.nutrientWise <-  d.directFlux_interConversion %>%  
  select(phenotype, targetCompound, sources, 
         nmol.min.animal.mean, nmol.min.animal.SEM)  %>% 
  group_by(phenotype, targetCompound) %>% 
  summarise(nmol.min.animal.mean.Fcirc.atom = sum(nmol.min.animal.mean),
            nmol.min.animal.SEM.Fcirc.atom = sum(nmol.min.animal.SEM^2) %>% sqrt())

plt.totalProduction.atom.3pheno <-  
  d.totalProduction.nutrientWise %>% 
  arrange(phenotype, targetCompound) %>% 
  mutate(targetCompound = factor(targetCompound, levels = ordered.Compound %>% rev(), ordered = F)) %>% 
  arrange(targetCompound) %>% 
  mutate(SEM.yAxis = func.cumulatedSum(nmol.min.animal.mean.Fcirc.atom)) %>% 
  mutate(phenotype = factor(phenotype, levels = ordered.phenotype)) %>% 
  
  ggplot(aes(x = phenotype, y = nmol.min.animal.mean.Fcirc.atom, fill = targetCompound)) + 
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = SEM.yAxis - nmol.min.animal.SEM.Fcirc.atom,
                    ymax = SEM.yAxis),
                width = .2) +
  scale_fill_manual(values = color.Compounds) +
  theme.myClassic +
  scale_x_discrete(expand = expansion(add = .8),
                   labels = function(x){str_replace(x, "WT", "control")}) +
  scale_y_continuous(limits = c(0, NA),
                     breaks = seq(0, 200000, 25000),
                     expand = expansion(mult = c(0, .1)) , 
                     labels = function(x){x / 1000}) + # here convert the unit from nmol to µmol
  labs(y = "µmol C atoms / min / animal\n",
       title = "Total production flux of circulating nutrients")

plt.totalProduction.atom.3pheno

ggsave(
  filename = "total production 3 pheno.pdf",
  path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures", 
  width = 5.5, height = 5)






# Mark the production with numbers
d.pct <- d.directFlux_interConversion %>% 
  group_by(phenotype, targetCompound) %>% 
  mutate(contr.fraction = (nmol.min.animal.mean / sum(nmol.min.animal.mean)) %>% round(2) )

plt.flx.interConversion.animal +
  # label with contribution fraction
  geom_text(data = d.pct,
            # if contribution  < 5%, not label it
            aes(label = ifelse(contr.fraction > .05, paste(contr.fraction * 100, "%"), NA)),
            position = position_stack(vjust = .4), 
            size = 2, color = "blue") +
  # absolute source flux
  geom_text(data = d.pct,
            # if flux smaller than 1000, not label it
            aes(label = ifelse(nmol.min.animal.mean > 1000, round(nmol.min.animal.mean), NA)),
            position = position_stack(vjust = .6),
            color = "black", 
            size = 2)  +
  # total production
  geom_text(
    data = d.totalProduction.nutrientWise,
    aes(x = phenotype, 
        y = nmol.min.animal.mean.Fcirc.atom * 1.05,
        label = nmol.min.animal.mean.Fcirc.atom %>% round()),
    inherit.aes = F, 
    size = 3, fontface = "bold", color = "cyan3"
  )

ggsave(
  filename = "production source flux_numbers.pdf",
  path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures", 
  width = 12, height = 10)


# display the fraction percentage
d.directFlux_interConversion %>% 
  filter(targetCompound %in% c("Glucose", "Lactate")) %>% 
  ggplot(aes(x = phenotype, y = nmol.min.animal.mean, fill = sources)) +
  geom_col(position = "fill", color = "black") +
  scale_fill_manual(values = color.Compounds) +
  scale_x_discrete(labels = function(x){str_replace(x, "WT", "control")}) +
  facet_wrap(~targetCompound) +
  theme.myClassic 



# glycerol quick math
# glycerol release from TG: 5441 C / min  = 1813 glycerol or TG molecules nmol/min
# total FFA released expected = 1813 * 3 = 5439 FFA
# Total TG molecules :
# (10896/16 + 18579/18 + 26122/18)  = 3163 molecules / min

(1813 * 3 - 3163) / 5439

# Compare the 3 ways of Fcirc calculation

# calculate atom Fcirc by molecule Fcirc times carbon number
d.Fcirc.timesCarbonNumber = d.Fcirc_standard.atom.summary %>% 
  select(Compound, phenotype, Fcirc_g.BW.mean, Fcirc_g.BW.sd, n.rep) %>% 
  # add total carbon number of each metabolite
  left_join(d.normalized.tidy %>% select(Compound, C_Label.max) %>% distinct()) %>% 
  mutate(Fcirc_gBW.timesCarbonNumber.mean = Fcirc_g.BW.mean * C_Label.max,
         Fcirc_gBW.timesCarbonNumber.sem = Fcirc_g.BW.sd * C_Label.max / sqrt(n.rep)) %>% 
  ungroup() %>% 
  select(phenotype, Compound, contains("timesCarbonNumber")) %>% 
  rename(targetCompound = Compound,
         nmol.min.gBW = Fcirc_gBW.timesCarbonNumber.mean,
         nmol.min.g.SEM = Fcirc_gBW.timesCarbonNumber.sem) %>% 
  mutate(method = "Fcirc.molecule x C #")
d.Fcirc.timesCarbonNumber


# function: calculate SD from summed elements
func.square.sum.sqrt = function(vector){vector^2 %>% sum() %>% sqrt()}

d.Firc.linearEq = d.directFlux_interConversion %>% 
  select(phenotype, targetCompound, sources, contains("min.g")) %>% 
  select(-contains("position")) %>% 
  group_by(phenotype, targetCompound) %>% 
  summarise(nmol.min.gBW = sum(nmol.min.g.mean),
            nmol.min.g.SEM = func.square.sum.sqrt(nmol.min.g.SEM)) %>% 
  mutate(method = "linear Eq. sets")

d.Ficrc.atom.methodComparison = d.Firc.linearEq %>% rbind(d.Fcirc.timesCarbonNumber)  %>% 
  ungroup() %>% 
  left_join(d.BW %>% rename(targetCompound = Compounds), 
            by = c("phenotype", "targetCompound")) %>%
  mutate(nmol.min.animal = nmol.min.gBW * BW.mean,
         nmol.min.animal.SEM = nmol.min.g.SEM * BW.mean)


# plot : calculate the relative difference of the two methods
d.Ficrc.atom.methodComparison.spread <- 
  d.Ficrc.atom.methodComparison %>% 
  filter(phenotype != "db/db") %>% 
  filter(targetCompound %in% ordered.Compound[c(1:6, 8:11)]) %>% 
  select(phenotype, targetCompound, nmol.min.gBW, method) %>% 
  spread(key = method, value = nmol.min.gBW) %>% 
  mutate(diff.pct = (`Fcirc.molecule x C #` - `linear Eq. sets`) / `linear Eq. sets` * 100) %>% 
  mutate(targetCompound = fct_reorder(targetCompound, diff.pct, .fun = mean,  .desc = T)) 

d.Ficrc.atom.methodComparison.spread %>% 
  ggplot(aes(x = targetCompound, y = diff.pct, fill = phenotype)) +
  geom_col(position = position_dodge(.7), color = "black",
           width = .7) +
  coord_flip() +
  scale_fill_manual(values = color.phenotype,
                    labels = function(x){str_replace(x, "WT", "control")}) +
  geom_hline(yintercept = 0) +
  theme.myClassic +
  labs(y = "(`Fcirc.molecule x C #` - `linear Eq. sets`) / `linear Eq. sets` * 100)",
       x = NULL)


# 
# d.diff.pct <- d.Ficrc.atom.methodComparison.spread %>%
#   select(phenotype, targetCompound, diff.pct)
# 
# 
# d.diff.pct.enrich <- d.enrich.atom.summary %>%
#   filter(infused.tracer == Compound) %>%
#   filter(phenotype != "db/db") %>%
#   select(phenotype, infused.tracer, enrich.original.mean, enrich.original.sd) %>%
#   rename(targetCompound = infused.tracer) %>%
#   left_join(d.diff.pct)
# #  
# # 
# d.diff.pct.enrich %>%
#   # filter(targetCompound %in% c("Valine", "Glutamine")) %>%
#   ggplot(aes(x = enrich.original.sd, diff.pct, color = phenotype)) +
#   geom_text(aes(label = targetCompound)) +
#   scale_color_manual(values = color.phenotype) +
#   facet_wrap(~targetCompound, scales = "free")
#  
# d.diff.pct.enrich %>%
#   # filter(targetCompound %in% c("Valine", "Glutamine")) %>%
#   ggplot(aes(x = enrich.original.mean, diff.pct, color = phenotype)) +
#   geom_text(aes(label = targetCompound)) +
#   scale_color_manual(values = color.phenotype)+
#   facet_wrap(~targetCompound, scales = "free")
# 
# lm(data = d.diff.pct.enrich, diff.pct ~ enrich.original.mean + enrich.original.sd) %>% summary()



# Extract the atom Fcirc 
d.Fcirc.atom.linear.Eq = d.Ficrc.atom.methodComparison %>% 
  filter(method == "linear Eq. sets") %>% 
  mutate(targetCompound = factor(targetCompound, levels = ordered.Compound, ordered = F)) %>% 
  arrange(targetCompound) %>%
  as_tibble()


# manual check
# d.directFlux_interConversion %>% filter(targetCompound == "Valine" ) %>% view()





# Compare overall OX + NOX (Fcirc atom) with total production
# total production 
x1 <- d.directFlux_interConversion %>% 
  group_by(phenotype, targetCompound) %>% 
  summarise(totalP = sum(nmol.min.animal.mean)) %>% 
  rename(Compound = targetCompound)

# Fcirc atom (sum of overall OX + NOS flux)
x2 <- d.ox.sink.overal.2 %>% 
  group_by(phenotype, Compound) %>% 
  summarise(Fcirc.atom = sum(nmol.min.animal))

# combine and calculate the difference between total P and Fcirc atom
x2 %>% left_join(x1) %>% 
  filter(phenotype != "db/db") %>% 
  mutate(delta.pct = (totalP - Fcirc.atom)/totalP * 100) %>% 
  ggplot(aes(Compound, y = delta.pct, fill = phenotype)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = color.phenotype,
                    labels = function(x){str_replace(x, "WT", "control")}) +
  coord_flip() +
  labs(y = "(Total P - Fcirc.atom) / total P × 100%")







# Nutrient metabolites destiny
d.directFlux_interConversion.Renamed = 
  d.directFlux_interConversion %>% 
  # think of sources as the target compounds, and original target as the destination of the sources
  rename(destiny = targetCompound, Compounds = sources) %>% 
  select(phenotype, destiny, Compounds, contains("nmol.min")) %>% select(-contains("axis"))



# dataset of the direct destiny of each nutrients, including conversion to other nutrients, to CO2 and to non-Ox sink
d.destiny =  d.flx.CO2.nonOx.sink %>% 
  select(phenotype, destiny, Compounds, contains("nmol.min")) %>% 
  mutate(destiny = factor(destiny, levels = ordered.Compound, ordered = F)) %>% 
  
  # combine with direct inter-conversion flux data between circulating nutrients (excluding storage)
  rbind(
    d.directFlux_interConversion.Renamed %>% filter(Compounds != "storage") %>% 
      mutate(destiny = factor(destiny, levels = ordered.Compound, ordered = F),
             Compounds = factor(Compounds, levels = ordered.Compound, ordered = F))
  ) %>% 
  
  mutate(destiny = factor(destiny, levels = ordered.Compound %>% rev(), ordered = F)) %>% 
  arrange(destiny, Compounds)

d.destiny



# combine with Fcirc atom value solved from linear equations

d.totalProduction <-  d.Fcirc.atom.linear.Eq %>% 
  select(-c(method, BW.mean, BW.sd)) %>% 
  rename(Compounds = targetCompound) %>% 
  rename(P.nmol.min.g = nmol.min.gBW, # total production rate
         P.nmol.min.g.SEM = nmol.min.g.SEM,
         P.nmol.min.animal = nmol.min.animal,
         P.nmol.min.animal.SEM = nmol.min.animal.SEM)

d.destiny$destiny %>% unique()
d.destiny$Compounds %>% unique()

# join 
d.destiny.totalProduction <- d.destiny %>% 
  left_join(d.totalProduction, by = c("phenotype", "Compounds"))  %>% 
  
  # Calculate y-axis position for error bar
  arrange(Compounds, destiny) %>% 
  group_by(phenotype, Compounds) %>% 
  mutate(errorPosition.gBW = func.cumulatedSum(nmol.min.g.mean),
         errorPosition.animal = func.cumulatedSum(nmol.min.animal.mean)) %>% 
  mutate(Compounds = factor(Compounds, levels = ordered.Compound),
         phenotype = factor(phenotype, levels = ordered.phenotype))


# plot
# per g BW
plt.destiny.g.BW = d.destiny.totalProduction %>% 
  
  ggplot(aes(x = phenotype, y = nmol.min.g.mean, fill = destiny)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = errorPosition.gBW - nmol.min.g.SEM,
                    ymax = errorPosition.gBW),width = .3) +
  facet_wrap(~Compounds, scales = "free", nrow = 2) +
  
  theme.myClassic + theme(axis.title.x = element_blank()) +
  scale_fill_manual(values = color.Compounds) +
  scale_x_discrete(expand = expansion(add = 1),
                   labels = function(x){str_replace(x, "WT", "control")}) +
  labs(title = "Fate of circulating nutrients") +
  guides(fill = guide_legend(title = "Direct destiny",
                             title.theme = element_text(face  = "bold"), reverse = T)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  # add Fcirc level
  geom_point(aes(y = P.nmol.min.g), show.legend = F,  color = "grey") +
  labs(y = "nmol C-atoms / min / g BW") +
  guides(fill = guide_legend(reverse = F)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10))

plt.destiny.g.BW



# per animal - plot 1: WT
plt.WT.consumption <-  d.destiny.totalProduction %>% 
  filter(phenotype == "WT") %>% ungroup() %>% 
  mutate(Compounds = fct_reorder(Compounds, nmol.min.animal.mean, .fun = sum, .desc = T)) %>% 
  
  # due to small analytically error, the compounds are ordered slightly different from the production flux 
  # i.e., glucose and C16:0 is swapped
  # here glucose is manually positioned behind C16:0
  mutate(Compounds = fct_relevel(Compounds, "Glucose", after = 7)) %>% 
  
  ggplot(aes(x = Compounds, y = nmol.min.animal.mean, fill = destiny)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = errorPosition.animal - nmol.min.animal.SEM,
                    ymax = errorPosition.animal),width = .3) +
  
  theme.myClassic + theme(axis.title.x = element_blank()) +
  scale_fill_manual(values = color.Compounds,
                    labels = c("non-Ox sink" = "Storage")) +
  labs(title = "Fate of circulating nutrients") +
  guides(fill = guide_legend(title = "Direct destiny",
                             title.theme = element_text(face  = "bold"), reverse = T)) +
  scale_x_discrete(expand = expansion(add = .8)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)),
                     n.breaks = 8,
                     labels = function(x){x / 1000}) +
  # add Fcirc level
  labs(y = "µmol C-atoms / min / animal\n") +
  geom_point(aes(y = P.nmol.min.animal), show.legend = F,  color = "grey") +
  
  guides(fill = guide_legend(reverse = F)) +
  theme.myClassic +
  theme(axis.title.x = element_blank(),
        axis.title.y.left = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(margin = margin(r =10, "pt")),
        panel.spacing = unit(20, "pt"),
        plot.title = element_blank(),
        legend.position = c(.86, .65)) 

plt.WT.consumption

ggsave(filename = "consumption fate WT.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       width = 6.3, height = 5.5)



# normalize WT flux by body weight -_-_-_-_-_-_-_-_-_-_-_-_-_
plt.WT.consumption +
  scale_y_continuous(
    expand = expansion(mult = c(0, .1)),
    n.breaks = 8,
    position = "right",
    labels = function(x){x / 1000},
    sec.axis = sec_axis(trans = ~ . / 29, breaks = seq(0, 1000, 100),
                        name = "C nmol / min / g BW")) +
  coord_cartesian(ylim = c(0, 26000))

ggsave(filename = "consumption fate WT_gBW.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       width = 7, height = 5.5)




# per animal - plot 2: all phenotypes
plt.destiny.animal = d.destiny.totalProduction %>% 
  
  ggplot(aes(x = phenotype, y = nmol.min.animal.mean, fill = destiny)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  geom_errorbar(aes(ymin = errorPosition.animal - nmol.min.animal.SEM,
                    ymax = errorPosition.animal),width = .3) +
  # add text: relative fraction
  geom_text(data = (u <- d.destiny.totalProduction %>% mutate(frac = nmol.min.animal.mean / sum(nmol.min.animal.mean))),
            aes(label = ifelse(frac < .01, NA, round(frac*100, 1) %>% paste("%"))),
            position = position_stack(vjust = .6), size = 2.5, color = "cyan") +
  # add text: absolute flux
  geom_text(data = u,
            aes(label = ifelse(frac < .01, NA, round(nmol.min.animal.mean))),
            position = position_stack(vjust = .4), size = 2.5, color = "green") +
  #
  facet_wrap(~Compounds, scales = "free", nrow = 2) +
  
  theme.myClassic + theme(axis.title.x = element_blank()) +
  scale_fill_manual(values = color.Compounds, labels = c("non-Ox sink" = "Storage")) +
  labs(title = "Fate of circulating nutrients") +
  guides(fill = guide_legend(title = "Direct destiny",
                             title.theme = element_text(face  = "bold"), reverse = T)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)),
                     labels = function(x){x/1000} ) + # to umol /min/animal
  # add Fcirc level
  labs(y = "µmol C-atoms / min / animal\n") +
  geom_point(aes(y = P.nmol.min.animal), show.legend = F,  color = "grey") +
  
  guides(fill = guide_legend(reverse = F)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
        panel.spacing = unit(20, "pt")) +
  scale_x_discrete(expand = c(.4, 0)) 

plt.destiny.animal


ggsave(filename = "destiny numbered.pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 10, width = 15)



# # - <>-trouble shoot abberantly high glycerol enrichment under palmitate infusion
# m <- d.enrich.atom %>% 
#   filter( str_detect(infused.tracer,":" ) & (Compound %in% c("Glycerol", "C16:0", "C18:1", "C18:2")) ) %>% 
#   select(Compound, infused.tracer, enrich.atom.corrected, infusion_mouseID, phenotype, infusion_round)
# 
# d.tracerEnrich <- m %>%  
#   mutate(istracer = Compound == infused.tracer) %>% 
#   filter(istracer == T) %>% 
#   select(-istracer) %>% 
#   rename(tracer.enrich =  enrich.atom.corrected) %>% 
#   select(infusion_mouseID, tracer.enrich)
# 
# 
# 
# (ddd <- m %>% left_join(d.tracerEnrich) %>% 
#     filter(Compound != infused.tracer ) %>% 
#     filter(phenotype != "db/db")) %>% 
#   # filter(infused.tracer == "C18:1") %>% 
#   ggplot(aes(x = tracer.enrich, y = enrich.atom.corrected, color = Compound)) +
#   geom_point() +
#   geom_text_repel(aes(label = infusion_mouseID)) +
#   facet_grid(phenotype~infused.tracer) +
#   labs(y = "Glycerol labeling") 
# # geom_text_repel(# data = ddd %>% filter(enrich.atom.corrected > .03),
# #           aes(label = infusion_round))
# 
# 


d1 <- d.directFlux_interConversion %>% 
  select(phenotype, targetCompound, sources, nmol.min.animal.mean, nmol.min.animal.SEM) %>% 
  rename(to = targetCompound, from = sources) %>% 
  mutate(index = "production",
         to = factor(to, ordered = F),
         from = factor(from, ordered = F),
         x.axis = str_c(phenotype, "_production"),
         myfill = from,
         about = to) 


d2 <- d.destiny.totalProduction %>% 
  select(phenotype, destiny, Compounds, nmol.min.animal.mean, nmol.min.animal.SEM) %>% 
  rename(to = destiny, from = Compounds) %>% 
  mutate(index = "consumption",
         x.axis = str_c(phenotype, "_consumption"),
         myfill = to, 
         about = from)

d.production.consumption <-  rbind(d1, d2) %>% 
  mutate(
    phenotype = factor(phenotype, levels = ordered.phenotype),
    myfill = factor(myfill, levels = ordered.Compound %>% rev(), ordered = T),
    x.axis = factor(x.axis, 
                    levels = c("WT_production", "HFD_production", "ob/ob_production",
                               "WT_consumption", "HFD_consumption", "ob/ob_consumption"))) %>% 
  # add error position  
  arrange(about, x.axis, desc(myfill)) %>% 
  ungroup() %>% 
  group_by(about, x.axis) %>% 
  mutate(err.Y = cumsum(nmol.min.animal.mean))

plt.production.consumption <- 
  d.production.consumption %>% 
  ggplot(aes(x = x.axis, y = nmol.min.animal.mean, fill = myfill)) +
  geom_col(color = "black") +
  facet_wrap(~ about, scales = "free", ncol = 4) +
  
  scale_fill_manual(
    values = color.Compounds, 
    name = "source of production, \nor destiny of consumption",
    labels = c("non-Ox sink" = "Storage (recycle)",
               "storage" = "Storage (release)")) +
  
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  
  scale_x_discrete(
    labels = function(x){
      x %>% str_remove(or("_production", "_consumption")) %>% 
        str_replace("WT", "control")},
    limits = c("WT_production", "HFD_production", "ob/ob_production", "",
               "WT_consumption", "HFD_consumption", "ob/ob_consumption"),
    expand = expansion(add = 1), 
    name = NULL)  +
  
  geom_errorbar(aes(ymax = err.Y,
                    ymin = err.Y - nmol.min.animal.SEM),
                width = .4) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(color = "black"),
        panel.spacing = unit(20, "pt"),
        legend.title = element_text(size = 10)
  ) 

plt.production.consumption

ggsave(filename = "production_consumption.pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 9, width = 15)


# visualize the production and consumption as mirror images 

d.production.consumption.mirror <- 
  d.production.consumption %>% 
  mutate(nmol.min.animal.mean = ifelse(str_detect(x.axis, "consumption"), 
                                       -nmol.min.animal.mean, nmol.min.animal.mean),
         err.Y.animal = ifelse(str_detect(x.axis, "consumption"), 
                               -err.Y, err.Y)) %>% 
  select(-err.Y) %>% # later need to specify the error to basis of animal, gBW, gFat, or gLean
  mutate(x.axis = str_remove(x.axis, or( "_production", "_consumption")),
         x.axis = factor(x.axis, levels = ordered.phenotype))

# set the boundary (max values) of the plot
d.edge <- d.production.consumption.mirror %>% 
  filter(phenotype == "ob/ob") %>% 
  mutate(side = ifelse(nmol.min.animal.mean > 0, "up", "down")) %>% 
  group_by(side, about) %>% 
  summarise(edge = sum(abs(nmol.min.animal.mean))) %>% 
  group_by(about) %>% 
  summarise(edge = max(edge))


# # visualization 
# # Per animal 
# plt.production.consumption.mirror <-  
#   d.production.consumption.mirror %>% 
#   ggplot(aes(x = x.axis, y = nmol.min.animal.mean, fill = myfill)) +
#   geom_col(color = "black") +
#   facet_wrap(~ about, scales = "free_y", nrow = 2) +
#   geom_errorbar(aes(ymax = err.Y.animal,
#                     ymin = ifelse(nmol.min.animal.mean > 0, 
#                                   err.Y.animal - nmol.min.animal.SEM,
#                                   err.Y.animal + nmol.min.animal.SEM)),
#                 width = .4) +
#   geom_hline(yintercept = 0, linewidth = 1) +
#   scale_fill_manual(
#     values = color.Compounds, 
#     name = "source of production, \nor destiny of consumption",
#     labels = c("non-Ox sink" = "Storage (recycle)",
#                "storage" = "Storage (release)")) +
#   theme_minimal(base_size = 14) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         strip.background = element_blank(),
#         strip.text = element_text(face = "bold", size = 14),
#         axis.text = element_text(color = "black", size = 13),
#         axis.ticks = element_line(color = "snow4"),
#         panel.spacing = unit(20, "pt"),
#         panel.grid = element_blank(),
#         # panel.grid.minor.y = element_line(linetype = "dashed", linewidth = .1),
#         # panel.grid.major.y = element_line(linetype = "dotted", linewidth = .1),
#         legend.title = element_text(size = 10))  +
#   scale_x_discrete(
#     expand = expansion(add = 1), 
#     labels = function(x){str_replace(x, "WT", "control")},
#     name = NULL) +
#   scale_y_continuous(labels = function(x){abs(x)/1000},
#                      name = "µmol C / min / animal\n",
#                      n.breaks = 7,
#                      expand = expansion(mult = c(0.05, 0.05))) +
#   
#   # as the consumption and production calculated values are not perfectly the same
#   # panels thus are not perfectly aligned in the center. To fix this:
#   # add point to mark the boundary of the plots
#   # so that the y = 0 is all aligned up across panels
#   geom_point(data = d.edge, 
#              aes(x = factor("ob/ob"), y = edge * 1.02), 
#              inherit.aes = F, size = 0, alpha = 0) +
#   geom_point(data = d.edge, 
#              aes(x = factor("ob/ob"), y = -edge * 1.02), 
#              inherit.aes = F, size = 0, alpha = 0) 
# 
# plt.production.consumption.mirror



# normalized based on per body weight
d.BW2 <- d.BW %>% group_by(phenotype) %>% summarise(BW = mean(BW.mean))
# composition of fat percent of body weight
d.composition <- tibble(
  phenotype = c("WT", "HFD", "ob/ob", "db/db"),
  fat.frac = c(0.17, 0.47, 0.51, .51))
# lean and fat mass
d.BW2 <- d.BW2 %>% left_join(d.composition) %>% 
  mutate(fat = BW * fat.frac, lean = BW * (1-fat.frac)) %>% 
  select(-fat.frac)

d.BW2


# calculate fluxes based on per g BW, per g fat mass, or per g of lean mass
d.production.consumption.mirror.normalized <- d.production.consumption.mirror %>%
  left_join(d.BW2) %>%
  mutate(
    # BW
    nmol.min.gBW = nmol.min.animal.mean / BW,
    nmol.min.gBW.SEM = nmol.min.animal.SEM / BW,
    err.Y.gBW = err.Y.animal / BW,
    
    # fat
    nmol.min.gFat = nmol.min.animal.mean / fat,
    nmol.min.gFat.SEM = nmol.min.animal.SEM / fat,
    err.Y.gFat = err.Y.animal / fat,
    
    # lean
    nmol.min.gLean = nmol.min.animal.mean / lean,
    nmol.min.gLean.SEM = nmol.min.animal.SEM / lean,
    err.Y.gLean = err.Y.animal / lean
  )


func.mirror <- function(
    y = "nmol.min.gFat",  # y axis
    err.Y = "err.Y.gFat", # cumulated y-axis position for error bars
    errorBar = "nmol.min.gFat.SEM", # the SEM for each conversion flux 
    norm.basis = "g Fat", # the basis of normalization
    myEdge = 35, # trial-error parameter to align all plots along the center
    fig.width = 11, fig.height = 7 # dimension for per animal basis; g basis has longer y-axis texts, needed wider width
){
  p <- d.production.consumption.mirror.normalized %>% 
    ggplot(aes(x = x.axis, y = .data[[y]], fill = myfill)) +
    geom_col(color = "black") +
    geom_errorbar(aes(ymax = .data[[err.Y]],
                      ymin = ifelse(.data[[y]] > 0, 
                                    .data[[err.Y]] - .data[[errorBar]],
                                    .data[[err.Y]] + .data[[errorBar]])),
                  width = .4) +
    geom_hline(yintercept = 0, linewidth = 1) +
    facet_wrap(~ about, scales = "free_y", nrow = 2) +
    scale_fill_manual(
      values = color.Compounds, 
      name = "source of production, \nor destiny of consumption",
      labels = c("non-Ox sink" = "Storage (recycle)",
                 "storage" = "Storage (release)")) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 14),
          axis.text = element_text(color = "black", size = 13),
          axis.ticks = element_line(color = "snow4"),
          panel.spacing = unit(20, "pt"),
          panel.grid = element_blank(),
          # panel.grid.minor.y = element_line(linetype = "dashed", linewidth = .1),
          # panel.grid.major.y = element_line(linetype = "dotted", linewidth = .1),
          legend.title = element_text(size = 10))  +
    scale_x_discrete(
      expand = expansion(add = 1), 
      labels = function(x){str_replace(x, "WT", "control")},
      name = NULL) +
    scale_y_continuous(name = paste("nmol C / min /", norm.basis, "\n"),
                       n.breaks = 7,
                       labels = function(x)abs(x),
                       expand = expansion(mult = c(0.05, 0.05))) +
    
    # as the consumption and production calculated values are not perfectly the same
    # panels thus are not perfectly aligned in the center. To fix this:
    # add point to mark the boundary of the plots
    # so that the y = 0 is all aligned up across panels
    geom_point(data = d.edge, 
               aes(x = factor("ob/ob"), y = edge * myEdge/1000), 
               inherit.aes = F, size = 0, alpha = 0) +
    geom_point(data = d.edge, 
               aes(x = factor("ob/ob"), y = -edge * myEdge/1000), 
               inherit.aes = F, size = 0, alpha = 0) 
  
  ggsave(filename = paste0("mirror ", norm.basis, ".pdf"), 
         path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
         height = fig.height, width = fig.width)
  
  return(p)
  
}

# per animal basis
func.mirror(y = "nmol.min.animal.mean",  # y axis
            err.Y = "err.Y.animal", # accumulated y-axis position for error bars
            errorBar = "nmol.min.animal.SEM", # the SEM for each conversion flux 
            norm.basis = "animal", # the basis of normalization
            myEdge = 1) +
  # update y axis: convert nmol to µmol
  scale_y_continuous(labels = ~ abs(.x) / 1000, name = "µmol C / min / animal")

ggsave(filename = paste0("mirror ", "animal", ".pdf"), 
       height = 7, width = 11,
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures")



# normalize by g BW -------------------------------------------------
func.mirror( y = "nmol.min.gBW",  # y axis
             err.Y = "err.Y.gBW", # accumulated y-axis position for error bars
             errorBar = "nmol.min.gBW.SEM", # the SEM for each conversion flux 
             norm.basis = "per g body weight", # the basis of normalization
             myEdge = 15, fig.width = 13)




# normalize by lean mass -------------------------------------------------
func.mirror( y = "nmol.min.gLean",  # y axis
             err.Y = "err.Y.gLean", # accumulated y-axis position for error bars
             errorBar = "nmol.min.gLean.SEM", # the SEM for each conversion flux 
             norm.basis = "per g lean mass", # the basis of normalization
             myEdge = 35, fig.width = 13)



# normalize by fat mass
func.mirror( y = "nmol.min.gFat",  # y axis
             err.Y = "err.Y.gFat", # accumulated y-axis position for error bars
             errorBar = "nmol.min.gFat.SEM", # the SEM for each conversion flux 
             norm.basis = "per g Fat", # the basis of normalization
             myEdge = 35, fig.width = 13)



# 
# # g BW
# d.production.consumption.mirror.normalized %>% 
#   ggplot(aes(x = x.axis, y = nmol.min.gBW, fill = myfill)) +
#   geom_col(color = "black") +
#   geom_errorbar(aes(ymax = err.Y.gBW,
#                     ymin = ifelse(nmol.min.gBW > 0, 
#                                   err.Y.gBW - nmol.min.gBW.SEM,
#                                   err.Y.gBW + nmol.min.gBW.SEM)),
#                 width = .4) +
#   geom_hline(yintercept = 0, linewidth = 1) +
#   facet_wrap(~ about, scales = "free_y", nrow = 2) +
#   scale_fill_manual(
#     values = color.Compounds, 
#     name = "source of production, \nor destiny of consumption",
#     labels = c("non-Ox sink" = "Storage (recycle)",
#                "storage" = "Storage (release)")) +
#   theme_minimal(base_size = 14) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         strip.background = element_blank(),
#         strip.text = element_text(face = "bold", size = 14),
#         axis.text = element_text(color = "black", size = 13),
#         axis.ticks = element_line(color = "snow4"),
#         panel.spacing = unit(20, "pt"),
#         panel.grid = element_blank(),
#         # panel.grid.minor.y = element_line(linetype = "dashed", linewidth = .1),
#         # panel.grid.major.y = element_line(linetype = "dotted", linewidth = .1),
#         legend.title = element_text(size = 10))  +
#   scale_x_discrete(
#     expand = expansion(add = 1), 
#     labels = function(x){str_replace(x, "WT", "control")},
#     name = NULL) +
#   scale_y_continuous(name = "nmol C / min / g BW\n",
#                      n.breaks = 7,
#                      expand = expansion(mult = c(0.05, 0.05))) +
#   
#   # as the consumption and production calculated values are not perfectly the same
#   # panels thus are not perfectly aligned in the center. To fix this:
#   # add point to mark the boundary of the plots
#   # so that the y = 0 is all aligned up across panels
#   geom_point(data = d.edge, 
#              aes(x = factor("ob/ob"), y = edge * 20/1000), 
#              inherit.aes = F, size = 0, alpha = 0) +
#   geom_point(data = d.edge, 
#              aes(x = factor("ob/ob"), y = -edge * 20/1000), 
#              inherit.aes = F, size = 0, alpha = 0) 
# 
# 
# 
# 
# # g lean mass
# d.production.consumption.mirror.normalized %>% 
#   ggplot(aes(x = x.axis, y = nmol.min.gLean, fill = myfill)) +
#   geom_col(color = "black") +
#   geom_errorbar(aes(ymax = err.Y.gLean,
#                     ymin = ifelse(nmol.min.gLean > 0, 
#                                   err.Y.gLean - nmol.min.gLean.SEM,
#                                   err.Y.gLean + nmol.min.gLean.SEM)),
#                 width = .4) +
#   geom_hline(yintercept = 0, linewidth = 1) +
#   facet_wrap(~ about, scales = "free_y", nrow = 2) +
#   scale_fill_manual(
#     values = color.Compounds, 
#     name = "source of production, \nor destiny of consumption",
#     labels = c("non-Ox sink" = "Storage (recycle)",
#                "storage" = "Storage (release)")) +
#   theme_minimal(base_size = 14) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         strip.background = element_blank(),
#         strip.text = element_text(face = "bold", size = 14),
#         axis.text = element_text(color = "black", size = 13),
#         axis.ticks = element_line(color = "snow4"),
#         panel.spacing = unit(20, "pt"),
#         panel.grid = element_blank(),
#         # panel.grid.minor.y = element_line(linetype = "dashed", linewidth = .1),
#         # panel.grid.major.y = element_line(linetype = "dotted", linewidth = .1),
#         legend.title = element_text(size = 10))  +
#   scale_x_discrete(
#     expand = expansion(add = 1), 
#     labels = function(x){str_replace(x, "WT", "control")},
#     name = NULL) +
#   scale_y_continuous(name = "nmol C / min / g lean mass\n",
#                      n.breaks = 7,
#                      expand = expansion(mult = c(0.05, 0.05))) +
#   
#   # as the consumption and production calculated values are not perfectly the same
#   # panels thus are not perfectly aligned in the center. To fix this:
#   # add point to mark the boundary of the plots
#   # so that the y = 0 is all aligned up across panels
#   geom_point(data = d.edge, 
#              aes(x = factor("ob/ob"), y = edge * 35/1000), 
#              inherit.aes = F, size = 0, alpha = 0) +
#   geom_point(data = d.edge, 
#              aes(x = factor("ob/ob"), y = -edge * 35/1000), 
#              inherit.aes = F, size = 0, alpha = 0) 
# 
# 
# 
# # fat mass 
# d.production.consumption.mirror.normalized %>% 
#   ggplot(aes(x = x.axis, y = nmol.min.gFat, fill = myfill)) +
#   geom_col(color = "black") +
#   geom_errorbar(aes(ymax = err.Y.gFat,
#                     ymin = ifelse(nmol.min.gLean > 0, 
#                                   err.Y.gFat - nmol.min.gFat.SEM,
#                                   err.Y.gFat + nmol.min.gFat.SEM)),
#                 width = .4) +
#   geom_hline(yintercept = 0, linewidth = 1) +
#   facet_wrap(~ about, scales = "free_y", nrow = 2) +
#   scale_fill_manual(
#     values = color.Compounds, 
#     name = "source of production, \nor destiny of consumption",
#     labels = c("non-Ox sink" = "Storage (recycle)",
#                "storage" = "Storage (release)")) +
#   theme_minimal(base_size = 14) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         strip.background = element_blank(),
#         strip.text = element_text(face = "bold", size = 14),
#         axis.text = element_text(color = "black", size = 13),
#         axis.ticks = element_line(color = "snow4"),
#         panel.spacing = unit(20, "pt"),
#         panel.grid = element_blank(),
#         # panel.grid.minor.y = element_line(linetype = "dashed", linewidth = .1),
#         # panel.grid.major.y = element_line(linetype = "dotted", linewidth = .1),
#         legend.title = element_text(size = 10))  +
#   scale_x_discrete(
#     expand = expansion(add = 1), 
#     labels = function(x){str_replace(x, "WT", "control")},
#     name = NULL) +
#   scale_y_continuous(name = "nmol C / min / g fat mass\n",
#                      n.breaks = 7,
#                      expand = expansion(mult = c(0.05, 0.05))) +
#   
#   # as the consumption and production calculated values are not perfectly the same
#   # panels thus are not perfectly aligned in the center. To fix this:
#   # add point to mark the boundary of the plots
#   # so that the y = 0 is all aligned up across panels
#   geom_point(data = d.edge, 
#              aes(x = factor("ob/ob"), y = edge * 35/1000), 
#              inherit.aes = F, size = 0, alpha = 0) +
#   geom_point(data = d.edge, 
#              aes(x = factor("ob/ob"), y = -edge * 35/1000), 
#              inherit.aes = F, size = 0, alpha = 0) 






# combine the consumption data with total production
plot_grid(plt.totalProduction.atom.3pheno + 
            theme(legend.position = "none",
                  axis.title.x = element_blank(),
                  axis.text.x = element_text(angle = 60, hjust = 1)),
          
          ggplot() + theme_void(),
          
          plt.CO2.sink.3pheno, 
          rel_widths = c(1.2, .2, 3), nrow = 1) 


ggsave(filename = "production CO2 sink stacked nutrients 3 pheno.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 4.3, width = 9)


d.production.consumption %>% 
  filter(from == "Valine" & to == "CO2")


# Export key quantification results into Excels
o <- d.production.consumption %>% ungroup() %>% 
  select(phenotype, from, to, nmol.min.animal.mean, nmol.min.animal.SEM) %>% 
  distinct()

write.xlsx(o %>% filter(phenotype == "WT") %>% as.data.frame(),
           file = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/flux summary.xlsx",
           sheetName = "WT")

write.xlsx(o %>% filter(phenotype == "HFD") %>% as.data.frame(),
           file = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/flux summary.xlsx",
           sheetName = "HFD", append = T)

write.xlsx(o %>% filter(phenotype == "ob/ob") %>% as.data.frame(),
           file = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/flux summary.xlsx",
           sheetName = "ob.ob", append = T)





#- # export key data object as project

# whole set of labeling matrix
m.L <- d.enrich.atom.summary %>%
  filter(phenotype == "WT") %>% ungroup() %>% 
  select(infused.tracer, Compound, enrich.original.mean) %>% 
  spread(Compound, enrich.original.mean)

save(m.L, file = "/Users/boyuan/Desktop/Harvard/Research/db db mice/Infusion data/Flux analysis script/labeling matrix.RData")


# infusion rate
k2 <- d.normalized.tidy %>% filter(phenotype == "WT") %>% 
  mutate(infused.tracer = factor(infused.tracer, levels = ordered.Compound)) %>% 
  arrange(infused.tracer) %>% 
  select(infused.tracer, infusion_nmol13C.atoms_perMin.gBW) %>% 
  distinct()
R <- k2$infusion_nmol13C.atoms_perMin.gBW
names(R) <- k2$infused.tracer

save(R, file = "/Users/boyuan/Desktop/Harvard/Research/db db mice/Infusion data/Flux analysis script/infusion_nmolC.min.gBW.RData")


# fox
k3 <- d.CO2.fox %>% 
  filter(phenotype == "WT") %>% 
  mutate(infused.tracer = factor(infused.tracer, levels = ordered.Compound))
fox <- k3$fox.mean
names(fox) <- k3$infused.tracer

save(fox, file = "/Users/boyuan/Desktop/Harvard/Research/db db mice/Infusion data/Flux analysis script/fox.RData")



#><-#><-#><-#><-#><-#><-#><-#><-#><-#><-#><-#><-#><-#><-#><-#><-#><-#><-#><-#><-#><-#><-



# calcualte the ultimate contribution fluxes from storages
func.Ultimate.storage.frac <- function(pheno.i, storage.i){
  # pheno.i = "WT"
  # storage.i = "glycogen"
  
  X <- d.directFlux_interConversion %>% 
    select(phenotype, targetCompound, sources, nmol.min.animal.mean, nmol.min.animal.SEM) 
  
  # the mean values
  d.mean <- X %>% select(-nmol.min.animal.SEM) %>% filter(phenotype == pheno.i) %>% 
    spread(key = sources, value = nmol.min.animal.mean) 
  # adjust the order of columns of circulating nutrients as the order of `targetCompound`
  d.mean.ordered <- d.mean %>% select(phenotype, targetCompound, storage, d.mean$targetCompound)
  
  # the SEM
  d.SEM <- X %>% select(-nmol.min.animal.mean) %>% filter(phenotype == pheno.i) %>% 
    spread(key = sources, value = nmol.min.animal.SEM) 
  # adjust the order of columns of circulating nutrients as the order of `targetCompound`
  d.SEM.ordered <- d.SEM %>% select(phenotype, targetCompound, storage, d.SEM$targetCompound)
  
  # for each vector, replace the NA with the negative sum of numbers (total production) in this vector
  func.replace_NA.withSum <- function(x){
    replace_na(x, replace = -sum(x, na.rm = T))
  }
  
  mat.mean <- apply(d.mean.ordered[, -c(1:2)], MARGIN = 1, func.replace_NA.withSum) %>% t()
  mat.SEM  <- apply(d.SEM.ordered [, -c(1:2)], MARGIN = 1, func.replace_NA.withSum) %>% t()
  
  
  # add row names
  rownames(mat.mean) <- colnames(mat.mean)[-1]
  rownames(mat.SEM) <- colnames(mat.SEM)[-1]
  
  # set up the main matrix onthe left side of the equation, without the storage column
  M.mean <- mat.mean[, -1]
  M.SEM  <- mat.SEM[, -1]
  
  # split the storage vector to the right side of the equations (also with a negative sign)
  v.mean <- mat.mean[, 1]
  v.SEM  <- mat.SEM[, 1]
  
  # apply the binary `a` to the vector on the right side of the equation
  if (storage.i == "glycogen"){
    a <- str_detect(colnames(M.mean), pattern = or("Glucose", "Lactate")) %>% as.double()
  } else if (
    storage.i == "triglyceride"){
    a <- str_detect(colnames(M.mean), pattern = or("C16:0", "C18:1", "C18:2", "Glycerol", "3-HB")) %>% as.double()
  } else if (
    storage.i == "protein") {
    a <- str_detect(colnames(M.mean), pattern = or("Alanine", "Glutamine", "Valine")) %>% as.double()
  }
  
  # update the vector with the binary a
  v.mean.binaried <- v.mean * a
  v.SEM.binaried  <- v.SEM * a
  
  # do Monte Carlo simulation
  c.frac <- vector()
  
  for (i in 1:100){
    mat.MonteCarlo <- M.mean + rnorm(nrow(M.mean) * ncol(M.mean), mean = 0, sd = 1) * M.SEM
    # vector on the right side ( with negative signs)
    v.MonteCarlo <-  - (v.mean.binaried + rnorm(length(v.mean.binaried), mean = 0, sd = 1) * v.SEM.binaried)
    
    # solve the equation
    list.frac = lsei(A = mat.MonteCarlo, B = v.MonteCarlo, # this solves Ax = b
                     G = diag(1, nrow = length(v.MonteCarlo)), # G is the identity matrix, dimension = B-vector length
                     H = rep(-.002, length(v.MonteCarlo)), # vector of zero, length = B-vector length
                     # occasionally, using 0 as lower limit renders no solution due to constraint inconsistency. 
                     # so here use a small negative value (.1% error) to ensure solvability
                     verbose = FALSE)
    v = list.frac$X
    c.frac <- cbind(c.frac, v)
  }
  
  d.ultimate.S.frac <- 
    tibble(compound = rownames(c.frac),
           Ult.frac.mean = c.frac %>% as_tibble() %>% apply(MARGIN = 1, mean),
           Ult.frac.SEM = c.frac %>% as_tibble() %>% apply(MARGIN = 1, sd)) %>% 
    mutate(phenotype = pheno.i) %>% 
    mutate(store = storage.i)
  return(d.ultimate.S.frac)
}

# calculate the ultimate fraction for each storage and each phenotype
d.ultimate.S.frac <-  do.call(what = rbind, list(
  func.Ultimate.storage.frac(pheno.i = "WT",    storage.i = "glycogen"),
  func.Ultimate.storage.frac(pheno.i = "HFD",   storage.i = "glycogen"),
  func.Ultimate.storage.frac(pheno.i = "ob/ob", storage.i = "glycogen"),  
  
  func.Ultimate.storage.frac(pheno.i = "WT",    storage.i = "triglyceride"),
  func.Ultimate.storage.frac(pheno.i = "HFD",   storage.i = "triglyceride"),
  func.Ultimate.storage.frac(pheno.i = "ob/ob", storage.i = "triglyceride"),  
  
  func.Ultimate.storage.frac(pheno.i = "WT",    storage.i = "protein"),
  func.Ultimate.storage.frac(pheno.i = "HFD",   storage.i = "protein"),
  func.Ultimate.storage.frac(pheno.i = "ob/ob", storage.i = "protein")  
)) %>% 
  mutate(store = factor(store, levels = c("glycogen", "triglyceride", "protein") %>% rev())) %>% 
  mutate(compound = factor(compound, levels = ordered.Compound),
         phenotype = factor(phenotype, levels = ordered.phenotype)) %>% 
  group_by(phenotype, compound) %>% 
  # scale down to a total sum of 1 (optional)
  mutate(Ult.frac.mean = Ult.frac.mean / (sum(Ult.frac.mean)) ) %>% 
  mutate(y.error = cumsum(Ult.frac.mean))


# plot


func.plt.Ult.f <- function(mydata, myx = "compound"){
  mydata %>% 
    ggplot(aes_string(x = myx, y = "Ult.frac.mean", fill = "store")) +
    geom_col(position = position_stack(), color = "black", alpha = .4) +
    geom_errorbar(aes(ymax = y.error,
                      ymin = y.error - Ult.frac.SEM), width = .3) +
    facet_wrap(~phenotype) +
    scale_fill_brewer(palette = "Pastel1") +
    scale_x_discrete(expand = expansion(add = .8)) +
    scale_y_continuous(expand = expansion(mult = c(0, .0)),
                       breaks = seq(0, 1, .2)) +
    labs(y = "Ultimate contribution fraction\n", x = NULL) +
    theme.myClassic +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom") +
    scale_fill_manual(values = c("steelblue", "firebrick",  "black"))
}  

d.ultimate.S.frac %>% 
  filter(phenotype == "WT" & compound %in% c("Glucose", "Lactate", "Alanine", "Glutamine")) %>% 
  func.plt.Ult.f() + facet_wrap(~"")

ggsave(filename = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures/ultimate storage contr fraction_WT.pdf",
       height = 4.6, width = 4.3)

# all phenotypes
d.ultimate.S.frac %>% 
  filter(compound %in% c("Glucose", "Lactate")) %>% 
  func.plt.Ult.f(myx = "phenotype") +
  facet_wrap(~compound)

ggsave(
  filename = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures/ultimate storage contr fraction_3pheno selected nutrient.pdf",
  height = 3.5, width = 6)


# save project
save.image(file = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/7_core_production_flux.RData")


# export key data to calculate interconversion flux as a template for labmates' use
save(d.enrich.atom.summary, 
     d.normalized.tidy,
     file = "/Users/boyuan/Desktop/Harvard/Research/db db mice/Infusion data/data_for_flux_calculate_template.RData")




