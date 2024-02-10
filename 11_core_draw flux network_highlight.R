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
setwd("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/")

# create a function to specify the bin
b <- c(1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3)

func.createBins <- function(x){
  # Find the interval in which x falls in vector b
  interval <- findInterval(x, b)
  
  # Determine the two elements that sandwich x
  elements_sandwiching_x <- c(b[interval], b[interval + 1])
  
  paste0(elements_sandwiching_x[1], "-", elements_sandwiching_x[2]) %>% return()
}

func.createBins(2.3)

# max flux of each phenotype
d.production.consumption %>% group_by(phenotype) %>% 
  summarise(max = max(nmol.min.animal.mean))

# Plot flux network
d.directFlux_interConversion.selected = d.directFlux_interConversion %>% 
  select(phenotype, targetCompound, sources, nmol.min.animal.mean) %>% 
  filter(phenotype %in% c("WT", "HFD", "ob/ob")) %>% 
  mutate(nmol.min.animal.mean = round(nmol.min.animal.mean, 2))

# compound in order in network, clockwisely
# ordered.Compound.network =  
#   c("Glucose", "Lactate", "Alanine",   
#     "C16:0", "C18:1", "C18:2", "3-HB", "Valine", "Glutamine", "Glycerol") 

ordered.Compound.network =  
  c("Lactate","Glucose", 
    "Glycerol", "C16:0", "C18:1", "C18:2", "3-HB", "Glutamine", "Valine", "Alanine") 

# function to abbreviate the label name
func.abbreviate <-  function(x){
  x <- x %>% str_replace("Glucose", "Glc")
  x <- x %>% str_replace("Lactate", "Lac")
  x <- x %>% str_replace("Alanine", "Ala")
  x <- x %>% str_replace("3-HB", "3-HB")
  x <- x %>% str_replace("Valine", "Val")
  x <- x %>% str_replace("Glutamine", "Gln")
  x <- x %>% str_replace("Glycerol", "Gly")
  
}

d.directFlux_interConversion.selected = d.directFlux_interConversion.selected %>% 
  mutate(targetCompound = factor(targetCompound, levels = ordered.Compound.network, ordered = T))

d.directFlux_interConversion.selected$targetCompound %>% unique()
ordered.Compound.network



# calculate metabolite coordinate in the network 
n.compound = ordered.Compound.network %>% length()
theta.delta = 2 * pi / n.compound
theta = 2/pi + theta.delta * c(1:n.compound)

coord.x = cos(theta)   
coord.y = sin(theta)  

# assign position to each metabolite in clockwise order 
names(coord.x) = ordered.Compound.network
names(coord.y) = ordered.Compound.network

d.directFlux_interConversion.selected.coordinated = 
  d.directFlux_interConversion.selected %>% 
  mutate(coord.x.target = coord.x[targetCompound],
         coord.y.target = coord.y[targetCompound])

# calculate fold change mapped to color
# WT fluxes as the basis to calculate fold change
d.directFlux_interConversion.selected.coordinated <- 
  d.directFlux_interConversion.selected.coordinated %>%
  filter(phenotype == "WT") %>%
  mutate(WT.flx = nmol.min.animal.mean + 1) %>%
  ungroup() %>%
  select(targetCompound, sources, WT.flx) %>%
  # combine with original dataset to calculate fold change
  left_join(d.directFlux_interConversion.selected.coordinated, ., by = c("targetCompound", "sources")) %>%
  mutate(fold.change = nmol.min.animal.mean / WT.flx) %>% select(-WT.flx) %>% 
  mutate(fold.change = ifelse(fold.change > 3, 3, fold.change)) #  %>%
  # mutate(fold.change = ifelse(fold.change <= 1, 1, fold.change))


# Cleanup the flux significance data 
func.toStorage <- function(x){
  x %>% str_replace("Protein", "storage") %>% 
    str_replace("TAG", "storage") %>% 
    str_replace("Glycogen", "storage") %>% 
    str_replace("non-Ox sink", "storage") %>% 
    str_replace("Storage", "storage")
}

d.p2 <- d.p %>% mutate(to = func.toStorage(to), from = func.toStorage(from)) %>% 
  rename(targetCompound = to, sources = from) %>% 
  mutate(phenotype = str_extract(compare, pattern = or("HFD", "ob")) %>% str_replace("ob", "ob/ob")) %>% 
  # if adjusted p value is significant
  mutate(p.adj.significant = (p.adj < 0.05) %>% as.character()) %>% 
  select(-compare)


# combine with fold change
d.directFlux_interConversion.selected.coordinated <- 
  d.directFlux_interConversion.selected.coordinated %>% 
  left_join(d.p2, by = c("targetCompound", "sources", "phenotype")) %>% 
  # for wt fold change, set pvalue and adjusted p value to be 1
  # i.e, replace NA values with 1
  mutate(pvalue = replace_na(pvalue, 1),
         p.adj = replace_na(p.adj, 1)) %>% 
  # if adjusted p value is not significant, set fold change to 1
  mutate(fold.change = ifelse(p.adj <= 0.05, fold.change, 1)) %>% 
  mutate(fold.change.bin = cut_width(fold.change, width = .2))


# Anchor points of target compound circle
d.anchor_Metabolite = d.directFlux_interConversion.selected.coordinated %>% 
  select(targetCompound, coord.x.target, coord.y.target) %>% distinct()



# flux to CO2 and to non-Ox sink
d.destiny.CO2.sink = d.destiny %>% 
  select(phenotype, Compounds, destiny, nmol.min.animal.mean) %>% 
  mutate(nmol.min.animal.mean = round(nmol.min.animal.mean, 2)) %>% 
  filter(destiny %in% c("non-Ox sink", "CO2")) %>% 
  rename(targetCompound  = Compounds) %>% mutate(targetCompound = as.character(targetCompound)) %>% 
  
  filter(phenotype %in% c("WT", "ob/ob", "HFD")) %>% 
  
  left_join(d.anchor_Metabolite %>% mutate(targetCompound = as.character(targetCompound))) %>% # add nutrient coord position
  mutate(coord.x.CO2 = 0, coord.y.CO2 = 0)

# add fold change
d.destiny.CO2.sink <- d.destiny.CO2.sink %>% ungroup() %>% 
  # WT flux as basis to calculate fold change
  filter(phenotype == "WT") %>% 
  mutate(WT.flx = nmol.min.animal.mean + 1) %>% 
  select(targetCompound, destiny, WT.flx) %>% 
  # merge with original dataset
  left_join(d.destiny.CO2.sink, ., by = c("targetCompound", "destiny")) %>% 
  # calculate fold change
  mutate(fold.change = nmol.min.animal.mean / WT.flx) %>% 
  mutate(fold.change = ifelse(fold.change > 3, 3, fold.change)) %>%
  # mutate(fold.change = ifelse(fold.change <= 1, 1, fold.change)) %>% 
  select(-WT.flx) %>% 
  # update storage name before join
  mutate(destiny = func.toStorage(destiny)) %>% 
  # add significance
  left_join(
    # update column name
    d.p2 %>% rename(from = sources,to = targetCompound) %>% rename(destiny = to, targetCompound = from), 
    by = c("targetCompound", "destiny", "phenotype")) %>% 
  # for wt fold change, set pvalue and adjusted p value to be 1
  # i.e, replace NA values with 1
  mutate(pvalue = replace_na(pvalue, 1),
         p.adj = replace_na(p.adj, 1)) %>% 
  # if adjusted p value is not significant, set fold change to 1
  mutate(fold.change = ifelse(p.adj <= 0.05, fold.change, 1)) %>% 
  mutate(fold.change.bin = cut_width(fold.change, width = .2))

d.destiny.CO2.sink



# Inter-conversion flux among circulating nutrients
d.layer.interConvert = 
  d.directFlux_interConversion.selected.coordinated %>%
  mutate(sources = as.character(sources)) %>% 
  filter(sources != "storage") %>% 
  
  left_join(d.anchor_Metabolite %>% rename(sources = targetCompound, 
                                           coord.x.source = coord.x.target, 
                                           coord.y.source = coord.y.target))  # add coord of source 


# Define edge size 
nmol.min.animal = c(d.destiny.CO2.sink$nmol.min.animal.mean, # flux to CO2 and non-Ox sink
                    d.directFlux_interConversion.selected$nmol.min.animal.mean) %>% unique() # interconversion flux and flux from storage to nutrients

cutoff.flx = 5000 # nmol/min/animal
fold.maxFlx.cutOff = (max(nmol.min.animal) / cutoff.flx) %>% ceiling()

# ----- Big flx: 20,000 ~ 4,000  nmol/min/animal
nmol.min.animal.BigFlx = nmol.min.animal[nmol.min.animal >= cutoff.flx] %>% sort()
size.range_Big = seq(10/fold.maxFlx.cutOff, 12, length.out = length(nmol.min.animal.BigFlx)) %>% round(2)

# plot(size.range_Big, nmol.min.animal.BigFlx)

mdl.size_big = lm(size.range_Big ~  nmol.min.animal.BigFlx )
nmol.min.animal.BigFlx.fitted = predict(mdl.size_big)

names(nmol.min.animal.BigFlx.fitted) = nmol.min.animal.BigFlx
nmol.min.animal.BigFlx.fitted


# ----- Mid flx: 4,000 ~ 800 nmol/min/animal (five fold lower)
nmol.min.animal.mid_Flx = nmol.min.animal[nmol.min.animal < cutoff.flx & nmol.min.animal >= cutoff.flx/5 ] %>% sort()
size.range_Mid = seq(1, nmol.min.animal.BigFlx.fitted %>% min() - .5, # max of mid-flx size using the min of fitted value of large flux size
                     length.out = length(nmol.min.animal.mid_Flx)) %>% round(2) 

# plot(size.range_Mid, nmol.min.animal.mid_Flx)

mdl.size_Mid = lm(size.range_Mid ~  nmol.min.animal.mid_Flx )
nmol.min.animal.Mid.Flx.fitted = predict(mdl.size_Mid)

names(nmol.min.animal.Mid.Flx.fitted) = nmol.min.animal.mid_Flx
nmol.min.animal.Mid.Flx.fitted



# ----- small flx: 800 ~ 160 nmol/min/animal (another five fold lower)
nmol.min.animal.small_Flx = nmol.min.animal[nmol.min.animal < cutoff.flx/5 & nmol.min.animal >= cutoff.flx/5/5] %>% sort()
size.range_small = seq(0, min(nmol.min.animal.Mid.Flx.fitted) - .5,
                       length.out = length(nmol.min.animal.small_Flx)) %>% round(2) 

# plot(size.range_small, nmol.min.animal.small_Flx)

mdl.size_small = lm(size.range_small ~  nmol.min.animal.small_Flx )
nmol.min.animal.Small_Flx.fitted = predict(mdl.size_small)

names(nmol.min.animal.Small_Flx.fitted) = nmol.min.animal.small_Flx
nmol.min.animal.Small_Flx.fitted

# combine
nmol.min.animal.Flxes.fitted = c(nmol.min.animal.Small_Flx.fitted,
                                 nmol.min.animal.Mid.Flx.fitted,
                                 nmol.min.animal.BigFlx.fitted)


# Check flux vs line size
curveSizes = nmol.min.animal.Flxes.fitted
fluxes = names(nmol.min.animal.Flxes.fitted)
plot(fluxes, curveSizes)


# fit a model to predict the size based on flux - used for making legend manually
mdl <- smooth.spline(x = fluxes, y = curveSizes)
# line size for legend
legendFlux = c(1, 5, 10, 15, 20, 30, 40) * 1000
newsize <-  predict(mdl, x = legendFlux)$y
names(newsize) = legendFlux

# update the size scale
nmol.min.animal.Flxes.fitted = c(nmol.min.animal.Flxes.fitted, newsize)


# Size testing
# k = tibble(flx = names(nmol.min.animal.Flxes.fitted) %>% as.double() , 
#            mysize = nmol.min.animal.Flxes.fitted %>% as.double() ) %>% as_tibble()
# k %>% 
#   ggplot(aes(x = flx, y = mysize, size = as.character(flx))) + 
#   geom_point(show.legend = F, position = position_jitter(1, 1))  +
#   scale_size_manual(values = nmol.min.animal.Flxes.fitted)

# PLOT

curveColor.store = "snow3"
curveColor.interconvert = "snow3"
color.CO2 = "snow3"
curveAlpha = .7
curvature.S.S = .7


fold.changes.bins <- c(
  d.directFlux_interConversion.selected.coordinated$fold.change.bin,
  d.destiny.CO2.sink$fold.change.bin) %>% unique() %>% sort() %>% as.character()

func.plt.network = function(whichPhenotype = "WT"){
  
  # whichPhenotype = "ob/ob"
  
  d.directFlux_interConversion.selected.coordinated.i = 
    d.directFlux_interConversion.selected.coordinated %>% 
    filter(phenotype == whichPhenotype)
  
  d.destiny.CO2.sink.i = d.destiny.CO2.sink %>% 
    filter(phenotype == whichPhenotype)
  
  d.layer.interConvert.i = d.layer.interConvert %>% 
    filter(phenotype == whichPhenotype)
  
  d.anchor_Metabolite.i = d.anchor_Metabolite %>% 
    filter(phenotype == whichPhenotype)
  
  # layer 1
  # flux from storage to metabolites, defining the boundary of the plot
  p1 = ggplot(mapping = aes(color =  fold.change, # fold.change.bin # , alpha = p.adj.significant
                            )) + 
    geom_curve(data = 
                 d.directFlux_interConversion.selected.coordinated.i %>%
                 
                 mutate(sources = as.character(sources)) %>% 
                 filter(sources == "storage"),
               aes(x = coord.x.target * 1.8, 
                   y = coord.y.target * 1.8, 
                   xend = coord.x.target , 
                   yend = coord.y.target,
                   size = nmol.min.animal.mean %>% as.character()),
               curvature = curvature.S.S, alpha = curveAlpha) 
  
  
  p2 = p1 +
    # add flux back to non-Ox sink
    geom_curve(data = d.destiny.CO2.sink.i %>% filter(destiny == "storage"),
               aes( x = coord.x.target, 
                    y = coord.y.target, 
                    xend = coord.x.target * 1.8,
                    yend = coord.y.target * 1.8,  
                    size = nmol.min.animal.mean %>% as.character() ),
               curvature = curvature.S.S, alpha = curveAlpha) +
    
    # add flux to CO2
    geom_curve(data = d.destiny.CO2.sink.i %>% filter(destiny == "CO2"),
               aes(x = coord.x.target, y = coord.y.target, 
                   xend = 0, yend = 0, 
                   size = nmol.min.animal.mean %>% as.character() ),
               curvature = 0.4, alpha = curveAlpha) 
  
  
  p3 = p2 +  
    # large flux 
    geom_curve(data = d.layer.interConvert.i %>% filter(nmol.min.animal.mean >  cutoff.flx ),
               aes(x = coord.x.source,    y = coord.y.source,
                   xend = coord.x.target, yend = coord.y.target,
                   size = nmol.min.animal.mean %>% as.character() ),
               curvature = .25, alpha = curveAlpha) +
    
    # Mid flux
    geom_curve(data = d.layer.interConvert.i %>% filter( nmol.min.animal.mean <= cutoff.flx & nmol.min.animal.mean > cutoff.flx / 4),
               aes(x = coord.x.source,    y = coord.y.source,
                   xend = coord.x.target, yend = coord.y.target,
                   size = nmol.min.animal.mean %>% as.character()),
               curvature = .25, alpha = curveAlpha) +
    # # small flux
    geom_curve(data = d.layer.interConvert.i %>% filter( nmol.min.animal.mean <= cutoff.flx/4 & nmol.min.animal.mean > cutoff.flx/4/4),
               aes(x = coord.x.source,    y = coord.y.source,
                   xend = coord.x.target, yend = coord.y.target,
                   size = nmol.min.animal.mean %>% as.character()),
               curvature = .25, alpha = curveAlpha) +
    
    # minimal flux with arbitrary line size starting from here
    geom_curve(data = d.layer.interConvert.i %>% filter( nmol.min.animal.mean <= cutoff.flx/4/4 & nmol.min.animal.mean > 100),
               aes(x = coord.x.source,    y = coord.y.source,
                   xend = coord.x.target, yend = coord.y.target),
               size = .15, alpha = curveAlpha,
               curvature = .2) +
    # minimal trace flux
    geom_curve(data = d.layer.interConvert.i %>% filter( nmol.min.animal.mean <= 100 & nmol.min.animal.mean > 10),
               aes(x = coord.x.source,    y = coord.y.source,
                   xend = coord.x.target, yend = coord.y.target),
               size = .15, alpha = curveAlpha,
               curvature = .2) 
  
  # layer 4: nutrient names
  p4 = p3 + 
    # storage / non-Ox sink
    geom_label( data = d.anchor_Metabolite.i, 
                aes(label = "S", x = coord.x.target* 1.8, y = coord.y.target*1.8), 
                size = 7, fontface = "bold", fill = "snow2", 
                color = "black", # replace the color aesthetic mapping
                label.padding = unit(2, "mm"), 
                alpha = 1) + # overwrite the aesthetic mapping inheritance
    
    # nutrients
    geom_label( data = d.anchor_Metabolite.i, 
                aes(label = targetCompound %>% func.abbreviate(),
                    x = coord.x.target, y = coord.y.target),
                size = 6.4, fontface = "bold", fill = "turquoise",
                color = "black", # replace the color aesthetic mapping
                label.padding = unit(2, "mm"), 
                # rounded corner
                label.r = unit(10, "pt"),
                label.size = NA, alpha = 1) +
    
    # CO2 background
    annotate(geom = "point", x = 0, y = 0, stroke = 1.3,
             alpha = .9, fill = "white", color = "snow4", size = 25, shape = 21)  +
    
    # CO2 logo
    annotate("text", x = -.04, y = .01, label = 'bold("CO")',
             colour = "black", size= 8, parse = T) +
    
    annotate("text", x = .15, y = -.06 , label = 'bold("2")',
             colour = "black",  size = 6, parse = T) 
  
  
  p5 = p4 + 
    scale_size_manual(values = nmol.min.animal.Flxes.fitted)  +
    coord_fixed() +
    theme_void() +
    theme(legend.title = element_blank(),
          # legend.position = "NA"
    ) +
    guides(size = "none")
  p5 
}


# function to create fux linewidth legend
func.add_Size_Legend <- function(p){
  
  legend.flx = legendFlux[c(1:5, 7)] %>% rev()
  legend.y = seq(-1.4, .7, length.out = length(legend.flx)) 
  legend.y <- c(-1.50,  -0.98,  -0.56,  -0.14,   0.28,  0.70)
  
  p + # add legend manually
    geom_segment(data = tibble(x = 2.4, 
                               xend = 2.9,
                               y = legend.y,
                               yend = legend.y,
                               flx = legend.flx),
                 aes(x = x, xend = xend, y = y, yend = yend,
                     size = flx %>% as.character()),
                 color = "snow4", alpha = 1) +
    # note flux
    annotate(geom = "text",
             x = 2.65, 
             # y = legend.y + seq(.25, .15, length.out = 6),
             y = legend.y + c(.25, .2, .15, .12, .11, .1),
             label = legend.flx/1000, hjust = .5) # convert nmol to µmol
}

# function to create color legend 
func.add_Color_alpha_Legend <- function(p){
  
  mycolors <- c("steelblue3", "steelblue1",  "snow3", "gold1", "orange", "red", "red3")
  
  p + scale_color_gradientn(
      colors = mycolors, 
      values = c(0, seq(.2, 1, length.out = 5)), 
      limits = c(0, 3),
      guide = guide_colorbar(barheight = unit(300, "pt"))) 
    # scale_color_manual(values = mycolors, breaks = fold.changes.bins,
    #                    guide = "none") +
    # annotate(geom = "tile", 
    #          x = seq(-2.3, 1, length.out = length(mycolors)) + 0.1 * 1:length(mycolors),
    #          y = -2.2, color = "white", # linewidth = 3,
    #          width = .5, height = .13,
    #          fill = mycolors) +
    # annotate(geom = "text", 
    #          x = seq(-2.3-.25, 1-.25, length.out = length(mycolors)) + 0.1 * 1:length(mycolors),
    #          y = -2.4, label = names(mycolors) %>% str_extract("\\d.\\d")) # +
    
    # alpha scale
    # adjusted p value significant in full color; otherwise in half transparency
    # scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = .2))
}

# WT
(func.plt.network(whichPhenotype = "WT") %>% func.add_Size_Legend()) %>% 
  func.add_Color_alpha_Legend()

# HFD
plt.network.HFD <- (func.plt.network(whichPhenotype = "HFD") %>% func.add_Size_Legend()) %>% 
  func.add_Color_alpha_Legend()
plt.network.HFD
ggsave(filename = "network HFD_color_scale.pdf", path = "R Figures", height = 6, width = 8.7)  

# ob/ob
plt.network.ob <- (func.plt.network(whichPhenotype = "ob/ob") %>% func.add_Size_Legend()) %>% 
  func.add_Color_alpha_Legend() 
plt.network.ob
ggsave(filename = "network ob_ob_color_scale.pdf", path = "R Figures", height = 6, width = 8.7)  


