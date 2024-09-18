rm(list = ls())

library(ComplexHeatmap)
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


# load workspace
load("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/4_core_serum_labeling_import data.RData")


# Calculate the adjusted labeling enrichment, standardized to the same infusion rate & same tracer concentration across animals and batches of experiments
# separately for each tracer, and each phenotype / genotype
# e.g., glucose concentration infused was around 200 mM but not exactly the same in different experiments; 
# here we adjust the labeling to 200 mM concentration for 13C-glucose 

d.normalized.tidy = d.normalized.tidy %>%
  mutate(infusion_mouseID = str_c(infusion_round, "_", mouse_ID)) %>%
  
  # total carbon number of each metabolite
  group_by(Compound) %>% 
  mutate(C_Label.max = max(C_Label)) %>% ungroup() %>% 
  mutate(infusion_uL_perMin_g.BW = infusion_uL_perMin / BW) %>% 
  
  # combine with the standard infusion parameter database
  left_join(d.standardInfusionParameters, by = c("infused.tracer", "phenotype")) %>% 
  
  # adjust the enrichment in proportion to the infusion rate & tracer concentration normalization/standardization
  mutate(enrichment = ifelse (C_Label > 0, 
                              (tracer_conc_mM_Normalize / tracer_conc_mM) *( infusion_uL_perMin_g.BW_Normalize/infusion_uL_perMin_g.BW) * enrichment,
                              enrichment)) %>% 
  group_by(infusion_round, sample, Compound) 

# calculate the sum labeling of isotopes M+1, M+2...of each metabolite in each sample,
# to be subtracted from 1 as the updated parent labeling
d.normalized.tidy = d.normalized.tidy %>% filter(C_Label > 0) %>% 
  summarise(enrichment.isotope.sum = sum(enrichment, na.rm = T)) %>% 
  right_join(d.normalized.tidy) %>% 
  mutate(enrichment = ifelse(C_Label == 0,  
                             1 - enrichment.isotope.sum, enrichment) ) %>% 
  
  # use the "standardized" concentration and infusion rate: replace the old two columns with new columns using same column name
  select(-c(tracer_conc_mM, infusion_uL_perMin_g.BW)) %>% 
  rename(tracer_conc_mM = tracer_conc_mM_Normalize, 
         infusion_uL_perMin_g.BW = infusion_uL_perMin_g.BW_Normalize) %>% # updated standardized tracer concentration and infusion rate
  
  # update into standardized tracer nmol 13C-atoms infused /min/animal
  mutate(infusion_nmol13C.atoms_perMin.gBW = tracer_conc_mM * infusion_uL_perMin_g.BW * C_Label.max_tracer,
         infusion_nmol13C.atoms_perMin = infusion_nmol13C.atoms_perMin.gBW * BW, 
         infusion_uL_perMin = infusion_uL_perMin_g.BW * BW) # update whole body infusion rate as well!



# Plot infusion parameters
# d.normalized.tidy %>%
#   ggplot(aes(x = infused.tracer, y = infusion_uL_perMin_g.BW, color = phenotype)) +
#   geom_point(position = "jitter", size = .1) + expand_limits(y = 0)
# 
# d.normalized.tidy %>%
#   ggplot(aes(x = infused.tracer, y = infusion_uL_perMin, color = phenotype)) +
#   geom_point(position = "jitter", size = .1) + expand_limits(y = 0)
# 
# d.normalized.tidy %>%
#   # filter(Compound == infused.tracer) %>%
#   ggplot(aes(x = infused.tracer, y = infusion_nmol13C.atoms_perMin, color = Compound)) +
#   geom_point(position = "jitter", size = .5, shape = 21) + expand_limits(y = 0) +
#   facet_wrap(~phenotype)




# color setup for C-label: same color assignment rule for all Compounds
# Specify the first seven colors (M+0, ...M + 6)
colors = c ("grey", "firebrick", "yellow", "turquoise2",  "skyblue2", "steelblue4", "tomato") 

# Labeling equal or higher than M+7 take colors interpolated from palette Dark2, based on max labeling from the input dataset
colors.more = colorRampPalette(brewer.pal(8, "Pastel1"))( (nmax <- d.corrected.tidy$C_Label %>% as.numeric() %>% max()) - length(colors) + 1)
colors = c(colors, colors.more)
names(colors) = 0:nmax %>% factor(ordered = F)

scales::show_col(colors)

# Plot labeling enrichment of metabolites in circulation ====
flx.plot_labeling_enrichment = function(listed.tidyData = listed.tidyData,
                                        mylabeled.compound.layer1 = NULL, # labeled metabolite
                                        enrichment_lower_bound = 0,
                                        facet.nrow = 1,
                                        show.TIC = T) {
  
  # Check there is listed data input
  if (is.null(listed.tidyData)) {
    stop("Please specify the dataset: the output from function \"flx.correct_natural_abundance\", which is in the format of a list.\n\n")
  }
  
  # Extract dataset from the list
  d.corrected.tidy = listed.tidyData$Corrected
  d.normalized.tidy = listed.tidyData$Normalized
  
  
  cmpd.all = d.corrected.tidy$Compound %>% unique()
  
  # Check there is Compound input
  if (is.null(mylabeled.compound.layer1 )) {
    stop("A Compound name must be specified in the argument of \"mylabeled.compound.layer1 = ... \" to make a plot.\n ")
  }
  
  # Check Compound is of length of one
  if (length(mylabeled.compound.layer1) > 1) {
    stop("Please input only one Compound in the argument of \"mylabeled.compound.layer1 = ... \" when making a plot.\n\n")
  }
  # Check Compound matches Compound names in the input listed.data
  if (! mylabeled.compound.layer1 %in% cmpd.all) {
    stop("Compound", " \"", mylabeled.compound.layer1, "\"", " is not found in the input dataset.\n" )
  }
  
  # Convert C-label into factor to visualize in order of labeling 
  ordered.C_label = d.normalized.tidy $ C_Label %>% unique() %>% rev()
  d.normalized.tidy$C_Label = factor(d.normalized.tidy $ C_Label, levels = ordered.C_label, ordered = F)
  d.corrected.tidy$C_Label = factor(d.corrected.tidy$C_Label, levels = ordered.C_label, ordered = F)
  
  
  # Select specified Compound and tidy up
  d.corrected.tidy = d.corrected.tidy %>% filter(Compound == mylabeled.compound.layer1) 
  
  d.TIC.tidy = d.corrected.tidy %>% group_by(sample, Compound) %>%
    summarise(TIC = sum(intensity, na.rm = T),
              phenotype = unique(phenotype))
  
  d.normalized.tidy = d.normalized.tidy %>%  filter(Compound == mylabeled.compound.layer1) 
  
  
  # max TIC intensity.
  TIC.max = (d.TIC.tidy)$TIC %>% max(na.rm = T)
  # Number of color
  C_label.i = d.normalized.tidy$C_Label %>% unique() %>% as.character()
  
  
  # Define plotting function
  theme_set(theme_bw() +
              theme(#axis.text.x = element_text(angle = 45, hjust = 1),
                strip.background = element_blank(),
                strip.text = element_text(size = 14, face = "bold"),
                panel.grid = element_blank(),
                panel.border = element_rect(colour = "black", size = 1),
                axis.text = element_text(colour = "black", size = 11),
                title = element_text(face = "bold", hjust = .5, size = 12),
                plot.title = element_text(hjust = .5, size = 16),
                legend.title = element_blank()) )
  
  
  # Label fraction bar plot
  plt.bar = d.normalized.tidy %>%
    ggplot(aes(x = sample)) +
    geom_bar(aes(y = enrichment, fill = C_Label, color = C_Label),
             stat = "identity", alpha = .8) +
    
    scale_color_manual(values = colors [C_label.i %>% as.character()] ) +
    scale_fill_manual(values = colors [C_label.i %>% as.character()] ) +
    labs(x = " ", y = "Labelling fraction\n", title = mylabeled.compound.layer1)  + # \n increase y axis - text gap
    scale_x_discrete(expand = c(0, 0)) +
    coord_cartesian(ylim = c(enrichment_lower_bound, 1)) +
    scale_y_continuous(expand = c(0, 0),
                       # Need this right side axis as place holder
                       sec.axis = sec_axis(~.*TIC.max, name = "Total isotopologues ion counts\n")) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, color = ifelse(show.TIC == F, "black", "white") ),
          axis.text.y.right = element_text(color = NA),  # turn off TIC axis text
          axis.title.y.right = element_text(colour = NA),
          axis.ticks.y.right = element_blank() ) # turn off TIC ticks
  
  # TIC line plot
  plt.TIC = d.TIC.tidy %>%
    ggplot(aes(x = sample)) +
    
    # set bar plot (transparent) so as to keep the legend as a place holder to allow overlay with real bar plot
    geom_bar(aes(y = TIC),
             stat = "identity", alpha = 0, color = NA, position = "fill") +
    
    labs(x = " ", y = "Labelling fraction\n", title = mylabeled.compound.layer1)  + # \n increase y axis - text gap
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0),
                       sec.axis = sec_axis(~.*TIC.max, name = "Total isotopologues ion counts\n")) +
    
    theme(panel.background = element_blank(),
          
          # turn on bar plot as place holder for legend, but not show legend here
          # so as to show legend of the bar plot
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.text = element_blank(),
          
          plot.background = element_blank(),
          axis.ticks.y.left = element_blank(), # turn off label fraction axis ticks
          axis.text.x = element_text(angle = 60, hjust = 1, colour = "black"),
          axis.ticks.x = element_blank(),
          axis.text.y.left = element_text(colour = NA), # turn off label fraction axis text
          axis.title.y.left = element_text(colour = NA) # turn off label fraction axis title
    )  +
    
    # TIC counts
    geom_point(aes(y = TIC  / TIC.max), fill = "black", color = "white", alpha = 1, size = 3, shape = 23, stroke = 1 ) 
  # geom_line(aes(y = TIC  / TIC.max, group = 1), color = "black", alpha = .5)
  
  # facet
  plt.TIC = plt.TIC + facet_wrap( ~ phenotype, scales = "free_x", nrow = facet.nrow) + 
    theme(panel.spacing = unit(1, "lines"))
  
  plt.bar = plt.bar + facet_wrap( ~ phenotype, scales = "free_x", nrow = facet.nrow) + 
    theme(panel.spacing = unit(1, "lines"))
  
  
  # Overlay and align up the two plots
  if (show.TIC == T) {
    aligned_plots <- align_plots(plt.bar, plt.TIC, align="hv", axis="tblr")
    ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
  } else { return(plt.bar)}
}

# Wrap the above function to visualize labeling of selected tracer, metabolite, and blood
flx.plot_labeling_enrichment.which.Tracer.Blood = 
  function(my.infused.tracer = "Glucose", 
           mylabeled.compound.layer2 = "Glucose",
           plotBlood = "tail", 
           enrichment_lower_bound = 0,
           showTIC = T){
    
    listed.tidyData = list(Normalized = d.normalized.tidy %>% filter(infused.tracer == my.infused.tracer & blood == plotBlood), 
                           Corrected = d.corrected.tidy %>% filter(infused.tracer == my.infused.tracer & blood == plotBlood))
    
    p = flx.plot_labeling_enrichment(listed.tidyData = listed.tidyData, 
                                     mylabeled.compound.layer1 = mylabeled.compound.layer2,
                                     enrichment_lower_bound = enrichment_lower_bound, show.TIC = showTIC)
    
    plot_grid(
      ggplot() + theme_void() + 
        ggtitle(paste("Infused tracer: ", my.infused.tracer, "; ", plotBlood, " blood")) +
        theme(plot.title = element_text(hjust = .4, face = "bold", size = 15, 
                                        color = ifelse(plotBlood == "art", "firebrick", "steelblue"))),
      p, 
      rel_heights = c(1, 15), nrow = 2 
    ) 
  }


# visualize: labeling of selected metabolite in venous or arterial blood during infusion of selected tracer
flx.plot_labeling_enrichment.which.Tracer.Blood(
  my.infused.tracer = "Glucose", mylabeled.compound = "Glucose",
  plotBlood = "tail", enrichment_lower_bound = .5
)


flx.plot_labeling_enrichment.which.Tracer.Blood(
  my.infused.tracer = "C18:2", mylabeled.compound = "3-HB",
  plotBlood = "tail", enrichment_lower_bound = .85
)


flx.plot_labeling_enrichment.which.Tracer.Blood(
  my.infused.tracer = "3-HB", mylabeled.compound = "3-HB",
  plotBlood = "tail", enrichment_lower_bound = .6
)





# summarize and export averaged isotopomer distribution
d.isotopomer <-  d.normalized.tidy %>% 
  group_by(phenotype, infused.tracer, Compound, C_Label) %>% 
  summarise(enrich.pct.mean = mean(enrichment*100) %>% round(1),
            enrich.pct.sd = sd(enrichment*100) %>% round(1)) %>% 
  # if the mean is zero, do not show SD
  mutate(enrich = ifelse(enrich.pct.mean == 0, 0, str_c(enrich.pct.mean, " ± ", enrich.pct.sd)),
         .keep = "unused") %>% 
  mutate(infused.tracer = factor(infused.tracer, levels = ordered.Compound),
         Compound = factor(Compound, levels = ordered.Compound)) %>% 
  arrange(infused.tracer, Compound)

# convert to wider format
d.isotopomer.spread <- d.isotopomer %>% 
  pivot_wider(names_from = C_Label, values_from = enrich)

# export as excel
d.isotopomer.spread %>% 
  rename(`infused tracer` = infused.tracer,
         `labeled compound` = Compound) %>% 
  as.data.frame() %>%  
  write.xlsx(
    showNA = F,
    file = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/isotopologues.xlsx"
  )




# Calculate Fcirc of molecules
d.Fcirc = d.normalized.tidy %>% ungroup() %>% 
  filter(Compound == infused.tracer, C_Label == C_Label.max)


# Calculate artery-vein blood ratio of full labeled metabolite
d.art.vs.tail.ratio = d.Fcirc %>% 
  select(Compound, phenotype, blood, enrichment, infusion_mouseID, infusion_round) %>%
  spread(key = blood, value = enrichment)

# all glycerol tracing removed as artery blood has been removed after data import
d.art.vs.tail.ratio = d.art.vs.tail.ratio[complete.cases(d.art.vs.tail.ratio), ] %>% 
  mutate(correctFactor = art/tail) %>% 
  mutate(Compound = factor(Compound, levels = ordered.Compound, ordered = F))

# visualize correction factor
plt.art.tail.correctFactor = d.art.vs.tail.ratio %>%
  ggplot(aes(x = phenotype, y = correctFactor, color = phenotype)) + 
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) + 
  geom_beeswarm(size = 3, alpha = .5, cex = 5) +
  
  # geom_text(aes(label = infusion_mouseID)) +
  # expand_limits(y = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~Compound, scales = "free_y", nrow = 2) +
  stat_summary(fun.data = mean_se, fun.args = list(mult = 1),
               geom = "errorbar", width = .3, size = .6) +
  stat_summary(fun = mean, fun.args = list(mult = 1),
               geom = "crossbar", width = .3, size = .5) +
  
  labs(y = "Correction factor\n",
       title = "Vein uniform enrich x correct factor = artery enrich\nmean ± SEM") + 
  theme.mybw +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.title.x = element_blank(), legend.position = "None") +
  scale_color_manual(values = color.phenotype)

plt.art.tail.correctFactor 





# Check A-V ratio (uniform labeled) in detail
func.plt.tail.vs.art = function(data = d.normalized.tidy,  
                                myCompound.infused, myCompound.labeled, 
                                yAxis.input = "enrichment.normalized") {
  data %>%
    filter(infused.tracer == myCompound.infused) %>%
    filter(Compound == myCompound.labeled) %>%
    filter(C_Label == C_Label.max) %>% 
    ggplot(aes_string(x = "blood", y = yAxis.input)) +
    geom_point() +
    geom_line(aes(group = infusion_mouseID)) +
    facet_wrap(~phenotype) +
    labs(title = paste("Infused tracer:", myCompound.infused), 
         y = paste( ifelse(yAxis.input == "enrichment.normalized", 
                           "Normalized enrichment (%) of", "Enrichment (%) of"), myCompound.labeled) ) +
    scale_y_continuous(labels = function(x){x*100})
}

func.plt.tail.vs.art(data = d.normalized.tidy,
                     myCompound.infused = "Lactate", 
                     myCompound.labeled = "Lactate", yAxis.input = "enrichment") +
  expand_limits(y = 0)

# func.plt.tail.vs.art(data = d.normalized.tidy ,
#                      myCompound.infused = "Glucose", myCompound.labeled = "Glucose", yAxis.input = "enrichment")
# 
# 
# func.plt.tail.vs.art(data = d.normalized.tidy ,
#                      myCompound.infused = "Glutamine", myCompound.labeled = "Glutamine", yAxis.input = "enrichment")
# 
# func.plt.tail.vs.art(data = d.normalized.tidy ,
#                      myCompound.infused = "Alanine", myCompound.labeled = "Alanine", yAxis.input = "enrichment")
# 
# func.plt.tail.vs.art(data = d.normalized.tidy ,
#                      myCompound.infused = "Glycerol", myCompound.labeled = "Glycerol", yAxis.input = "enrichment")
# 
# func.plt.tail.vs.art(data = d.normalized.tidy ,
#                      myCompound.infused = "3-HB", myCompound.labeled = "3-HB", yAxis.input = "enrichment")
# 
# func.plt.tail.vs.art(data = d.normalized.tidy ,
#                      myCompound.infused = "C16:0", myCompound.labeled = "C16:0", yAxis.input = "enrichment")



# A-V correction factor mean value
d.art.vs.tail.ratio.summary = d.art.vs.tail.ratio %>%
  group_by(Compound, phenotype) %>%
  summarise(correctFactor.mean = mean(correctFactor),
            n.replicate = n_distinct(infusion_mouseID),
            correctFactor.SEM = sd(correctFactor) / sqrt(n.replicate) )
d.art.vs.tail.ratio.summary 

# export summary statistics as table: full label ratio of artery vs. vein
d.art.vs.tail.ratio.summary.output =  d.art.vs.tail.ratio.summary %>%
  mutate(ratio = str_c(correctFactor.mean %>% round(2), " ± ",
                       correctFactor.SEM %>% round(2), " (", n.replicate, ")")) %>%
  select(Compound, phenotype, ratio) %>%
  spread(Compound, ratio) 

# output organized molecular A-V ratio
# write.xlsx(d.art.vs.tail.ratio.summary.output %>% as.data.frame(), 
#            file = "/Users/boyuan/Desktop/Harvard/Research/db db mice/Infusion data/A-V ratio_raw.xlsx", 
#            sheetName = "Molecule ratio", append = T)


# read manually organized molecular A-V ratio
d.art.vs.tail.ratio.summary = read_excel("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/A-V ratio.xlsx", sheet = "Molecular ratio") %>% 
  gather(-phenotype, key = Compound, value = correctFactor.mean)  



# Use labeling preferentially only from artery blood; 
# if art not available, use tail blood, corrected with A-V ratio
d.Fcirc = d.Fcirc %>% select(-sample) %>% ungroup() %>% 
  spread(key = blood, value = enrichment) %>% 
  mutate(enrich = ifelse(is.na(art), tail, art), blood = ifelse(is.na(art), "tail", "art")) %>%  # if art not available, use tail blood (later corrected with A-V ratio)
  mutate(enrich = ifelse(is.na(enrich), art, enrich), blood = ifelse(is.na(blood), "art", blood)) %>%  # a few samples tail blood not available, and therefore use artery blood
  select(-c(art, tail))

# Calculate Fcirc using artery blood, and corrected-tail blood if artery blood is not available
d.Fcirc = d.Fcirc %>% 
  left_join(d.art.vs.tail.ratio.summary, by = c("Compound", "phenotype"))


# add new column of corrected tail blood for full labeled isotopomer; # for art blood: enrichment = enrich.corrected
d.Fcirc = d.Fcirc %>% 
  mutate(enrich.corrected = ifelse(blood == "art", enrich, enrich * correctFactor.mean)) %>% 
  select(-enrich)



# now calculate Fcirc, in unit of nmol / min; 
# at this step, both A-V sample pairs are included in this dataset; 
# later only artery blood if available is shown, or use corrected tail blood if artery blood is not available
d.Fcirc = d.Fcirc %>%
  mutate(Fcirc_animal = infusion_uL_perMin * tracer_conc_mM * (1-enrich.corrected) / enrich.corrected,
         Fcirc_g.BW = Fcirc_animal / BW) 




#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~





# Scrambled carbon atom enrichment =========
d.normalized.tidy = d.normalized.tidy %>% ungroup() %>% 
  mutate(enrich.isotopomer.contri = C_Label / C_Label.max * enrichment) 

d.enrich.atom = d.normalized.tidy %>% 
  group_by(Compound, infused.tracer,  sample, mouse_ID, infusion_round, blood, phenotype, infusion_mouseID) %>%
  summarise(enrich.atom = sum(enrich.isotopomer.contri, na.rm = T))



# calculate A-V enrichment ratio
d.enrich.atom = d.enrich.atom %>% 
  ungroup() %>% select(-sample) %>% 
  # since both A-V columns are from the same mouse - infusion round, sample column needs to be removed
  spread(key = blood, value = enrich.atom) %>% 
  mutate(correctFactor.atom = art / tail)

d.art.vs.tail.ratio.atom =  d.enrich.atom[complete.cases(d.enrich.atom), ] %>% 
  filter(is.finite(correctFactor.atom)) %>% # remove outliers
  filter(correctFactor.atom >= 0)


# remove outliers from 3 tracers in lean and obese mice
d.art.vs.tail.ratio.atom = d.art.vs.tail.ratio.atom %>% 
  filter(! (infusion_mouseID == "L_L6" & Compound == "Alanine")) %>%
  filter(! (infusion_mouseID == "m_L5" & Compound == "Alanine")) %>%
  filter(! (infusion_mouseID == "h_O2" & Compound == "Alanine")) %>%
  filter(! (infusion_mouseID == "h_O3" & Compound == "Alanine")) %>% 
  filter(! (infusion_mouseID == "m_d5" & Compound == "Glutamine")) %>%
  filter(! (infusion_mouseID == "m_L5" & Compound == "Glutamine")) %>%
  filter(! (infusion_mouseID == "k_O3" & Compound == "Glutamine")) %>%
  filter(! (infusion_mouseID == "k_d4" & Compound == "Glucose")) %>% 
  filter(! (infusion_mouseID == "ad_L8" & Compound == "Lactate"))


# Visualize A-V ratio in atomized labeling enrichment
func.plt.ratio.atom = function(labeled.Compound = "Glucose"){
  
  d.art.vs.tail.ratio.atom.selectedCompounds =  d.art.vs.tail.ratio.atom %>% 
    filter(infused.tracer %in% c("Glucose", "Lactate", "Alanine", "Glutamine")) %>% 
    filter(Compound %in% c("Glucose", "Lactate", "Alanine", "Glutamine", "Glycerol", "3-HB", "Acetate", "C16:0")) %>% 
    # show facet panels in order of infused tracer
    mutate(infused.tracer = factor(infused.tracer, levels =  ordered.Compound, ordered = F))
  
  
  d.art.vs.tail.ratio.atom.selectedCompounds %>% 
    # plot designated Compound across different infused tracers
    filter(Compound == labeled.Compound) %>% 
    
    ggplot(aes(x = phenotype, y = correctFactor.atom, color = phenotype, fill = phenotype)) +
    geom_hline(yintercept = 1, linetype = "dashed",  color = "steelblue") +
    geom_beeswarm(alpha = .3, size = 2, cex = 6) +
    geom_beeswarm(alpha = 1, size = 2, cex = 6, fill = NA, shape = 21) +
    # geom_text(aes(label = infusion_mouseID)) +
    # stat_summary(fun = mean, geom = "crossbar", width = .8) +
    # stat_summary(fun.data = mean_se, fun.args = list(mult = 1),
    #              size = .4) +
    stat_summary(fun.data = mean_se, fun.args = list(mult = 1),
                 geom = "errorbar", width = .3, size = .6) +
    stat_summary(fun = mean, fun.args = list(mult = 1),
                 geom = "crossbar", width = .3, size = .5) +
    
    facet_wrap( ~infused.tracer, scales = "free", nrow = 1) +
    scale_y_continuous(expand = expansion(mult = c(.1, .1))) +
    theme.mybw + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                       axis.ticks.length.x = unit(-.1, "cm"),
                       plot.title = element_text(hjust = .5),
                       strip.text = element_text(size = 11),
                       panel.spacing = unit(1, "lines"),
                       strip.placement = "outside") +
    labs(y = labeled.Compound) +
    scale_color_manual(values = color.phenotype) +
    scale_fill_manual(values = color.phenotype) 
}

# Vein blood atomized enrich x correct factor = art atomized enrich, mean ± SEM"
# left : labeled Compound; top strip: infused tracer
plt.art.tail.correctFactor.atomized = 
  plot_grid(func.plt.ratio.atom("Glucose") + coord_cartesian(ylim = c(.9, 1.2)),
            func.plt.ratio.atom("Lactate") + coord_cartesian(ylim = c(1, 2.5)),
            func.plt.ratio.atom("Alanine") + coord_cartesian(ylim = c(.8, 2)),
            func.plt.ratio.atom("Glutamine") + coord_cartesian(ylim = c(.9, 1.4)),
            nrow = 4
  )
plt.art.tail.correctFactor.atomized





# Summary of atomized ratio * important dataset
d.art.vs.tail.ratio.atom.summary = d.art.vs.tail.ratio.atom %>% 
  group_by(infused.tracer, Compound, phenotype) %>% 
  summarise(correctFactor.atom.mean = mean(correctFactor.atom, na.rm = T),
            n.replicate = n_distinct(infusion_mouseID),
            correctFactor.atom.SEM = sd(correctFactor.atom, na.rm = T) / sqrt(n.replicate) )




#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-





# export summary statistics for publication
d.art.vs.tail.ratio.atom.summary.paper = d.art.vs.tail.ratio.atom.summary %>% 
  filter(infused.tracer %in% c("Glucose", "Lactate", "Alanine", "Glutamine" )) %>% 
  filter(Compound %in% c("Glucose", "Lactate", "Alanine", "Glutamine")) %>% 
  mutate(atomized.ratio = str_c(correctFactor.atom.mean %>% round(digits = 2), " ± ",  
                                correctFactor.atom.SEM %>% round(digits = 2),
                                " (", n.replicate, ")")) %>% 
  select(infused.tracer, Compound, phenotype, atomized.ratio) %>% 
  spread(Compound, atomized.ratio) %>% 
  select(infused.tracer, phenotype, Glucose, Lactate, Alanine, Glutamine) %>% 
  mutate(infused.tracer = factor(infused.tracer, levels = ordered.Compound, ordered = F)) %>% arrange(infused.tracer) 

d.art.vs.tail.ratio.atom.summary.paper 



# Visualize the A-V atomized labeling enrichment as a heatmap
kkk = d.art.vs.tail.ratio.atom.summary %>% 
  filter(Compound %in% c("Glucose", "Lactate", "Glutamine", "Alanine", "Glycerol")) %>% 
  mutate(infused.tracer = factor(infused.tracer, levels = ordered.Compound, ordered = F ),
         Compound = factor(Compound, levels = ordered.Compound, ordered = F),
         phenotype = factor(phenotype, levels = c("WT", "ob/ob", "db/db"), ordered = F)) %>% 
  arrange(Compound, phenotype) %>% 
  mutate(Compound_phenotype = str_c(Compound, ", ", phenotype))

kkk$Compound_phenotype = factor(kkk$Compound_phenotype, levels = unique(kkk$Compound_phenotype) %>% rev(), ordered = F)

kkk.selected = kkk %>% 
  filter(Compound %in% c("Glucose", "Lactate", "Glutamine", "Alanine")) %>% 
  filter(infused.tracer %in% c("Glucose", "Lactate", "Glutamine", "Alanine"))

plt.art.tail.correctFactor.atomized.heatmap = 
  kkk.selected %>% 
  ggplot(aes(x = infused.tracer, y = Compound_phenotype, fill = correctFactor.atom.mean)) +
  geom_tile(linetype = 1, color = "white", size = .7) +
  geom_text(data = kkk.selected, aes(label = round(correctFactor.atom.mean, 2)), color = "snow1",  fontface = "bold") +
  geom_text(data = kkk.selected %>% filter(correctFactor.atom.mean >= 1.9),
            aes(label = round(correctFactor.atom.mean, 2)), color = "black", fontface = "bold") +
  # scale_fill_gradient(low = "white", high = "firebrick") +
  # scale_fill_gradient2(low = "steelblue", mid = "white", high = "brown",
  #                      midpoint = .9, 
  #                      breaks = seq(1, 2, .2), name = "") +
  scale_fill_viridis(option = "A") +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(mult = c(0, 0)))   +
  theme.mybw  + theme(legend.key.height = unit(1.2, "cm"), legend.title = element_blank()) +
  labs(y = "labeled metabolite\n", x = "infused tracer",
       title = "A-V Correction factor\nvein enrichment * factor = artery enrichment") 

plt.art.tail.correctFactor.atomized.heatmap




# Output the A-V ratio of the gluconeogenic set compound; manually set ratio = 1 for fatty acid + acetate series
setwd("/Users/boyuan/Desktop/Harvard/Research/db db mice/Infusion data")

d.art.vs.tail.ratio.atom.summary.output = 
  d.art.vs.tail.ratio.atom.summary %>% 
  filter(Compound %in% ordered.gluconeogenic.Set[-3]) %>% 
  filter(infused.tracer %in% ordered.gluconeogenic.Set[-3]) %>%
  select(infused.tracer, Compound, phenotype, correctFactor.atom.mean) %>% 
  spread(Compound, correctFactor.atom.mean)

write.xlsx(d.art.vs.tail.ratio.atom.summary.output %>% as.data.frame(), 
           file = "A-V ratio_raw.xlsx", sheetName = "Gluconeogenic set")
write.xlsx(d.art.vs.tail.ratio.atom.summary.paper %>% as.data.frame(), 
           file = "A-V ratio_raw.xlsx", sheetName = "Gluconeogenic set_stats", append = T)



# import manually organized A-V ratio dataset
d.A.V.ratio.allSet = read_excel("/Users/boyuan/Desktop/Harvard/Research/db db mice/Infusion data/A-V ratio.xlsx")
d.A.V.ratio.allSet.tidy = d.A.V.ratio.allSet %>% 
  gather(-c(1:2), key = Compound, value = correctFactor.atom.mean )


# Fcirc based on carbon atoms
d.enrich.atom = d.enrich.atom %>% 
  
  # correctFactor.atom is specific to each A-V sample pair, 
  # already used above to calculate the mean A-V ratio, no longer useful, thus removed here
  select(-correctFactor.atom) %>% 
  
  left_join(d.A.V.ratio.allSet.tidy, 
            by = c("infused.tracer", "Compound", "phenotype")) %>% 
  mutate(tail.corrected = tail * correctFactor.atom.mean) %>% 
  
  # when art blood is available, only use artery blood enrichment without correction;
  # otherwise use A-V ratio-corrected tail venous blood; and integrate both A and V blood into single column
  # NOTE that enrich.atom.corrected includes: original arterial blood enrichment, or venous blood enrichment after A-V ratio correction
  mutate(enrich.atom.corrected = ifelse(!is.na(art), art, tail.corrected)) %>% 
  select(-c(art, tail, tail.corrected))



d.Fcirc.standard.atom = d.enrich.atom %>% 
  filter(infused.tracer == Compound) %>% 
  # left_join(d.infusion_rounds, by = c("infusion_round", "mouse_ID", "infused.tracer")) %>% 
  left_join(d.Fcirc, by = c("Compound",  "infused.tracer",  "mouse_ID", "infusion_round", "infusion_mouseID", "phenotype")) %>% 
  mutate(Fcirc_atom.animal = tracer_conc_mM * infusion_uL_perMin * C_Label.max * (1-enrich.atom.corrected) / enrich.atom.corrected, # nmol carbon atoms / min/animal
         Fcirc_atom.g.BW = Fcirc_atom.animal / BW) # nmol carbons atoms/min/g BW

# Mark mice ID that has both tail venous and art blood; use only art blood in this case
infusion_mouseID.havingBoth_art_tail = d.art.vs.tail.ratio$infusion_mouseID # (or use d.art.vs.tail.ratio.atom data)
d.Fcirc.standard.atom = d.Fcirc.standard.atom  %>% 
  mutate(both.art.tail = infusion_mouseID %in% infusion_mouseID.havingBoth_art_tail)

d.Fcirc.standard.atom_art.or.tail = d.Fcirc.standard.atom %>% 
  filter(! (both.art.tail == T & blood == "tail"))

# summary stats
d.Fcirc_standard.atom.summary = d.Fcirc.standard.atom_art.or.tail %>% # filter(Compound == "3-HB")
  group_by(Compound, phenotype) %>%
  summarise( 
    n.rep = n_distinct(infusion_mouseID),
    Fcirc_animal.mean = mean(Fcirc_animal, na.rm = T),
    Fcirc_animal.sd = sd(Fcirc_animal, na.rm = T),
    
    Fcirc_g.BW.mean = mean(Fcirc_g.BW, na.rm = T),
    Fcirc_g.BW.sd = sd(Fcirc_g.BW, na.rm = T),
    Fcirc_g.BW.sem = Fcirc_g.BW.sd / sqrt(n.rep),
    
    Fcirc_animal.atom.mean = mean(Fcirc_atom.animal, na.rm = T),
    Fcirc_animal.atom.sd = sd(Fcirc_atom.animal, na.rm = T),
    
    Fcirc_g.BW.atom.mean = mean(Fcirc_atom.g.BW, na.rm = T),
    Fcirc_g.BW.atom.sd = sd(Fcirc_atom.g.BW, na.rm = T) )

d.Fcirc_standard.atom.summary[, 1:5] %>% 
  mutate(Fcirc_animal.sd / Fcirc_animal.mean * 100)




# plot Fcirc
d.Fcirc.standard.atom_art.or.tail = d.Fcirc.standard.atom_art.or.tail %>%  # plot in order of infused tracer
  mutate(infused.tracer = factor(infused.tracer, levels = ordered.Compound, ordered = F ))

o <- d.Fcirc.standard.atom_art.or.tail %>% 
  filter(phenotype %in% c("WT", "ob/ob", "HFD"))

func.plt.Fcirc.Compound = function(yAxis = "Fcirc_animal") {
  o %>% 
    # filter(infused.tracer == "C18:1") %>% 
    mutate(phenotype = factor(phenotype, levels = ordered.phenotype, ordered = F)) %>% 
    ggplot(aes_string(x = "phenotype", y = yAxis, color = "phenotype", fill = "phenotype")) +
    
    stat_summary(fun = mean, fun.args = list(mult = 1),
                 geom = "bar", width = 1, color = "black") +
    
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
                 geom = "errorbar",width = .5, color = "black") +
    
    geom_beeswarm(size = 1.5, cex = 3, alpha = 1, show.legend = F,
                  fill = "white", color = "black", shape = 21 ) + # base layer
    geom_beeswarm(size = 1.5, cex = 3,  shape = 21,
                  alpha = .4, show.legend = F, color = "black" ) + # top layer
    
    # geom_text(aes(label = infusion_mouseID), size = 5, color = "black",
    #           position = position_jitter(.2, 0)) +
    
    facet_wrap(~infused.tracer, scales = "free", nrow = 2) +
    labs() + theme.myClassic +
    theme(axis.title.x = element_blank(),
          axis.ticks.length.x = unit(-0.1, "cm"),
          
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          panel.spacing = unit(10, "mm")
    ) +
    expand_limits(y = 0) +
    scale_y_continuous(expand = expansion(mult = c(0, .4))) +
    scale_color_manual(values = color.phenotype) +
    scale_fill_manual(values = color.phenotype,
                      labels = function(x){str_replace(x, "WT", "chow")})  +
    scale_x_discrete(expand = expansion(mult = c(.5, .5))) 
}


plt.Fcirc.perGram_BW = func.plt.Fcirc.Compound(yAxis = "Fcirc_g.BW") 
plt.Fcirc.perGram_BW

ggsave(filename = "molecule per gBW Fcirc.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 5.5 * 1.1, width = 9 * 1.1)

plt.Fcirc_animal = func.plt.Fcirc.Compound(yAxis = "Fcirc_animal") +
  labs(y = "nmol molecules / min /animal\n")

plt.Fcirc_animal 



# add statistic test
d.Tukey.animal <- o %>% 
  nest(-infused.tracer) %>% 
  mutate(model = map(data, ~aov(formula = Fcirc_animal ~ phenotype, data = .)),
         tukey = map(model, ~TukeyHSD(.))) %>% 
  # extract statistics 
  mutate(stats = map(tukey, ~pluck(., "phenotype") %>% 
                       as.data.frame() %>% rownames_to_column(var = "group"))) %>% 
  unnest(stats) %>% 
  select(infused.tracer, group, `p adj`)

d.Tukey.animal

# calculate y position
d.Tukey.animal2 <- o %>% 
  # calculate max Fcirc value per comparision group
  group_by(infused.tracer, phenotype) %>% 
  summarise(max = max(Fcirc_animal)) %>% spread(phenotype, max) %>% 
  mutate(`WT-HFD` = max(HFD, WT),
         `ob/ob-HFD` = max(HFD, `ob/ob`),
         `WT-ob/ob` = max(`ob/ob`, WT, HFD) * 1.15,
         .keep = "unused") %>% 
  pivot_longer(-infused.tracer, names_to = "group", values_to = "y") %>% 
  # combine with the tukey data
  left_join(d.Tukey.animal, ., by = c("group", "infused.tracer")) %>% 
  # calculate x axis
  separate(group, into = c("g1", "g2"), sep = "-") %>% 
  mutate(across(c(g1, g2), ~factor(., levels = ordered.phenotype, exclude = "db/db")),
         x.star = (as.numeric(g1) + as.numeric(g2))/2)

# calculate significant stars
func.generate_stars <- function(p_value) {
  if (p_value < 0.0001) {
    return("****")
  } else if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("")
  }
}

d.Tukey.animal3 <- d.Tukey.animal2 %>% 
  mutate(stars = map_chr(`p adj`, func.generate_stars))

# add stars to the plot
plt.Fcirc_animal +
  geom_segment(data = d.Tukey.animal3 %>% filter( `p adj` < 0.05), 
               aes(x = g1, xend = g2, y = y*1.13, yend = y*1.13), inherit.aes = F) +
  geom_text(data = d.Tukey.animal3 %>% filter( `p adj` < 0.05),
            aes(label = stars, 
                x = x.star, 
                y = y * 1.15),
            size = 5,
            inherit.aes = F) +
  scale_y_continuous(labels = ~.x / 1000, 
                     name = "circulatory turnover flux (µmol molecules / min / animal)\n",
                     expand = expansion(mult = c(0, .1))) +
  theme(legend.position = "bottom")

ggsave(filename = "molecule Fcirc.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 5.5 * 1.1, width = 9 * 1.1)

# plt.Fcirc_animal + theme(panel.spacing.y = unit(20, "pt"))


plt.Fcirc.perGram_BW.atom = func.plt.Fcirc.Compound(yAxis = "Fcirc_atom.g.BW") 
plt.Fcirc_animal.atom = func.plt.Fcirc.Compound(yAxis = "Fcirc_atom.animal") 


plt.Fcirc = plot_grid(plt.Fcirc.perGram_BW + labs(y = "nmol molecules / min / g" ),
                      plt.Fcirc_animal + labs(y = "nmol molecules / min / animal" ) , nrow = 2)
# plt.Fcirc

plt.Fcirc.atom = plot_grid(plt.Fcirc.perGram_BW.atom + labs(y = "nmol C-atoms / min / g" ), 
                           plt.Fcirc_animal.atom  + labs(y = "nmol C-atoms / min / animal" ), nrow = 2)
# plt.Fcirc.atom



# add statistic test
d.Tukey.animal.Fcirc_atom <- o %>% 
  nest(-infused.tracer) %>% 
  mutate(model = map(data, ~aov(formula = Fcirc_atom.animal ~ phenotype, data = .)),
         tukey = map(model, ~TukeyHSD(.))) %>% 
  # extract statistics 
  mutate(stats = map(tukey, ~pluck(., "phenotype") %>% 
                       as.data.frame() %>% rownames_to_column(var = "group"))) %>% 
  unnest(stats) %>% 
  select(infused.tracer, group, `p adj`)

d.Tukey.animal.Fcirc_atom

# calculate y position
d.Tukey.animal.Fcirc_atom2 <- o %>% 
  # calculate max Fcirc value per comparision group
  group_by(infused.tracer, phenotype) %>% 
  summarise(max = max(Fcirc_atom.animal)) %>% spread(phenotype, max) %>% 
  mutate(`WT-HFD` = max(HFD, WT),
         `ob/ob-HFD` = max(HFD, `ob/ob`),
         `WT-ob/ob` = max(`ob/ob`, WT, HFD) * 1.15,
         .keep = "unused") %>% 
  pivot_longer(-infused.tracer, names_to = "group", values_to = "y") %>% 
  # combine with the tukey data
  left_join(d.Tukey.animal.Fcirc_atom, ., by = c("group", "infused.tracer")) %>% 
  # calculate x axis
  separate(group, into = c("g1", "g2"), sep = "-") %>% 
  mutate(across(c(g1, g2), ~factor(., levels = ordered.phenotype, exclude = "db/db")),
         x.star = (as.numeric(g1) + as.numeric(g2))/2)


d.Tukey.animal.Fcirc_atom3 <- d.Tukey.animal.Fcirc_atom2 %>% 
  mutate(stars = map_chr(`p adj`, func.generate_stars))

# add stars to the plot
plt.Fcirc_animal.atom +
  geom_segment(data = d.Tukey.animal.Fcirc_atom3 %>% filter( `p adj` < 0.05), 
               aes(x = g1, xend = g2, y = y*1.13, yend = y*1.13), inherit.aes = F) +
  geom_text(data = d.Tukey.animal.Fcirc_atom3 %>% filter( `p adj` < 0.05),
            aes(label = stars, 
                x = x.star, 
                y = y * 1.15),
            size = 5,
            inherit.aes = F) +
  scale_y_continuous(labels = ~.x / 1000, 
                     name = "µmol C / min / animal\n",
                     expand = expansion(mult = c(0, .1))) +
  theme(legend.position = "bottom")

ggsave(filename = "Fcirc atom per animal.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 5.5 * 1.1, width = 9 * 1.1)



d.Fcirc.standard.atom_art.or.tail %>% 
  group_by(Compound, phenotype) %>% 
  summarise(Fcirc_g.BW.mean = mean(Fcirc_g.BW) %>% round(1),
            Fcirc_g.BW.SD = sd(Fcirc_g.BW) %>% round(1) ) %>% 
  mutate(Fcirc = paste0(Fcirc_g.BW.mean,  "±", Fcirc_g.BW.SD), .keep = "unused") %>% 
  spread(phenotype, Fcirc) %>% 
  select(Compound, WT)



# # add significant stars
# d.pairwise.T.Fcirc.molecule.animal <-  d.Fcirc.standard.atom_art.or.tail %>% 
#   filter(! infused.tracer %in% c("3-HB", "C18:2")) %>% 
#   # filter(infused.tracer == "Glycerol") %>% 
#   mutate(phenotype = factor(phenotype, levels = ordered.phenotype, ordered = F)) %>% 
#   filter(phenotype %in% c("WT", "ob/ob", "HFD")) %>% 
#   select(infused.tracer, phenotype, Fcirc_animal ) 
# 
# 
# d.pairwise.T.Fcirc.molecule.animal.stats <-  d.pairwise.T.Fcirc.molecule.animal %>% 
#   group_by(infused.tracer) %>% 
#   pairwise_t_test(Fcirc_animal ~ phenotype, paired = F, 
#                   p.adjust.method = "bonferroni") %>% 
#   select(-c(n1, n2, p, p.signif)) %>% 
#   add_xy_position(x = "phenotype") 
# 
# 
# d.pairwise.T.Fcirc.molecule.animal.stats %>% 
#   filter(p.adj < 0.05)
# 
# d.pairwise.T.Fcirc.molecule.animal.stats %>% 
#   filter(infused.tracer == "Glutamine")


# Plot bar plot of accumulated carbon atoms
# Calculate the y-axis position of accumulated error bar
func.cumulatedSum = function(vector){
  cc = vector()
  for (i in length(vector):1){ cc[i] <- sum(vector[length(vector):i]) }
  return(cc)
}


d.Fcirc_standard.atom.summary = d.Fcirc_standard.atom.summary %>% 
  mutate(Compound = factor(Compound, levels = ordered.Compound %>% rev(), ordered = F)) %>% 
  arrange(Compound) %>% 
  group_by(phenotype) %>% 
  mutate(position.y.error_Fcic.g.BW.atom = func.cumulatedSum(Fcirc_g.BW.atom.mean),
         position.y.error_Fcic.animal.atom = func.cumulatedSum(Fcirc_animal.atom.mean))



# plot 
func.plt.Fcirc.atom.bar = function(yAxis){
  p = d.Fcirc_standard.atom.summary %>% 
    mutate(phenotype = factor(phenotype, levels = ordered.phenotype)) %>% 
    filter(phenotype != "db/db") %>% 
    ggplot(aes_string(x = "phenotype", y = yAxis, fill = "Compound")) +
    geom_bar(stat = "identity", position = "stack", alpha = .7,color = "black") +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    # scale_color_nejm() +
    # scale_fill_nejm() +
    theme.mybw + 
    theme(axis.title.x = element_blank(), 
          legend.title = element_blank()) +
    scale_x_discrete(expand = expansion(add = .8),
                     labels = function(x){str_replace(x, "WT", "chow")})
  
  if (yAxis == "Fcirc_animal.atom.mean") { 
    p = p + geom_errorbar(aes(ymin =  position.y.error_Fcic.animal.atom - Fcirc_animal.atom.sd,
                              ymax =  position.y.error_Fcic.animal.atom ),
                          width = .2) +
      labs(y = "nmol C-atoms / min / animal") 
  }
  
  if (yAxis == "Fcirc_g.BW.atom.mean") { 
    p = p + geom_errorbar(aes(ymin =  position.y.error_Fcic.g.BW.atom - Fcirc_g.BW.atom.sd,
                              ymax =  position.y.error_Fcic.g.BW.atom ),
                          width = .2)  +
      labs(y = "nmol C-atoms / min / g BW") 
  } 
  return(p)
}
p1 = func.plt.Fcirc.atom.bar(yAxis = "Fcirc_g.BW.atom.mean" )
p2 = func.plt.Fcirc.atom.bar(yAxis = "Fcirc_animal.atom.mean")
plt.Fcirc.atom.bar = plot_grid(p1, p2, nrow = 1)
plt.Fcirc.atom.bar




#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~#~~~



# normalized labeling 
d.enrich.atom = d.enrich.atom %>% 
  mutate(infused.tracer = factor(infused.tracer, levels = ordered.Compound, ordered = F),
         Compound = factor(Compound, levels = ordered.Compound, ordered = F))

d.enrich.atom.selected = d.enrich.atom %>% 
  filter(Compound %in% ordered.Compound[c(1:6, 8:11)]) %>% 
  filter(infused.tracer %in% ordered.Compound[c(1:6, 8:11)]) %>% 
  select(-c(correctFactor.atom.mean))

myCompounds = d.enrich.atom.selected$Compound %>% unique()
names(myCompounds) = myCompounds


# myCompounds = "Glucose"
d.enrich.atom.normalized = tibble()

for (tracer in myCompounds) {
  # select each infused tracer
  d.i = d.enrich.atom.selected %>% 
    filter(infused.tracer == tracer) %>% 
    spread(Compound, enrich.atom.corrected)
  
  # add tracer enrich column
  d.i = d.i %>% mutate(tracer = d.i[[tracer]]) 
  
  # normalize each labeled metabolite
  for (Compound.i in myCompounds){
    c.i = d.i[[Compound.i]] / d.i[[tracer]]
    d.i = cbind(d.i, c.i) %>% as_tibble()
    colnames(d.i)[ncol(d.i)] = paste0(Compound.i, "_norm")
  }
  d.enrich.atom.normalized = rbind(d.enrich.atom.normalized, d.i)
}

d.enrich.atom.normalized = d.enrich.atom.normalized %>% 
  select(infused.tracer,phenotype, infusion_mouseID, contains("norm"))

d.enrich.atom.normalized.tidy = d.enrich.atom.normalized %>% 
  gather(-c(infused.tracer, phenotype, infusion_mouseID), key = Compound, value = enrich.normalized) %>% 
  mutate(Compound = str_remove(Compound, pattern = "_norm")) %>% 
  filter(enrich.normalized <=1) %>% 
  mutate(Compound = factor(Compound, levels = ordered.Compound, ordered = F))

# summary
d.enrich.atom.normalized.summary = d.enrich.atom.normalized.tidy %>% 
  group_by(phenotype, infused.tracer, Compound) %>% 
  summarise(enrich.normalized.mean = mean(enrich.normalized, na.rm = T ),
            enrich.normalized.sd = sd(enrich.normalized, na.rm = T),
            # calculate standard error
            n.replciate = n_distinct(infusion_mouseID),
            enrich.normalized.sem = enrich.normalized.sd / sqrt(n.replciate)) 
d.enrich.atom.normalized.summary

d.enrich.atom.normalized.summary$infused.tracer %>% unique()


# Normalized labeling
plt.labeling.normalized = 
  d.enrich.atom.normalized.summary %>% 
  # filter(phenotype != "db/db") %>% 
  
  # filter(Compound %in% c("Glucose", "Lactate", "Glycerol")) %>% 
  ggplot(aes(x = Compound, y = enrich.normalized.mean, fill = phenotype, color = phenotype)) +
  geom_col(alpha = .4, position = "dodge", width = .8) +
  facet_wrap(~infused.tracer, nrow = 4) +
  geom_errorbar(aes(ymin = enrich.normalized.mean - enrich.normalized.sd,
                    ymax = enrich.normalized.mean + enrich.normalized.sd), 
                size = 1, width = .5, position = position_dodge(.8)) +
  
  # geom_point(data = d.enrich.atom.normalized.tidy %>%
  #              filter(infused.tracer == Compound & phenotype != "db/db"), # & Compound != infusedTracer
  #            aes(y = enrich.normalized),  alpha = .6, position = position_dodge(0)) +
  # 
  # geom_point(data = d.enrich.atom.normalized.tidy %>%
  #              filter(infused.tracer != Compound& phenotype != "db/db"), # & Compound != infusedTracer
  #            aes(y = enrich.normalized),  alpha = .6,  position = position_dodge(.8)) +
  
  labs(title = "Normalized labeling by infused tracer", 
       x = "\nLabeled metabolite", 
       y = paste("Normalized labeling \n" )) +
  scale_color_manual(values = color.phenotype, labels = function(x){str_replace(x, "WT", "chow")}) +
  scale_fill_manual (values = color.phenotype, labels = function(x){str_replace(x, "WT", "chow")}) +
  scale_y_continuous(expand = expansion(mult = c(0, .05)), 
                     breaks = seq(0, 1, .2)) +
  coord_cartesian(ylim = c(0, NA)) +
  theme.myClassic + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.ticks.length.x = unit(-.1, "cm"),
        panel.spacing = unit(.5, "lines"),
        strip.text = element_text(size = 13),
        legend.position = c(.8, .04))


plt.labeling.normalized <- plt.labeling.normalized +
  geom_point(data = d.enrich.atom.normalized.tidy, # %>% filter(phenotype != "db/db"),
             aes(y = enrich.normalized,  x = Compound),
             position = position_dodge(.8)) # +
# coord_cartesian(ylim = c(0, .4)) +
# geom_hline(yintercept = .2) 

plt.labeling.normalized   
# scale_y_continuous(expand = expansion(mult = c(0, .01)),
#                    breaks = seq(0, .04, .005)) +
# coord_cartesian(ylim = c(0, .01))


# heatmap
d.enrich.atom.normalized.summary = d.enrich.atom.normalized.summary %>% 
  mutate(tracer.pheno = str_c(phenotype, " ", infused.tracer)) %>% 
  arrange(infused.tracer, phenotype)

# put in current order
d.enrich.atom.normalized.summary$tracer.pheno = 
  factor(d.enrich.atom.normalized.summary$tracer.pheno,
         levels = d.enrich.atom.normalized.summary$tracer.pheno %>% unique() %>% rev(), 
         ordered = F)

# plot
color.enrich = c(viridis::viridis(20, option = "F"), "ivory") %>% rev()
# color.enrich %>% scales::show_col()

plt.heatmap.normalized.enrich = d.enrich.atom.normalized.summary %>%  # good
  filter(phenotype != "db/db") %>% 
  ggplot(aes(x = Compound, y = tracer.pheno, fill = enrich.normalized.mean)) +
  geom_tile(color = "white") +
  # scale_fill_viridis_c(option = "F", direction = -1) +
  scale_fill_gradientn(colors = color.enrich, breaks = seq(0, 1, .2)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  coord_fixed(ratio = .8) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  guides(fill = guide_colorbar(barheight = unit(200, "pt")))

plt.heatmap.normalized.enrich






# -~# -~# -~# -~# -~# -~# -~# -~# -~# -~# -~# -~# -~# -~# -~# -~# -~# -~# -~
# Original enrichment before any normalization

d.enrich.atom.summary = d.enrich.atom.selected %>% 
  group_by(phenotype, infused.tracer, Compound) %>% 
  summarise(enrich.original.mean = mean(enrich.atom.corrected, na.rm = T ),
            enrich.original.sd = sd(enrich.atom.corrected, na.rm = T),
            n.replicate = n_distinct(infusion_mouseID),
            enrich.original.sem = enrich.original.sd / sqrt(n.replicate)) %>% 
  mutate(phenotype = factor(phenotype, levels = ordered.phenotype))
d.enrich.atom.summary


# Plot 1: WT

plt.original.labeling.WT <-  d.enrich.atom.summary %>% 
  filter(phenotype == "WT") %>% 
  ggplot(aes(x = Compound, y = enrich.original.mean, 
             fill = phenotype, color = phenotype)) +
  # bar
  geom_col(alpha = .6, position = position_dodge(.8), width = .8,
           color = "black", fill = "snow4") +
  # point
  geom_quasirandom(
    data = d.enrich.atom.selected %>% 
      # filter(infused.tracer == Compound)  %>% 
      filter(phenotype == "WT"), 
    aes(y = enrich.atom.corrected), size = .8, color = "black") +
  # geom_text(data = d.enrich.atom.selected %>% 
  #             # filter(infused.tracer == Compound)  %>% 
  #             filter(phenotype == "WT"), 
  #           aes(label = infusion_mouseID, y = enrich.atom.corrected)) +
  
  # errorbar
  geom_errorbar(
    aes(ymin = enrich.original.mean - enrich.original.sd,
        ymax = enrich.original.mean + enrich.original.sd),
    linewidth = .5, width = .5, position = position_dodge(.8),
    color = "black") +
  
  facet_wrap(~infused.tracer, ncol = 3, scales = "free") +
  
  scale_color_manual(values = color.phenotype) +
  scale_fill_manual(values = color.phenotype) +
  scale_y_continuous(expand = expansion(mult = c(0, .05)), 
                     n.breaks = 5) +
  labs(title = "Original labeling", 
       y = "Labeling enrichment in serum\n", 
       x = "Labeled metabolite") +
  theme.myClassic + 
  theme(axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 10),
        panel.spacing = unit(1, "lines"),
        legend.position = "none") 

plt.original.labeling.WT

ggsave(filename = "original labeling WT.pdf", 
       plot = plt.original.labeling.WT, 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 12, width = 10)



# Plot 2: WT vs. HFD vs. ob/ob

plt.original.labeling.3phenotype <- d.enrich.atom.summary %>% 
  filter(phenotype != "db/db") %>% 
  ggplot(aes(x = Compound, y = enrich.original.mean, fill = phenotype, color = phenotype)) +
  geom_col(alpha = .8, position = position_dodge(.8),
           width = .8, color = "black", linewidth = .3) +
  
  geom_errorbar(aes(ymin = enrich.original.mean - enrich.original.sd,
                    ymax = enrich.original.mean + enrich.original.sd),
                linewidth = .5, width = .5, position = position_dodge(.8),
                color = "black") +
  
  geom_quasirandom(
    data = d.enrich.atom.selected %>% filter(phenotype != "db/db") %>% 
      mutate(phenotype = factor(phenotype, levels = ordered.phenotype)), 
    aes(y = enrich.atom.corrected), 
    shape = 21, size = .7, alpha = .6, dodge.width = .8, color = "black",
    show.legend = F) +
  
  facet_wrap(~infused.tracer, ncol = 2, scales = "free") +
  
  scale_color_manual(values = color.phenotype, labels = function(x){str_replace(x, "WT", "chow")}) +
  scale_fill_manual (values = color.phenotype, labels = function(x){str_replace(x, "WT", "chow")}) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), 
                     n.breaks = 5) +
  theme.myClassic + 
  theme(axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 10),
        panel.spacing = unit(1, "lines"),
        legend.position = "bottom")  +
  labs(title = "Original labeling", 
       y = "Labeling enrichment in serum\n", 
       x = "Labeled metabolite")

plt.original.labeling.3phenotype

ggsave(filename = "original labeling 3 phenotypes.pdf", 
       plot = plt.original.labeling.3phenotype, 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 12, width = 9)




# Heatmap of original labeling

d.enrich.atom.summary = d.enrich.atom.summary %>% 
  mutate(tracer.pheno = str_c(phenotype, " ", infused.tracer)) %>% 
  arrange(infused.tracer, phenotype)

# put in current order
d.enrich.atom.summary$tracer.pheno = 
  factor(d.enrich.atom.summary$tracer.pheno,
         levels = d.enrich.atom.summary$tracer.pheno %>% unique() %>% rev(), 
         ordered = F)


color.enrich.original = c(viridis::viridis(20, option = "G"), "white") %>% rev()
# color.enrich.original  <- brewer.pal(11, "Spectral") %>% rev()
color.enrich.original <- colorRampPalette(color.enrich.original, bias = 1)(20)

plt.heatmap.original.enrich.WT = d.enrich.atom.summary %>% 
  filter(phenotype == "WT") %>% 
  ggplot(aes(x = Compound, y = tracer.pheno, fill = enrich.original.mean)) +
  geom_tile(color = "white") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12)) +
  coord_fixed(ratio = .9, expand = 0) +
  scale_fill_gradientn(colors = color.enrich.original, 
                       breaks = seq(0, .4, .05),
                       guide = guide_colorbar(barheight = unit(150, "pt"))) +
  scale_y_discrete(labels = function(x){str_remove(x, pattern = "WT ")})

plt.heatmap.original.enrich.WT


ggsave(filename = "heatmap interlabeling WT.pdf", 
       plot = plt.heatmap.original.enrich.WT, 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 4*1.2, width = 4*1.3)




# check number of replicates
r <- d.enrich.atom.selected %>% filter(phenotype != "db/db") %>% 
  group_by(phenotype, infused.tracer, Compound) %>% 
  summarise(n = n()) %>% 
  spread(infused.tracer, n)
r


# - Generalized linear regression analysis
d.Fcirc.BW <- d.Fcirc.standard.atom_art.or.tail %>% 
  select(Compound, phenotype, infusion_mouseID, infusion_round, BW, Fcirc_animal) %>% 
  filter(phenotype != "db/db") %>% 
  # arrange phenotype in order (WT as reference)
  mutate(phenotype = factor(phenotype, levels = ordered.phenotype)) %>% 
  arrange(phenotype) %>% 
  mutate(Compound = factor(Compound, levels = ordered.Compound))


# remove rounds (during paper revision) designed to compare fluxes between HFD and ob/ob of matched body weight  
# these rounds contain HFD with high body mass (50-70g) with 6-12 months on the high-fat diet
# these matched body weight fluxes are otherwise displayed as bars stratified by body weight range 
rounds.matched.BW <- c("cg", "da", "dc", "dd", "ci")

plt.Fcirc.BW <- d.Fcirc.BW %>% 
  
  filter(! infusion_round %in% rounds.matched.BW) %>%
  
  ggplot(aes(x = BW, y = Fcirc_animal, color = phenotype, fill = phenotype)) +
  geom_point(size = 3, stroke = 1, alpha = .4) +
  facet_wrap(~Compound, scales = "free", nrow = 2) +
  geom_smooth(method = "lm", se = F, show.legend = F) +
  scale_color_manual(values = color.phenotype) +
  scale_fill_manual(values = color.phenotype) +
  # convert nmol to µmol / min / animal
  scale_y_continuous(labels = function(x){x/1000}) +
  labs(y = "µmol molecules / min") +
  theme(legend.position = "bottom",
        panel.spacing = unit(7, units = "pt"),
        axis.text = element_text(size = 16),
        strip.text = element_text(size = 15))  +
  ggalt::geom_encircle(s_shape = 1, expand = 0, alpha = .15, show.legend = F)

plt.Fcirc.BW

ggsave(filename = "Fcirc vs body weight.pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 6, width = 12)


# plot body-weight matched HFD and ob/ob
d.Fcirc.BW.matched <- d.Fcirc.BW %>% 
  filter(phenotype %in% c("ob/ob", "HFD") & Compound %in% c("Glucose", "C16:0")) %>%
  filter(BW <= 60 & BW >= 40) %>% 
  mutate(BW.range = ifelse(BW <= 50, "40-50", "50-60")) 

d.Fcirc.BW.matched %>% 
  ggplot(aes(x = BW.range, y = Fcirc_animal, fill = phenotype)) +
  stat_summary(geom = "bar", fun = "mean", position = "dodge", 
               alpha = .5, color = "black", width = .7) +
  stat_summary(
    geom = "errorbar", fun.data = "mean_sdl", fun.args = list(mult = 1), 
    position = position_dodge(.7), width = .5 ) +
  geom_quasirandom(dodge.width = .7, width = .1, show.legend = F) +
  scale_y_continuous(labels = function(x) x/ 1000,
                     expand = expansion(mult = c(0, .1)), 
                     n.breaks = 6) +
  facet_wrap(~Compound, scales = "free") +
  labs(x = "body weight (g)", y = "µmol molecules / min /animal\n") + 
  scale_fill_manual(values = color.phenotype) +
  scale_x_discrete(expand = expansion(mult = c(.6, .6))) + 
  theme.myClassic 

ggsave(filename = "Fcirc matched body weight.pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 3, width = 7)

# calculate stats
d.Fcirc.BW.matched %>%
  nest(-c(Compound, BW.range)) %>% 
  mutate(model = map(data, ~t.test(Fcirc_animal ~ phenotype, data = .x, ))) %>% 
  mutate(glance = map(model, broom::glance)) %>%  
  unnest(glance) %>% 
  mutate(stars = map_chr(p.value, ~func.generate_stars(.x)), .before = 3)



# Generalized linear model, with interaction term of BW and phenotype
# none of the interaction terms are significant
model.glm <- glm(formula = Fcirc_animal ~ BW + phenotype + BW:phenotype, 
                 family = gaussian(link = "identity"),
                 data = d.Fcirc.BW %>% filter(Compound == "Valine"))
model.glm %>% summary()


# traditional ANCOVA without interaction term
model.glm <- glm(formula = Fcirc_animal ~ BW + phenotype, 
                 family = gaussian(link = "identity"),
                 data = d.Fcirc.BW %>% filter(Compound == "Glutamine"))
summary(model.glm)

# perform all GLM
d.models <- d.Fcirc.BW %>% 
  # group_by(Compound) %>% 
  nest(-Compound) %>% 
  mutate(
    # ANCOVA, including BW
    ANCOVA = map(data, ~glm(formula = Fcirc_animal ~ BW + phenotype, 
                            family = gaussian(link = "identity"), data = .)),
    # anova, without BW
    anova = map(data, ~glm(formula = Fcirc_animal ~ phenotype, 
                           family = gaussian(link = "identity"), data = .))) %>% 
  mutate(tidied.ANCOVA = map(ANCOVA, tidy),
         tidied.anova = map(anova, tidy))

# ANCOVA unnested
d.ANCOVA <- d.models %>% unnest(tidied.ANCOVA) %>% 
  rename(p.ANCOVA = p.value) %>% 
  select(Compound, term, p.ANCOVA) 

# anova unnested
d.anova <- d.models %>% unnest(tidied.anova) %>% 
  rename(p.anova = p.value) %>% 
  select(Compound, term, p.anova) 

# combine the unnested dataset (note the terms are more in GLM, due to additional BW term)
BW.signif.compounds <- (d.ANCOVA %>% filter(term == "BW") %>% filter(p.ANCOVA < 0.05))$Compound

d.GLM <- d.ANCOVA %>% left_join(d.anova) %>% 
  mutate(p = ifelse(Compound %in% BW.signif.compounds, p.ANCOVA, p.anova)) %>% 
  filter(term != "(Intercept)")

# add stars
func.addStars <- function(x) {
  sapply(x, function(value) {
    if (is.na(value)) return("") else
      if (value <= 0.0001) return("****") else
        if (value <= 0.001) return("***") else
          if (value <= 0.01) return("**") else
            if (value <= 0.05) return("*") else
              return("")
  })
}

# add stars, prepare for visualization
d.GLM <- d.GLM %>% mutate(stars = func.addStars(p)) %>% 
  mutate(Compound = factor(Compound, levels = rev(ordered.Compound))) %>% 
  mutate(term = str_remove(term, "phenotype"))


plt.Fcirc.BW.pvalues <- d.GLM %>% 
  ggplot(aes(x = term, y = Compound, fill = stars)) +
  geom_tile(color = "snow4") +
  scale_fill_manual(values = c("grey100", "grey80", "grey40", "grey20")) +
  labs (y = NULL) +
  coord_cartesian(expand = 0) +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 15))
plt.Fcirc.BW.pvalues

plot_grid(plt.Fcirc.BW,
          ggplot() + theme_void(),
          plt.Fcirc.BW.pvalues, nrow = 1, 
          rel_widths = c(3.9, .1, 1))

ggsave(filename = "Fcirc_BW.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 6, width = 14)



# -<># -<># -<># -<># -<># -<># -<># -<># -<># -<># -<># -<># -<># -<># -<># -<># -<># -<># -<>
# Check nutrients level in serum

# analysis of Fcirc with TIC 
# define function to normalize a vector to [0, 1]
func.norm = function(x){
  x / max(x)
}

d.corrected.tidy2 <- d.corrected.tidy %>% 
  filter(! sample %in% c("bt-tail-O101")) %>%  # peculiar 3HB outlier that is 2-times high than the 2nd largest value
  filter(phenotype != "db/db") %>% 
  filter(C_Label == 0 & blood == "tail") %>% 
  # phenotype and compunds in order
  mutate(phenotype = factor(phenotype, levels = ordered.phenotype),
         Compound = factor(Compound, levels = ordered.Compound))

# glycerol measured as derivatized G3P
d.glycerol <- d.corrected.tidy2 %>% 
  filter(Compound == "Glycerol" & str_detect(MS.run,or( "10_G3P_ar_bd", "7_G3P_a-n_q-z-aa-ae", "9_G3P_HILIC_am-ao")))
# other nutrients
d.otherNutrients <- d.corrected.tidy2 %>% 
  filter(Compound != "Glycerol" & (!str_detect(MS.run, "G3P"))) %>% 
  filter(MS.run %in% c("11_ar_bd_EarlyElution", "12_ar_bd_LateElution", "18_bt-bw", "20_ce", "21_cf"))



# d.corrected.tidy2 %>%
#   filter(Compound == "Glucose") %>%
#   ggplot(aes(x = phenotype, y = intensity, color = phenotype)) +
#   geom_quasirandom(width = .1, alpha = .4) +
#   facet_wrap(~MS.run, nrow = 2, scales = "free") +
#   theme(legend.position = "none") +
#   expand_limits(y = 0)


# combine glycerol and other nutrients
d.intensity.norm.eachMouse <- d.glycerol %>% bind_rows(d.otherNutrients)  %>% 
  # normalize based on each separate MS run
  group_by(MS.run, Compound) %>%
  mutate(intens.norm = func.norm(intensity)) %>%
  # take the average for individual mouse
  group_by(phenotype, mouse_ID, Compound) %>%
  summarise(intens.norm.mean = mean(intens.norm))

d.intensity.norm.eachMouse %>% filter(Compound == "3-HB") %>% arrange(desc(intens.norm.mean))

# plot
plt.serum.metabolites <-  d.intensity.norm.eachMouse %>% 
  ggplot(aes(x = phenotype, y = intens.norm.mean, fill = phenotype)) + 
  stat_summary(fun = mean, geom = "bar", color = "black", alpha = .5) +
  stat_summary(fun.data = mean_sdl, geom = "errorbar",
               fun.args = list(mult = 1), 
               width = .3, color = "black") +
  geom_quasirandom(width = .3, alpha = .5) +
  # stat_summary(fun = mean, geom = "crossbar", width = .5, color = "black") +
  facet_wrap(~Compound, nrow = 2, scales = "free_x") +
  expand_limits(y = 0) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), breaks = seq(0, 1.2, .2)) +
  scale_x_discrete(expand = expansion(add = .8)) +
  labs(x = NULL, y = "relative abundance\n") +
  theme.myClassic +
  theme(panel.spacing = unit(20, "pt"), 
        legend.position = "none",
        axis.text = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 16)) +
  scale_fill_manual(values = color.phenotype)
plt.serum.metabolites

# perform significant analysis on serum metabolites
library(rstatix)
d.intensity.norm.eachMouse

d.intensity.signif <- d.intensity.norm.eachMouse %>% 
  group_by(Compound) %>% 
  pairwise_t_test(intens.norm.mean ~ phenotype, p.adjust.method = "bonferroni") %>% 
  add_xy_position(x = "phenotype") %>% 
  filter(p < 0.05)

plt.serum.metabolites +
  # crossbar
  geom_segment(
    data = d.intensity.signif,
    aes(x = xmin, xend = xmax, y = y.position, yend = y.position), 
    inherit.aes = F) +
  # stars
  geom_text(
    data = d.intensity.signif, size = 6,
    aes(x = (xmin +  xmax)/2, y = y.position + 0.02, label = p.signif), 
    inherit.aes = F)

ggsave(filename = "serum metabolomics.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 6, width = 10)


# Check concentration relationship with Fcirc
# use concentration in that specific infusion experiment, instead of using pooled data
d.intensity.normalized <- d.corrected.tidy %>%
  filter(C_Label == 0, blood == "tail") %>%
  select(MS.run, infused.tracer, Compound, phenotype, mouse_ID, infusion_round, intensity) %>%
  group_by(MS.run, infused.tracer, Compound) %>%
  mutate(intens.normalized = func.norm(intensity)) %>%
  mutate(infusion_mouseID = paste0(infusion_round, "_", mouse_ID))

# combine Fcirc data with normalized intensity data
d.intensity.normalized3 <- d.intensity.normalized %>% 
  filter(Compound == infused.tracer) %>% ungroup() %>% 
  select(phenotype, infusion_mouseID, Compound, intens.normalized)

d.Fcirc.BW.intensity <- d.Fcirc.BW %>% 
  left_join(d.intensity.normalized3, by = c("Compound", "infusion_mouseID", "phenotype"))

# !!! NOTE that the normalization is NOT perfect. 
# e.g., for measurement of lactate, the HFD samples were ran alone, and normalized by itself, without WT and ob/ob
# and thus cann't compare their concentration directly
d.Fcirc.BW.intensity %>% 
  ggplot(aes(x = intens.normalized, y = Fcirc_animal, color = phenotype)) + 
  geom_point(size = 3, alpha = .6) +
  geom_smooth(method = "lm", se = F) +
  scale_y_continuous(labels = function(x){x/1000}, name = "µmol / min") +
  # scale_x_continuous(breaks = seq(0, 1, .2)) +
  facet_wrap(~Compound, scales = "free", nrow = 2) +
  scale_color_manual(values = color.phenotype) +
  theme.mybw +
  theme(axis.text.x = element_text(angle = 45),
        legend.position = "bottom")

ggsave(filename = "serum metabolomics_Fcirc.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 6, width = 10)



# remove infusion rounds containing heavy HFD mice (> 50g; 8-12 months on diet) with matched body weight as ob/ob
# these heavy HFD mice do NOT affect the whole-body level flux; 
# For flux normalized by body weight, we use ~ 45 g HFD as the typical body mass (3-4 months on fat diet)
d.normalized.tidy <- d.normalized.tidy %>%  filter(! infusion_round %in% rounds.matched.BW) 


save.image(file = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/5_core_labeling_analysis.RData")






