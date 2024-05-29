rm(list = ls())

# library(plyr)
library(rebus)
library(viridis)
library(lubridate)
library(readxl)
library(purrr)
library(broom)
library(RColorBrewer)
library(splines)
library(cowplot)
library(gridExtra)
library(ggrepel)
library(ggbeeswarm)
library(scales)
library(ggsci)
library(tidyverse)

load("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/2_core_13CO2_bolus_injection.RData")

options(dplyr.print_min = 4)

# --------------------------------------------------------------------

# Infusion data
d.infusion.normalize = read_excel(pathID, sheet = "infusion.normalization")

# infusion data
d.infusion = d.all.treated %>%
  filter(inj == "infusion" & treatment == "13C") %>% 
  filter(! round %in% c(43, 61 ))

# combine with normalization dataset: normalize to 0.1 umol 13C atoms infused / min /animal
d.infusion = d.infusion %>%
  left_join(d.infusion.normalize %>% select(cage, infusion.norm.factor), by = "cage") 

# normalize to infusion of 0.1 umol 13C atoms /min/animal; overwrite the normalized column (prior for bolus normalization)
d.infusion = d.infusion %>% 
  mutate(P1.umol.min.normalized = P1.umol.min / infusion.norm.factor) %>% 
  filter(is.finite(P1.umol.min.normalized))


# all time.h subtract ca 5 min dead time due to the ca. 10 uL dead volume with 2 uL/min infusion rate; this is needed as the exponential model is forced to pass (0, 0)
d.infusion = d.infusion %>% 
  mutate(time.h = time.h - 5/60,
         infusion.length = infusion.length - 5/60) %>%  # infusion length also adjusted as it determines the pump turnoff time (marking the start of decay phase)
  filter(P1.umol.min.normalized > -.02) # remove those sudden abnormal downward signal spikes

# remove outliers
d.infusion = d.infusion %>% 
  filter(! (P1.umol.min.normalized < .02 & tracer == "Glutamine" & time.h > 1)) %>% # abnormal downward spikes in glutamine after 2 h
  filter(! cage %in% c(
    435, # L64, lactate infusion, low recovery 
    458, # O58, C16:0
    463,# O51, Alanine, ob/ob, much lower signal than others
    476, 471, # O53, O56, lactate infusion, no tracer signal
    490, 494, # L63, L64, 3-HB, abnormal signal 
    503, 508, 512, # O58, L62, L64, Acetate, no exponential curve
    526,# O61, 3HB infusion, no signal
    535, # O60, 3HB infusion, low signal  # tether # 3 has problem, leading to low signal to 490, 508, 526 and 535
    552, # HFD, glucose infusion, line twisted, no 13CO2 signal
    678, # HFD, H6,C16:0 infusion, no 13CO2 signal since the start of infusion, found dead at the end of infusion,
    658, # HFD, valine infusion,  13CO2 start to drop after 2 h
    991, 1036 # ob/ob, glucose infusion (6h group), too lower recovery
  ) ) %>% 
  filter(! (cage == 499 & time.h > 2.5) ) %>%  # signal abnormally drops, possibly tracer ran out in syringe
  filter(! (cage %in% c(658, 664, 662) & time.h > 2) )   # WT, HFD(H-20), valine infusion, 13CO2 signal start to drop slightly but continuously after 2 h; appeared very weak after infusion


# exponential - steady state curve
d.infusion.expCurve = d.infusion %>% filter(time.h*60 <= infusion.length)


func.plot.infusion.kinetics.modeling <- 
  function(mydata = d.infusion.expCurve) {
    
    # plot kinetic curve
    plt.infusion.expCurve.kinetics = mydata %>% 
      # filter(round == 109) %>% 
      filter(time.h >= 0) %>% 
      ggplot(aes(x = time.h, y = P1.umol.min.normalized, color = phenotype )) +
      geom_point(alpha = .4, size = .8) +
      geom_line(aes(group = cage), alpha = .3, size = .2) +
      facet_wrap(~tracer, nrow = 2, scales = "free_x") +
      labs(y = "13CO2 exhaled (µmol / min / animal)\n")  
    
    # fit exponential stage with a exponential model 
    d.infusion.expCurve.modeled.all = mydata %>% 
      filter(time.h >=0) %>% 
      nest(-c(tracer, phenotype)) %>% 
      # mutate(model = map(data, ~nls(P1.umol.min ~ a - b * exp(-c * time.h), start = list(a = .05, b = 10, c= 1), data = .))) %>% 
      # force to pass (0, 0)
      mutate(model = map(data, 
                         ~nls(P1.umol.min.normalized ~ a - a * exp(-c * time.h), 
                              start = list(a = .05,  c= 1), data = .))) %>% 
      mutate(fitted = map2(.x = data, .y = model,
                           .f = ~ tibble(time.h = .x$time.h, fitted = predict(.y)  ) ) ) %>% 
      mutate(tidied = map(model, ~tidy(.)))
    
    # fitted dataset
    d.infusion.expCurve.modeled = d.infusion.expCurve.modeled.all %>% unnest(fitted)
    
    # model summary
    d.model.summary = d.infusion.expCurve.modeled.all %>% unnest(tidied) %>% 
      mutate(error.pct = std.error / estimate * 100) %>% 
      select(-c(data, model, fitted)) %>% arrange(tracer) 
    d.model.summary
    
    
    # add exponential curve fitted to the kinetics  
    plt.kinetics.infusion.expCurve.modeled = 
      plt.infusion.expCurve.kinetics +
      geom_line(
        data = d.infusion.expCurve.modeled,
        aes(x = time.h, y = fitted, color = phenotype), size = 1.2) +
      expand_limits(y = 0) + 
      scale_y_continuous(expand = expansion(mult = c(0, .02)),
                         breaks = seq(0, .1, .02), 
                         position = "right",
                         sec.axis = sec_axis( breaks = seq(0, 1, .20), 
                                              trans = ~./.1, # * 100,
                                              name="recovery"),
      ) + 
      coord_cartesian(xlim = c(0, NA), ylim = c(0, .12)) +
      labs(title = "normalized to 0.1 µmol 13C atoms infused /min/animal") +
      theme_bw(base_size = 14) +
      theme(legend.title = element_blank(),
            legend.key.size = unit(.8, 'cm'),
            axis.line = element_line(linewidth = .5),
            legend.position = "bottom",
            panel.grid = element_blank(),
            strip.text = element_text(face = "bold"),
            strip.background = element_blank()
      ) +
      scale_color_manual(
        values = c("HFD" = "chocolate1", "ob/ob" =  "deepskyblue3", "WT" = "grey30"),
        limits = c("WT", "HFD", "ob/ob"),
        labels = function(x){str_replace(x, "WT", "control")}) +
      # legend line size 
      guides(color = guide_legend(override.aes = list(linewidth = 3), 
                                  nrow = 1)) 
    
    myoutput <- list(
      plot = plt.kinetics.infusion.expCurve.modeled, 
      model.summary = d.model.summary)
    
    return(myoutput)
  }

# all 3 phenotypes
plt.kinetics.infusion.expCurve.modeled_3phenotypes <- 
  func.plot.infusion.kinetics.modeling(
    mydata = d.infusion.expCurve %>% filter(! tracer %in% c("Acetate", "Methionine", "Tryptophan"))
  )[[1]]

plt.kinetics.infusion.expCurve.modeled_3phenotypes

ggsave(filename = "13CO2 recovery_infusion_3phenotypes.pdf", 
       plot = plt.kinetics.infusion.expCurve.modeled_3phenotypes, 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 6, width = 10)

# only WT
plt.kinetics.infusion.expCurve.modeled_WT <- 
  (func.plot.infusion.kinetics.modeling(
    mydata = d.infusion.expCurve %>% 
      filter(phenotype == "WT" & tracer != "Acetate")
  ))[[1]]

plt.kinetics.infusion.expCurve.modeled_WT +
  coord_cartesian(expand = 0, ylim = c(0, .13))

ggsave(filename = "13CO2 recovery_infusion_WT.pdf", 
       plot = plt.kinetics.infusion.expCurve.modeled_WT + 
         theme(text = element_text(size = 20),
               axis.title.y.left = element_text(margin = margin(r = 15)),
               axis.title.y.right = element_text(margin = margin(l = 5))) + labs(x = "time (h)"), 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 7, width = 13)


# plot infusion of glucose to WT
plt.kinetics.infusion.expCurve.modeled_WT_glucose <- 
  (func.plot.infusion.kinetics.modeling(
    mydata = d.infusion.expCurve %>% filter(phenotype == "WT" & tracer == "Glucose")
  ))[[1]] +
  # overwrite the default with x-limits of 200 min
  coord_cartesian(expand = 0) 

plt.kinetics.infusion.expCurve.modeled_WT_glucose

ggsave(filename = "13CO2 recovery_infusion_WT_palm.pdf", 
       plot = plt.kinetics.infusion.expCurve.modeled_WT_glucose, 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 4, width = 4)


# Calculate recovery
d.recovery.infusion = d.infusion.expCurve %>% 
  # take average over time points after 2.6 h 
  filter((cage !=561 & time.h >= 2.6) | 
           (cage == 561 & time.h >= 2.6)) %>% # cage 561, C18:1, HFD,  started infusion late at 4:11 pm, and less than 3 h infusion was conducted
  # for glucose, calculate recovery after 4h
  filter(tracer == "Glucose" & time.h >= 4 |
           tracer != "Glucose" & time.h > 0) %>% 
  group_by(phenotype, tracer, inj.umol.13C.atoms, cage, round) %>% 
  dplyr::summarise(P1.umol.min.mean = mean(P1.umol.min),
                   P1.umol.min.SD = sd(P1.umol.min)) %>% 
  dplyr::mutate(recovery.infuse = P1.umol.min.mean / inj.umol.13C.atoms,
                recovery.infuse.SD = P1.umol.min.SD / inj.umol.13C.atoms)


# plot recovery of infusion
plt.recovery.infusion = d.recovery.infusion %>% 
  mutate(phenotype = factor(phenotype, levels = ordered.phenotype)) %>% 
  ggplot(aes(x = phenotype, y = recovery.infuse, color = phenotype)) +
  stat_summary(geom = "crossbar", fun = "mean", size = 1, show.legend = F) +
  geom_quasirandom(size = 2, alpha = .7) +
  # geom_text(aes(label = cage), color = "black") +
  # geom_errorbar(aes(ymax = recovery.infuse + recovery.infuse.SD,
  #                   ymin = recovery.infuse - recovery.infuse.SD)) +
  facet_wrap(~tracer, nrow = 2, scales = "free_x")+
  scale_y_continuous(expand = expansion(mult = c(0, .1)),
                     n.breaks = 7,
                     name = "infusion recovery") +
  expand_limits(y = 0) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.length.x = unit(-1, "mm")) +
  # geom_text_repel(aes(label = cage)) +
  # geom_label(aes(label = recovery.infuse %>% round(2)), size = 2, 
  #            position = position_jitter(.3, 0)) +
  scale_color_manual(values = color.phenotype,
                     labels = function(x){str_replace(x, "WT", "control")})

plt.recovery.infusion 



# -<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-


# Combine recovery of infusion and bolus 
d.recovery.infuse.bolus = d.recovery.infusion %>% ungroup() %>% 
  select(phenotype,tracer, recovery.infuse)  %>% rename(recovery = recovery.infuse) %>% 
  mutate(inj = "infusion") %>% 
  # combine with bolus recovery
  rbind(CO2.13.bolus.recovered %>% ungroup() %>% 
          select(phenotype,tracer, recovery, inj)) %>% 
  # phenotype in ordered factor
  mutate(phenotype = factor(phenotype, levels = ordered.phenotype, ordered = T)) %>% 
  mutate(inj = str_replace(inj, "tail.vein_anesthesia", "bolus inj."))

# remove outliers: 
d.recovery.infuse.bolus <- d.recovery.infuse.bolus %>% 
  filter(! (phenotype == "ob/ob" & tracer == "Glycerol" & recovery > .7)) %>%  # remove a single outlier
  filter(! (phenotype == "HFD" & tracer == "Alanine" & recovery < .5)) %>%  # remove a single outlier from HFD alanine, too low
  filter(! (phenotype == "WT" & tracer == "Glucose" & recovery > .8)) # remove a single outlier from WT glucose way too high

# summary stats of each of the method
d.recovery.infuse.bolus.summary <- d.recovery.infuse.bolus %>% 
  group_by(phenotype, tracer, inj) %>% 
  summarise(recovery.mean = mean(recovery, na.rm = T),
            recovery.SD = sd(recovery, na.rm = T),
            n.rep = length(tracer),
            recovery.SEM = recovery.SD / sqrt(n.rep))


# summary stats of pooled data from both methods
d.recovery.infuse.bolus.pooled.summary <- d.recovery.infuse.bolus %>% 
  # remove an ob/ob - glycerol outlier that has aberrant too high recovery
  filter(!(phenotype == "ob/ob" & tracer == "Glycerol" & recovery > 0.798)) %>% 
  group_by(phenotype, tracer) %>% 
  summarise(recovery.mean = mean(recovery, na.rm = T),
            recovery.SD = sd(recovery, na.rm = T),
            n.rep = length(tracer),
            recovery.SEM = recovery.SD / sqrt(n.rep)) %>% 
  mutate(method = "pooled") %>% 
  # adjust column order to be aliged with "d.recovery.infuse.bolus.summary"
  select(phenotype, tracer, method, recovery.mean, recovery.SD, n.rep, recovery.SEM)

# update - for glycerol, only use infusion data
d.recovery.infuse.bolus.pooled.summary <- 
  d.recovery.infuse.bolus.pooled.summary %>% 
  filter(tracer != "Glycerol") %>% 
  # combine with glycerol infusion method  
  rbind(d.recovery.infuse.bolus.summary %>% 
          filter(tracer == "Glycerol" & inj == "infusion") %>% rename(method = inj)
  )

# Output recovery
library(xlsx)

# sheet 1: recovery from bolus injection and infusion, respectively
write.xlsx(d.recovery.infuse.bolus.summary %>% as.data.frame(), 
           file = "nutrients 13CO2 recovery and fox.xlsx", sheetName = "recovery_each_method")
# sheet 2: recovery from pooled data of bolus injection and infusion
write.xlsx(d.recovery.infuse.bolus.pooled.summary %>% as.data.frame(), 
           file = "nutrients 13CO2 recovery and fox.xlsx", sheetName = "recovery_pooled", append = T)
# note that in the pooled dataset, glycerol is using infusion data




# Calculate fraction of oxidation fox after correction by bicarbonate recovery
x <- (d.recovery.infuse.bolus %>% filter(tracer == "Bicarbonate"))
recovery.bicarbonate <- x$recovery %>% mean()

# fox of bolus & infusion, separately 
d.fox.infuse.bolus.summary <- d.recovery.infuse.bolus.summary %>% 
  mutate(across(contains("recovery"), .fns = ~.x / recovery.bicarbonate))
# update names from recovery into fox
names(d.fox.infuse.bolus.summary) <- 
  names(d.fox.infuse.bolus.summary) %>% str_replace("recovery", "fox")


# fox of pooled bolus & infusion
d.fox.infuse.bolus.pooled.summary <- d.recovery.infuse.bolus.pooled.summary %>%
  mutate(across(.col = c(recovery.mean, recovery.SD, recovery.SEM),
                .fns = ~.x / recovery.bicarbonate)) 
# update names from recovery into fox
names(d.fox.infuse.bolus.pooled.summary) <- 
  names(d.fox.infuse.bolus.pooled.summary) %>% str_replace("recovery", "fox")




# Output fox data into excel
# sheet 3: fox from bolus injection and infusion, separately
write.xlsx(d.fox.infuse.bolus.summary %>% as.data.frame(), 
           file = "nutrients 13CO2 recovery and fox.xlsx", sheetName = "fox_each_method", append = T)

# sheet 4: fox of pooled data
write.xlsx(d.fox.infuse.bolus.pooled.summary %>% as.data.frame(), 
           file = "nutrients 13CO2 recovery and fox.xlsx", sheetName = "fox_pooled", append = T)
# note that in the pooled dataset, glycerol is using infusion data


# Visualization of recovery

# only WT
dg <- .7
plt.recovery.WT <-  d.recovery.infuse.bolus %>% 
  mutate(tracer = fct_reorder(tracer, recovery, .desc = T)) %>% 
  filter(phenotype == "WT") %>% 
  ggplot(aes(x = tracer, y = recovery, fill = inj, color = inj)) +
  stat_summary(geom = "bar", fun = mean, alpha = .4, 
               position = position_dodge(dg),
               width = dg) +
  stat_summary(geom = "errorbar", fun.data = mean_sdl, 
               fun.args = list(mult = 1),
               width = .4,
               position = position_dodge(dg),
               show.legend = F) +
  geom_quasirandom(dodge.width = dg, 
                   show.legend = F,
                   size = 1) +
  theme.myClassic +
  scale_y_continuous(expand = expansion(mult = c(0, .1)),
                     n.breaks = 6) +
  scale_x_discrete(expand = expansion(add = .8)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(25, "pt")) +
  labs(x = NULL, y = "recovery\n") +
  theme(legend.position = c(.9, .9)) +
  # label the average recovery
  geom_text(data = d.recovery.infuse.bolus.summary %>% filter(phenotype == "WT"),
            aes(x = tracer, y = recovery.mean + .15, 
                label = round(recovery.mean, 2),
                color = inj),
            angle = 90, position = position_dodge(dg), fontface = "bold",
            show.legend = F)

plt.recovery.WT

ggsave(filename = "recovery_separate_WT.pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 4, width = 7)


# ---- style 2 -----averaged bars --------
k <- d.recovery.infuse.bolus %>% 
  filter(phenotype == "WT") %>% 
  filter(! (tracer == "Glycerol" & inj == "bolus inj.")) %>% 
  mutate(tracer = fct_reorder(tracer, recovery, .fun = mean, .desc = F, .na_rm = T))


k %>% 
  filter(!(recovery > .8 & tracer == "Glucose")) %>% 
  # filter(tracer != "Acetate") %>% 
  ggplot(aes(x = tracer, y = recovery, label = recovery)) +
  stat_summary(geom = "bar", fun = mean, alpha = .3, 
               color = "black", fill = "coral", width = .8) +
  # errorbar
  stat_summary(geom = "errorbar", fun.data = mean_sdl, 
               fun.args = list(mult = 1),
               width = .4,
               position = position_dodge(dg),
               show.legend = F) +
  # point
  geom_quasirandom(dodge.width = dg, show.legend = F,
                   size = 1.5, shape = 21, alpha = 1, fill = "grey40") +
  # text labels  
  geom_text(data = k %>% 
              # filter(tracer != "Acetate") %>% 
              group_by(tracer) %>% summarise(recovery.mean = mean(recovery)),
            aes(x = tracer, 
                # y = .1, 
                y = recovery.mean + .2, 
                label = recovery.mean %>% round(2) %>% as.character()), #  %>% str_remove(START %R% "0")),
            inherit.aes = F, size = 6, color = "turquoise4", fontface = "bold") +
  theme.myClassic +
  scale_y_continuous(expand = expansion(mult = c(0, .1)),
                     breaks = seq(0, 1.2, .2)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(axis.text.x = element_text(hjust = 1),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(25, "pt"),
        legend.position = "none") +
  labs(x = NULL, y = "\nrecovery (pooled) of WT mice") +
  coord_flip() 

ggsave(filename = "recovery_mean_WT.pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 7, width = 4.1)



# ---- style 3 -----scatterplot --------
a <- d.recovery.infuse.bolus %>% 
  group_by(phenotype, tracer, inj) %>% 
  summarise(recovery.mean = mean(recovery),
            recovery.sd = sd(recovery, na.rm = T)) 

a1 <- a  %>% filter(inj == "bolus inj.") %>% 
  rename(bolus = recovery.mean, bolus.SD = recovery.sd) %>% select(-inj)

a2 <- a  %>% filter(inj == "infusion") %>% 
  rename(infusion = recovery.mean, infusion.SD = recovery.sd) %>% select(-inj)

A <- left_join(a1, a2, by = c("phenotype", "tracer"))

func.plt.recovery.scatter.2methods <- function(
    mydata, whichColor = "tracer", textsize = 3){
  mydata %>% 
    ggplot(aes_string(x = "bolus", y = "infusion", color = whichColor)) +
    # error - X
    geom_errorbar(aes(xmin = bolus - bolus.SD,
                      xmax = bolus + bolus.SD)) +
    # error - Y
    geom_errorbar(aes(ymin = infusion - infusion.SD,
                      ymax = infusion + infusion.SD)) +
    # point
    geom_point(shape = "diamond", size = 5) +
    # text
    ggrepel::geom_text_repel(aes(label = tracer), max.overlaps = Inf, size = textsize) +
    # anti-diagonal line
    geom_abline(slope = 1, intercept = 0) +
    geom_abline(slope = .9, intercept = 0, linetype = "dashed", color = "snow4") +
    geom_abline(slope = 1.1, intercept = 0, linetype = "dashed", color = "snow4") +
    # theme
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1),
                expand = 0, clip = "off") +
    scale_x_continuous(breaks = seq(0, 1, .2)) +
    scale_y_continuous(breaks = seq(0, 1, .2)) +
    scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(11)) +
    theme.myClassic +
    theme(legend.position = "none",
          plot.margin = margin(rep(10, 4), unit = "pt")) +
    labs(x = "bolus recovery", y = "infusion recovery")
}

# WT only
A %>% filter(phenotype == "WT") %>% 
  func.plt.recovery.scatter.2methods(textsize = 4) +
  
  # annotate bolus values
  geom_segment(aes(x = bolus, y = infusion, xend = bolus, yend = 0),
               linetype = "dotted") +
  geom_text_repel(aes(label = round(bolus, 2), y = 0), 
                  angle = 0, max.overlaps = Inf, fontface = "bold") +
  
  # annotate infusion values
  geom_segment(aes(x = bolus, y = infusion, xend = 0, yend = infusion), 
               linetype = "dotted") +
  geom_text_repel(aes(label = round(infusion, 2), x = 0),
                  angle = 0, max.overlaps = Inf, fontface = "bold") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15))


ggsave(filename = "recovery 2 methods WT.pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 5, width = 5)

# 3 phenotypes
A %>% func.plt.recovery.scatter.2methods(whichColor = "phenotype", textsize = 3) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "right") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18))

ggsave(filename = "recovery 2 methods 3pheno.pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 7, width = 7)
# ----


# plot all 3 phenotypes, both bolus & infusion methods

plt.recovery.bolus.infusion_3pheno =  
  d.recovery.infuse.bolus %>% 
  mutate(tracer = fct_reorder(tracer, recovery, .desc = T)) %>% 
  ggplot(aes(x = phenotype, y = recovery, color = phenotype, fill = phenotype,
             shape = inj)) + 
  # average bars
  stat_summary(aes(group = phenotype), 
               fun = mean, 
               fun.args = list(mult = 1),
               geom = "bar", 
               width = 1, color = "black", alpha = .7) +
  # error bars
  stat_summary(aes(group = phenotype), 
               fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar",
               width = .5, color = "black") +
  facet_wrap(~tracer, scales = "free_x", nrow = 2) +
  # points
  # bolus
  geom_quasirandom(
    data = d.recovery.infuse.bolus %>% filter(inj == "bolus inj."),
    dodge.width = .9, aes(shape = inj), 
    color = "black", size = 2, show.legend = F) +  
  # infusion
  # not including WT-glucose
  geom_quasirandom(
    data = d.recovery.infuse.bolus %>% filter(inj == "infusion" & (!(tracer == "Glucose" & phenotype == "WT"))),
    dodge.width = .9, aes(shape = inj),
    fill = "white", color = "black", size = 2) +
  
  # only WT-glucose
  geom_quasirandom(
    data = d.recovery.infuse.bolus %>% 
      filter(inj == "infusion" & tracer == "Glucose" & phenotype == "WT"),
    dodge.width = .9, aes(shape = inj), width = .5,
    fill = "white", color = "black", size = 1.6) +
  
  # theme
  theme.myClassic + 
  theme(axis.ticks.length.x = unit(-1, "mm"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())  +
  expand_limits(y = 0) + 
  scale_y_continuous(breaks = seq(0, 1, .2), 
                     expand = expansion(mult = c(0, .2))) +
  scale_x_discrete(expand = expansion(mult = c(.5, .5))) +
  scale_color_manual(values = c("grey30" ,   "chocolate1" , "deepskyblue3"))  +
  scale_fill_manual(values = c("grey30" , "chocolate1" , "deepskyblue4"))  +
  labs(y = "13CO2 recovery fraction\n") +
  scale_shape_manual(values = c(21, 23)) 
# guides(shape = guide_legend(override.aes = list(
#   size = 4, stroke = 1, fill = c("black", "white"))))

plt.recovery.bolus.infusion_3pheno


ggsave(filename = "13CO2_recovery_bolus_infuse_3phenotypes.pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 5, width = 12)


# plot recovery differences between the two methods
d.recovery.bolus.infusion = d.recovery.infuse.bolus %>% 
  group_by(inj, tracer, phenotype) %>% 
  summarise(recovery.mean = mean(recovery, na.rm = T)) %>% 
  spread(inj, recovery.mean) %>% 
  mutate(error.percent = (infusion - `bolus inj.`)/infusion * 100 )

plt.recovery.bolus.vs.infusion = d.recovery.bolus.infusion[complete.cases(d.recovery.bolus.infusion), ] %>%
  filter(phenotype %in% c("ob/ob", "WT", "HFD")) %>% 
  ggplot(aes(x = tracer, y  = error.percent, fill = phenotype)) + 
  geom_bar(stat = "identity", position = "dodge", width = .5) +
  coord_flip() +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = color.phenotype)
plt.recovery.bolus.vs.infusion



# plot the fox data -<>--<>--<>--<>--<>--<>-

d.fox.infuse.bolus <- d.recovery.infuse.bolus %>% 
  mutate(fox = recovery / recovery.bicarbonate)

f <- d.fox.infuse.bolus %>% 
  filter(phenotype == "WT") %>% 
  filter(tracer != "Bicarbonate") %>% 
  filter(! (tracer == "Glycerol" & inj == "bolus inj.")) %>% 
  mutate(tracer = fct_reorder(tracer, fox, .fun = mean, .desc = T, .na_rm = T))

f %>% 
  ggplot(aes(x = tracer, y = fox, label = recovery)) +
  stat_summary(geom = "bar", fun = mean, alpha = .3, 
               color = "black", fill = "coral", width = .8) +
  # errorbar
  stat_summary(geom = "errorbar", fun.data = mean_sdl, 
               fun.args = list(mult = 1),
               width = .4,
               position = position_dodge(dg),
               show.legend = F) +
  # point
  geom_quasirandom(dodge.width = dg, show.legend = F,
                   size = 1.5, shape = 21, alpha = 1, fill = "grey40") +
  # text labels  
  geom_text(data = f %>% 
              group_by(tracer) %>% summarise(fox.mean = mean(fox)),
            aes(x = tracer, 
                y = fox.mean + .2, 
                label = fox.mean %>% round(2) %>% as.character(),
                angle = 60), 
            inherit.aes = F, size = 6, color = "turquoise4", fontface = "bold") +
  theme.myClassic +
  scale_y_continuous(expand = expansion(mult = c(0, .1)),
                     breaks = seq(0, 1.2, .2)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        axis.text = element_text(size = 16),
        axis.text.y = element_text(size = 17),
        legend.position = "none") +
  labs(x = NULL, y = "fraction of oxidation (fox)\n") 


ggsave(filename = "fox_mean_WT.pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 4.1, width = 7.2)



# 1) MERGE the data from the bolus and infusion methods
# 2) but only use infusion data for glycerol
# 3) Show fraction of oxidation (fox)
# 4) remove bicarbonate
# 5) remove acetate, methionie, and tryptophan
dg <- .7

d.fox.infuse.bolus.subset <- 
  d.fox.infuse.bolus %>% 
  filter(! tracer %in% c("Bicarbonate", "Acetate", "Methionine", "Tryptophan")) %>% 
  filter(! (tracer == "Glycerol" & inj == "bolus inj.")) %>% 
  mutate(tracer = fct_reorder(tracer, recovery, .desc = T)) 
# mutate(tracer = factor(tracer, levels = c(
#   "Lactate", "Alanine", "Glucose", "3-HB", "Glutamine", "Glycerol", "C16:0","C18:1", "C18:2", "Valine"))) 

plt.fox.3phenotypes <- d.fox.infuse.bolus.subset %>% 
  ggplot(aes(x = tracer, y = fox, color = phenotype, fill = phenotype)) + 
  # average bars
  stat_summary(aes(group = phenotype), fun = mean, 
               color = "black", width = dg,
               fun.args = list(mult = 1), position = position_dodge(dg),
               geom = "bar", width = dg, color = "black", alpha = .7) +
  # error bars
  stat_summary(aes(group = phenotype), 
               fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", position = position_dodge(dg),
               width = .3, color = "black") +
  # points
  geom_quasirandom(shape = 21, color = "black", alpha = .6, width = .1,
                   dodge.width = dg, size = 2, show.legend = F) +
  expand_limits(y = 0) +
  scale_y_continuous(breaks = seq(0, 1, .2), 
                     expand = expansion(mult = c(0, .1))) +
  scale_x_discrete(expand = expansion(mult = c(.08, .05))) +
  scale_fill_manual(values = c("grey40" , "chocolate1" , "deepskyblue3"),
                    labels = function(x){str_replace(x, "WT", "control")})  +
  labs(y = "13CO2 fox\n") +
  scale_shape_manual(values = c(21, 23)) +
  coord_cartesian(clip = "off") +
  guides(shape = guide_legend(override.aes = list(
    size = 4, stroke = 1, fill = c("black", "white")))) +
  # theme
  theme.myClassic + 
  theme(
    # legend.position = c(1, .85),
    axis.ticks.length.x = unit(-1, "mm"),
    axis.text.x = element_text(angle = 40, hjust = 1),
    axis.title.x = element_blank(),
    plot.margin = margin(r = 34, unit = "pt"))  

plt.fox.3phenotypes

ggsave(filename = "13CO2_fox_bolus_infuse_3phenotypes.pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 2.5, width = 12)



# add label
plt.fox.3phenotypes +
  geom_text(data = d.fox.infuse.bolus.subset %>%
              group_by(phenotype, tracer) %>% 
              summarise(fox.mean = mean(fox, na.rm = T)),
            aes(label = round(fox.mean, 2), y = fox.mean + .14, color = phenotype),
            angle = 60, hjust = 0, fontface = "bold",
            position = position_dodge(.9)) +
  scale_color_manual(values = color.phenotype)

ggsave(filename = "13CO2_fox_bolus_infuse_3phenotypes_numbered.pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 3.7, width = 12) 


# Check number of replicate
d.fox.infuse.bolus %>% group_by(tracer, phenotype, inj) %>% 
  summarise(n = n()) %>% spread(tracer, n) %>% 
  arrange(inj, phenotype)


# signifance test
library(rstatix)

d.fox.infuse.bolus.pooled.summary %>% 
  select(phenotype, tracer, fox.mean) %>% 
  spread(phenotype, fox.mean) %>% 
  mutate(ratio.WT_HFD = WT / HFD - 1,
         ratio.WT_ob = WT / `ob/ob` - 1,
         ratio.HFD_ob = HFD / `ob/ob` - 1) %>% 
  filter(abs(ratio.WT_HFD) > .1 | abs(ratio.WT_ob) > .1 | ratio.HFD_ob > .1 )


p <- d.fox.infuse.bolus.subset %>% 
  group_by(tracer) %>% 
  pairwise_t_test(fox ~ phenotype, p.adjust.method = "bonferroni") %>% 
  add_xy_position(x = "tracer") %>% 
  arrange(p)

crosswidth = .3

p2 <- p %>% 
  mutate(xmin = ifelse(xmin < x,  x - crosswidth, xmin),
         xmax = ifelse(xmax > x,  x + crosswidth, xmax),
         # lower the height of all bars
         y.position = ifelse(tracer != "Glutamine", y.position - .15, y.position)) %>% 
  filter(p.adj < 0.05) %>% 
  filter(! (tracer == "C18:2" & group1 == "WT" & group2 == "HFD")) %>% 
  group_by(tracer) %>% 
  # switch the label position
  mutate(y.position = ifelse(y.position == max(y.position), min(y.position), max(y.position))) %>% 
  # adjust bar height of glutamine
  mutate(y.position = ifelse(tracer == "Glutamine" & group1 == "WT", 1, y.position)) %>% 
  mutate(y.position = ifelse(tracer == "Glutamine" & group1 == "HFD", 1.072, y.position))

# add significant stars
plt.fox.3phenotypes.stars <- plt.fox.3phenotypes +
  # crossbar
  geom_segment(
    data = p2,
    aes(x = xmin - 4, xend = xmax - 4, y = y.position, yend = y.position), 
    inherit.aes = F) +
  # stars
  geom_text(
    data = p2, size = 6,
    aes(x = (xmin +  xmax)/2 - 4, y = y.position + 0.02, label = p.adj.signif), 
    inherit.aes = F)

plt.fox.3phenotypes.stars

ggsave(filename = "13CO2_fox_bolus_infuse_3phenotypes.pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 3.7, width = 12)


# style II
plt.fox.3phenotypes +
  facet_wrap(~tracer, scales = "free_x", nrow = 2) +
  theme(legend.position = "bottom") +
  scale_x_discrete(expand = expansion(mult = c(.5, .5))) +
  theme(axis.text.x = element_blank())

ggsave(filename = "13CO2_fox_bolus_infuse_3phenotypes_faceted.pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 5.5 * 1.1, width = 9 * 1.1)



# # Pool the recovery from both I.V. bolus injection and infusion
# # (for glycerol, only use the infusion data for network computation)
# d.recovery.infuse.bolus.pooled.summary = d.recovery.infuse.bolus %>% 
#   group_by(phenotype, tracer) %>% 
#   summarise(recovery.mean = mean(recovery, na.rm = T),
#             recovery.SD = sd(recovery, na.rm = T),
#             n.rep = n_distinct(recovery),
#             recovery.SEM = recovery.SD / sqrt(n.rep) ) %>% 
#   mutate(method = "pooled")
# d.recovery.infuse.bolus.pooled.summary
# 
# d.recovery.infuse.bolus$phenotype %>% unique()

d.fox.infuse.bolus %>% 
  filter(phenotype %in% c("WT", "HFD", "ob/ob")) %>% 
  filter(! tracer %in% c("Methionine", "Tryptophan")) %>% 
  group_by(tracer) %>% 
  rstatix::pairwise_t_test(
    paired = F, fox ~ phenotype , p.adjust.method = "bonferroni"
  ) %>% 
  filter(p.adj < .05) %>% 
  select(-c(p, p.signif)) 


# ------------------------------------------------------------



# Compare body weight used in experiments of 13CO2 measurement
# We'll use the body weight from the serum labeling infusion experiments
# that data is directly involved in flux calculation
d.all %>% select(phenotype, cage, tracer, BW) %>% 
  filter(phenotype %in% c("WT", "ob/ob", "HFD")) %>% distinct() %>% 
  mutate(phenotype = factor(phenotype, levels = ordered.phenotype))  %>% 
  ggplot(aes(x = phenotype, y = BW, fill = phenotype)) +
  stat_summary(geom = "bar", fun = mean, alpha = .7, color = "black") +
  geom_quasirandom(width = .2, size = 1, show.legend = F)+
  # facet_wrap(~tracer, nrow = 2, scales = "free_x") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)),
                     breaks = seq(0, 80, 10),
                     name = "BW (in 13CO2 experiment)\n") +
  scale_x_discrete(expand = expansion(add = .8),
                   labels = function(x){str_replace(x, "WT", "control")}) +
  expand_limits(y = 0) +
  scale_fill_manual(values = color.phenotype,
                    labels = function(x){str_replace(x, "WT", "control")}) +
  theme.myClassic +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_blank()) 



#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>
#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>#<>



# Validation: radioactive 14C carcass signal +13CO2 recovery = 100%, confirming accurate measurement
d.14C = read_excel(pathID, sheet = "14C")
C14.rec <- (d.14C$carc.recovery)[c(1:3)]
C14.recovery.mean <- mean(C14.rec)
C14.recovery.sd <- sd(C14.rec)

# WT glucose bolus-based recovery
kkk <- d.recovery.infuse.bolus %>% filter(inj == "bolus inj." & phenotype == "WT" & tracer == "Glucose") %>% 
  group_by(phenotype) %>% 
  summarise(recovery.mean = mean(recovery),
            recovery.SD = sd(recovery))

plt.14C <- tibble(phenotype = "WT", 
                  recovery.mean = C14.recovery.mean, 
                  recovery.SD = C14.recovery.sd) %>% 
  rbind(kkk) %>% 
  mutate(treatment = factor(c( "14C-carcass", "13CO2"), ordered = T)) %>% 
  mutate(error.y = c( C14.recovery.mean + kkk$recovery.mean, 
                      kkk$recovery.mean)) %>% 
  ggplot(aes(x = 1, y = recovery.mean, fill = treatment)) +
  geom_col(position = position_stack(reverse = T), 
           alpha = .6, color = "black") +
  geom_errorbar(aes(ymin = error.y - recovery.SD ,
                    ymax = error.y + recovery.SD),
                width  = .3) +
  theme.myClassic +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.y = element_text(face = "bold"),
        text = element_text(size = 19, color = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 23),
        plot.margin = margin(b = 10, unit = "pt")) +
  scale_y_continuous(breaks = seq(0, 1, .2),
                     expand = expansion(mult = c(0, .1))) +
  scale_x_discrete(expand = expansion(add = .3)) +
  # geom_hline(yintercept = 1, linetype = "dashed")  +
  labs(y = "recovery of injected dose\n") +
  scale_fill_manual(values = c("red4", "orange"))

plt.14C

ggsave(filename = "14C.pdf", 
       plot = plt.14C, 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 4.5, width = 2.9)




##-<>-##-<>-##-<>-##-<>-##-<>-##-<>-##-<>-##-<>-##-<>-##-<>-##-<>-##-<>-

# Calculate total CO2 production, total O2 consumption, and RER, 
# based on the last hour of the infusion experiments 

d.infusion.stable.0 <- d.infusion %>% 
  # filter(time.h >= 1.5 & time.h <= 2.5)
  filter(time.h >= 2.5 & time.h <= 3) %>% 
  # select experiments done before 2013 Xmax
  filter(round <= 109 )

# Calcuate total CO2 output rate from animal
d.infusion.stable.1 <- d.infusion.stable.0 %>% 
  mutate(P_total.mL.min = (ppm.t.out - ppm.in) * 2000 / 10^6,
         P_total.umol.min = P_total.mL.min/1000 / 22.4 * 10^6) 


# calculate the baseline of O2 for each single round of experiment
d.O2.baseline <- d.all %>% 
  filter(treatment == "baseline") %>% 
  filter(time.h > 3) %>% 
  group_by(round) %>% 
  summarise(O2.baseline.mean = mean(O2_A))

# calculate O2 consumption rate
d.infusion.stable.2 <- d.infusion.stable.1 %>% 
  left_join(d.O2.baseline, by = "round") %>% 
  mutate(O2.consume.mL.min = (O2.baseline.mean - O2_A)/100 * 2000,
         O2.consume.umol.min = O2.consume.mL.min/1000/22.4 * 10^6) %>% 
  select(-c(O2.baseline.mean, O2_A, contains("mL"))) %>% 
  # phenotype in order
  mutate(phenotype = factor(phenotype, levels = ordered.phenotype))

# for ID, if NA, assign a random number
d.infusion.stable.2 <- d.infusion.stable.2 %>% 
  mutate(mouseID = ifelse(is.na(mouseID), sample(1:1000, 1) %>% as.character(), mouseID))


# summarize for each mouse in each experiment
d.infusion.stable.summary1 <- d.infusion.stable.2 %>% 
  group_by(mouseID, round.cage, round, phenotype) %>% 
  # group_by(mouseID, phenotype) %>% 
  summarise(CO2.umol.min = mean(P_total.umol.min),
            O2.umol.min = mean(O2.consume.umol.min),
            RER = CO2.umol.min / O2.umol.min,
            BW = mean(BW)) %>% 
  
  # remove obvious 1-3 outliers in each phenotype that has RER > 1.5
  filter(RER < 1.5) %>% 
  filter(RER < .95) %>% 
  # remove 3 outliers in HFD with very low CO2 output,
  filter(! (phenotype == "HFD" & CO2.umol.min < 50)) %>% 
  
  # remove 1 outliers in HFD with very high O2 consumption
  filter(! (phenotype == "HFD" & O2.umol.min > 125)) 

# # summarize 
# d.infusion.stable.summary1 %>% ungroup() %>% 
#   select(-c(mouseID, round.cage, round)) %>% 
#   group_by(phenotype) %>% 
#   summarise(CO2 = mean(CO2.umol.min))


# significant analysis
library(rstatix)

# define function making bar plots
func.plt.basicsPhysiology <- function(
    mydata = d.infusion.stable.summary1,
    whichY = "CO2.umol.min", mytitle = whichY){
  
  mydata %>% 
    ggplot(aes_string(x = "phenotype", y = whichY, fill = "phenotype")) + 
    # mean
    stat_summary(fun = mean, geom = "bar", alpha = .7, color = "black") +
    # standard deviation
    stat_summary(fun.data = mean_sdl, geom = "errorbar",
                 fun.args = list(mult = 1),
                 width = .3) +
    # points
    geom_quasirandom(show.legend = F ) +
    theme.myClassic +
    scale_fill_manual(values = color.phenotype,
                      labels = function(x){str_replace(x, "WT", "control")}) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 19),
          axis.title.y.left = element_text(margin = margin(r = 8, unit = "pt")),
          axis.title.y.right = element_text(margin = margin(l = 8, unit = "pt")),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 17),
          axis.ticks.x = element_blank(),
          legend.position = "bottom") +
    scale_x_discrete(expand = expansion(add = 1)) +
    ggtitle(label = mytitle)
}

# define function processing significant analysis
func.sig1 <- function(
    dataset, # dataset output from pairwise_t_test()
    y.max.up = 5,  # upper crossbar shift further upward
    y.min.down = 5, # bottom crossbar further shifted downward
    xmin.max.shift = 0 # sfhit both xmin and xmax to align with bar position; negative to shift to the left
){
  dataset %>% 
    add_xy_position(x = "phenotype") %>% 
    filter(p.adj < 0.05) %>% arrange(xmax) %>% 
    mutate(y.position = ifelse(y.position == max(y.position), y.position + y.max.up, y.position)) %>%
    mutate(y.position = ifelse(y.position == min(y.position), y.position - y.min.down, y.position)) %>%
    mutate(xmin = xmin + xmin.max.shift, xmax = xmax + xmin.max.shift) %>% 
    return()
}

# define function to label with stars
func.sig2 <- function(plot){
  plot + 
    geom_text(data = s, # in global environment
              aes(x = (xmin + xmax)/2 , y = y.position, label = p.adj.signif), 
              inherit.aes = F, size = 8)
}

# plot RER
plt.RER <- func.plt.basicsPhysiology(whichY = "RER") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)),
                     breaks = seq(.7, 1, .05)) +
  coord_cartesian(ylim = c(0.7, NA)) 

plt.RER <- (plt.RER + 
              geom_segment(data = s <- d.infusion.stable.summary1 %>% ungroup() %>% 
                             pairwise_t_test(RER ~ phenotype, p.adjust.method = "bonferroni") %>% 
                             func.sig1(y.max.up = .02, y.min.down = .00),
                           aes(x = xmin, xend = xmax, y = y.position, yend = y.position), inherit.aes = F)) %>% 
  func.sig2()
plt.RER



# plot total CO2 output
plt.CO2 <- func.plt.basicsPhysiology(whichY = "CO2.umol.min", mytitle = "CO2 production") +
  scale_y_continuous(
    breaks = seq(0, 120, 25),limits = c(0, NA),
    expand = expansion(mult = c(0, .1)), name = "µmol / min" ,
    sec.axis = sec_axis(trans = ~. / 10^6 * 22.4 * 1000 * 60 ,
                        breaks = seq(0, 150, 30),
                        name = "mL / h")) 
plt.CO2 <- (plt.CO2 + 
              geom_segment(data = s <- d.infusion.stable.summary1 %>% ungroup() %>% 
                             pairwise_t_test(CO2.umol.min ~ phenotype, p.adjust.method = "bonferroni") %>% 
                             func.sig1(y.max.up = 5, y.min.down = 5),
                           aes(x = xmin, xend = xmax, y = y.position, yend = y.position), inherit.aes = F))  %>% 
  func.sig2()
plt.CO2


# plot total O2 consumption
plt.O2 <- func.plt.basicsPhysiology(whichY = "O2.umol.min", mytitle = "O2 consumption") +
  scale_y_continuous(
    breaks = seq(0, 150, 25), #  limits = c(0, 135), 
    expand = expansion(mult = c(0, .1)), name = "µmol / min",
    sec.axis = sec_axis(trans = ~. / 10^6 * 22.4 * 1000 * 60 ,
                        breaks = seq(0, 180, 30),
                        name = "mL / h"))
plt.O2 <- (plt.O2 +
             geom_segment(data = s <- d.infusion.stable.summary1 %>% ungroup() %>% 
                            pairwise_t_test(O2.umol.min ~ phenotype, p.adjust.method = "bonferroni") %>% 
                            func.sig1(y.max.up = 7, y.min.down = 7),
                          aes(x = xmin, xend = xmax, y = y.position, yend = y.position), inherit.aes = F))  %>% 
  func.sig2()
plt.O2


# func.plt.basicsPhysiology(whichY = "BW", mytitle = "BW")


# load basic animal information from the ordinary infusion experiments
p <- c("db/db", "WT", "HFD", "ob/ob")
names(p) <- c("d", "L", "H", "O")                     

mypath = "/Users/boyuan/Desktop/Harvard/Research/db db mice/Infusion data/infusion rounds data.xlsx"
d.animal  <-  read_excel(path = mypath, sheet = "infusion rounds", skip = 7) %>% 
  mutate(phenotype = str_extract(mouse_ID, "[a-zA-Z]+")) %>% 
  mutate(phenotype = p[phenotype] %>% factor(levels = p) ) 

# plot body weight using the serum-labeling experiment infusion 
d.BW2 <- d.animal %>% select(phenotype, mouse_ID, BW) %>% 
  group_by(mouse_ID, phenotype) %>% 
  summarise(BW = mean(BW)) %>% 
  filter(phenotype != "db/db") 

plt.BW <- func.plt.basicsPhysiology(
  mydata = d.BW2 %>% group_by(mouse_ID, phenotype) %>% summarise(BW = mean(BW)), 
  whichY = "BW", mytitle = "body weight") +
  scale_y_continuous(breaks = seq(0, 80, 15),limits = c(0, NA),
                     name = "gram",
                     expand = expansion(mult = c(0, .1)))
plt.BW <- (plt.BW +
             geom_segment(data = s <- d.BW2 %>% ungroup() %>% 
                            pairwise_t_test(BW ~ phenotype, p.adjust.method = "bonferroni") %>% 
                            add_xy_position(x = "phenotype") %>% 
                            func.sig1(y.max.up = 4, y.min.down = 4, xmin.max.shift = -1),
                          aes(x = xmin, xend = xmax, y = y.position, yend = y.position), inherit.aes = F))  %>% 
  func.sig2()
plt.BW



# Glycemia (from the 13CO2 measurement experiments)
d.glycemia <- d.mouseID %>% 
  filter(phenotype %in% c("WT", "ob/ob", "HFD")) %>% 
  mutate(phenotype = factor(phenotype, levels = c("WT", "HFD", "ob/ob", "db/db"))) %>% 
  mutate(glycemia = as.numeric(`glycemia.mg/dL-basal`)) 

x <- d.glycemia %>% filter(phenotype == "ob/ob") %>% group_by(round) %>% 
  summarise(glycemia = mean(glycemia, na.rm = T), phenotype = "ob/ob")

d.glycemia <- d.glycemia %>% filter(phenotype != "ob/ob") %>% bind_rows(x) %>% 
  mutate(phenotype = factor(phenotype, levels = ordered.phenotype))

# Glycemia [from the ordinary infusion (serum labeling)  experiments]
# d.glycemia <- d.animal %>% select(phenotype, mouse_ID, tail_blood_glucose_1) %>% 
#   group_by(mouse_ID, phenotype) %>% 
#   summarise(glycemia = mean(tail_blood_glucose_1, na.rm = T)) %>% 
#   filter(phenotype != "db/db") 
# d.glycemia <- d.glycemia[complete.cases(d.glycemia), ] %>% filter(between(glycemia, 100, 300))

plt.glycemia <- d.glycemia %>% 
  func.plt.basicsPhysiology(whichY = "glycemia") +
  scale_y_continuous(
    breaks = seq(0, 600, 100), limits = c(0, NA),
    name = "mg / dL", expand = expansion(mult = c(0, .1)),
    sec.axis = sec_axis(trans = ~ . / 180 * 10, name = "mM"))

plt.glycemia <- (plt.glycemia +
                   geom_segment(data = s <- d.glycemia %>% ungroup() %>%
                                  pairwise_t_test(glycemia ~ phenotype, p.adjust.method = "bonferroni") %>% 
                                  add_xy_position(x = "phenotype") %>% 
                                  func.sig1(),
                                aes(x = xmin, xend = xmax, y = y.position, yend = y.position), inherit.aes = F))  %>% 
  func.sig2()
plt.glycemia


# - total energy expenditure
d.energyExpenditure <- d.infusion.stable.summary1 %>% 
  mutate(O2.L.min = O2.umol.min / 10^6 * 22.4,
         CO2.L.min = CO2.umol.min / 10^6 * 22.4) %>% 
  # energy expenditure using Weir's formula
  mutate(kcal.h = (3.941 * O2.L.min + 1.106 * CO2.L.min) * 1440 / 24,
         cal.min = kcal.h / 60 * 1000 )

plt.calorie <- d.energyExpenditure %>% 
  func.plt.basicsPhysiology(whichY = "cal.min", mytitle = "energy expenditure") +
  scale_y_continuous(
    name = "cal / min",
    breaks = seq(0, 12, 3), 
    limits = c(0, NA), 
    expand = expansion(mult = c(0, .1)),
    sec.axis = sec_axis(trans = ~. * 60 / 1000,
                        name = "kcal / h"))
plt.calorie <- (plt.calorie +
                  geom_segment(data = s <- d.energyExpenditure %>% ungroup() %>% 
                                 pairwise_t_test(cal.min ~ phenotype, p.adjust.method = "bonferroni") %>% 
                                 func.sig1(y.max.up = 1, y.min.down = 1),
                               aes(x = xmin, xend = xmax, y = y.position, yend = y.position), inherit.aes = F))  %>% 
  func.sig2()
plt.calorie

# plot all together
b <- ggplot() + theme_void()
p1 <- plot_grid(plt.BW , b, plt.glycemia , b, plt.CO2, 
                nrow = 1, rel_widths = c(1, .1, 1, .1, 1))

p2 <- plot_grid(plt.O2 , b, plt.RER , b, plt.calorie, 
                nrow = 1, rel_widths = c(1, .1, 1, .1, 1))

plot_grid(p1, p2, nrow = 2, align = "v")

ggsave(filename = "basic physiology.pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 8, width = 14)



# ANCOVA analysis
d.energyExpenditure.selected <- d.energyExpenditure %>% ungroup() %>% 
  select(phenotype, mouseID, CO2.L.min, O2.L.min, cal.min, BW, RER) 


# plot
d.energyExpenditure.selected %>%
  ggplot(aes(BW, CO2.L.min, color = phenotype)) + 
  geom_point(size = 3, alpha = .5) +
  geom_smooth(method = "lm", se = F) +
  theme.myClassic +
  scale_color_manual(values = color.phenotype)

ggsave(filename = "BW_energy_expenditure_ANCOVA.pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 4, width = 5)

# ANCOVA test
model.glm <- glm(CO2.L.min ~ BW + phenotype, 
                 data = d.energyExpenditure.selected, 
                 family = gaussian(link = "identity"))
summary(model.glm)


# summarize infusion parameters - <>-- <>-- <>-- <>-- <>-- <>-- <>-- <>-- <>-- <>-
d.inj.info <-  d.all.raw2 %>% 
  filter(inj == "tail.vein_anesthesia" & treatment == "13C") %>% 
  group_by(tracer, treatment) %>% 
  summarise(inj.mM.max = max(`tracer mM`, na.rm = T),
            inj.mM.min = min(`tracer mM`, na.rm = T)) %>% 
  mutate(inj.vol.uL = 200) %>% 
  select(-treatment)

d.inf.info <-  d.mouseID %>% 
  filter(inj == "infusion" & treatment == "13C") %>% 
  group_by(tracer, treatment) %>% 
  summarise(inj.mM = mean(`tracer mM`, na.rm = T),
            inf = mean(inj.uL)) %>% ungroup() %>% 
  select(-treatment) %>% 
  rename(inf.mM = inj.mM) 

d.13C.info <- left_join(d.inj.info, d.inf.info) %>% as.data.frame()

# output useful information



# <>-# <>-# <>-# <>-# <>-# <>-# <>-# <>-# <>-# <>-# <>-# <>-# <>-# <>-# <>-

# BODY COMPOSITION

# dataset of body composition by MRI
d.bodyComposition <- read_excel(
  "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/body composition.xlsx") %>% 
  mutate(phenotype = factor(phenotype, levels = ordered.phenotype))


# mean and SD of mass of fat and lean tissues
d.bodyComposition.tidy <- d.bodyComposition %>%  
  group_by(phenotype) %>% 
  summarise(across(c(fat, lean),
                   .fn = list(mean = ~mean(.x), sd = ~sd(.x)))
  ) %>% 
  pivot_longer(-phenotype, names_sep = "_",
               names_to = c("part", ".value")) %>% 
  # error bar position
  group_by(phenotype) %>% 
  arrange(desc(part)) %>% 
  mutate(y.error = cumsum(mean))

# plot
plt.fat.lean <- d.bodyComposition.tidy %>% 
  ggplot(aes(x = phenotype, y = mean, fill = part)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = y.error - sd, ymax = y.error),
                width = .2) +
  theme.myClassic +
  scale_y_continuous(expand = expansion(mult = 0),
                     breaks = seq(0, 70, 10)) +
  scale_fill_manual(values = c("fat" = "snow2", "lean" = "snow4")) +
  labs(y = "mass (g)")

plt.fat.lean

# add labels
plt.fat.lean.labeled <- plt.fat.lean + geom_text(
  aes(label = round(mean, 1)), 
  position = position_stack(vjust = .5), color = "cyan") 
plt.fat.lean.labeled


# calculate fraction of fat and lean - ------------------------
d.bodyComposition.summary <- d.bodyComposition %>% 
  mutate(sum = fat + lean,
         fat.frac = fat/sum,
         lean.frac = lean / sum) %>% 
  group_by(phenotype) %>% 
  summarise(across(contains("frac"), 
                   .fn = list(
                     mean = ~ mean(.x, na.rm = T), # calcualte the mean
                     sd = ~ sd(.x, na.rm = T))))  # calculate standard deviation
# tidy up
d.bodyComposition.summary.tidy <- d.bodyComposition.summary %>% 
  pivot_longer(-phenotype, names_to = c("part", ".value"), names_sep = "_") %>% 
  group_by(phenotype) %>% 
  arrange(desc(part)) %>% 
  mutate(y.error = cumsum(mean)) %>% 
  mutate(part = str_remove(part, ".frac"))

# plot of distribution fraction
plt.fat.lean.frac <- d.bodyComposition.summary.tidy %>% 
  ggplot(aes(x = phenotype, y = mean, fill = part)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = y.error - sd, ymax = y.error),
                width = .2) +
  theme.myClassic +
  scale_y_continuous(expand = expansion(mult = 0),
                     breaks = seq(0, 1, .2)) +
  scale_fill_manual(values = c("fat" = "snow2", "lean" = "snow4")) +
  labs( y = "fraction")

plt.fat.lean.frac

# add label
plt.fat.lean.frac.labeled <- plt.fat.lean.frac + geom_text(
  aes(label = round(mean*100, 1)), 
  position = position_stack(vjust = .5), color = "cyan") 

plt.fat.lean.frac.labeled

# plot all
plot_grid(plt.fat.lean.labeled + theme(legend.position = "bottom"), 
          plt.fat.lean.frac.labeled + theme(legend.position = "bottom"))


# add labels



# export the data
save(d.energyExpenditure, func.plt.basicsPhysiology, d.13C.info, d.infusion.expCurve, 
     # body composition
     d.bodyComposition.summary.tidy, 
     d.bodyComposition.tidy,
     file = "3_core_13CO2_infusion_physiology_basics.RData")



