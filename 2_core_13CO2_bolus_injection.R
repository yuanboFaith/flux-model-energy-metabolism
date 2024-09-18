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

load("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/1_core_13CO2_import_data.RData")

theme.mybw <- theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = .5, face = "bold", size = 14),
        axis.title = element_text(size = 12, face = "bold"),
        panel.spacing = unit(10, "mm"))

# remove the clamp groups in rounds #99-103
d.all.raw3 <- d.all.raw2 %>% 
  filter(!(`glycemia.mg/dL-basal` %in% "clamp")) %>% 
  mutate(`glycemia.mg/dL-basal` = as.numeric(`glycemia.mg/dL-basal`))

#remove "old" glucose infusion data (3-4 infusion), and only keep 6h infusion
d.all.raw3 <- d.all.raw3 %>%
  filter(! (tracer == "Glucose" & inj == "infusion" & cage <= 972 ))


# define time = 0 upon tracer administration
d.all = d.all.raw3 %>% 
  mutate(inj.start = mdy_hms(inj.start)) %>% 
  # Now the time.h, starting from the point of tracer-injection
  group_by(round.cage) %>% 
  mutate(time.h = difftime(Time, inj.start, units = "h") %>% as.double()) %>% 
  filter(! (time.h >=  integrate.cutoff.time.h & inj == "tail.vein_anesthesia"))



d.all <- d.all %>% 
  # correct for water vapor pressure. This is critical for O2 and RER measurement
  # but not important for CO2
  mutate(O2_A = O2_A * 101.325 / (BP_A - WVP_A),
         CO2_A = CO2_A * 101.325 / (BP_A - WVP_A),
         SI_CO2_A = SI_CO2_A * 101.325 / (BP_A - WVP_A))

# # this section is optional, good for visualization -@#$%^&*()∑´®†¥¨∂ƒ©˙˜∆µ≈-@#$%^&*()∑´®†¥¨∂ƒ©˙˜∆µ≈
# 
# # for 13C-treated cage, BEFORE time = 0, change baseline delta 13C (empty cage measurement) to control level 
# # to help visualize 13C signal rise from control (delta 13C = -10) instead of ambient air (delta 13C = -30)
# d.control = d.all %>% filter(time.h >=0 & treatment == "control")
# ctrl.delta.13C.mean = d.control$SI_13C_A %>% mean()
# ctrl.delta.13C.sd = d.control$SI_13C_A %>% sd()
# 
# func.create.control.delta13C = function(){
#   rnorm(1, mean = ctrl.delta.13C.mean, ctrl.delta.13C.sd) 
# }
# 
# # For 13C-treated mice, before injection (time < 0 ), set the delta 13C to control level to help visualize
# # this does not affect integration, as integration has time >=0
# for(i in 1:nrow(d.all)){
#   if(d.all$treatment[i] == "13C" & d.all$time.h[i] <=0){
#     d.all$SI_13C_A[i] = rnorm(1, mean = ctrl.delta.13C.mean, ctrl.delta.13C.sd) 
#   }  
# }
# -@#$%^&*()∑´®†¥¨∂ƒ©˙˜∆µ≈-@#$%^&*()∑´®†¥¨∂ƒ©˙˜∆µ≈-@#$%^&*()∑´®†¥¨∂ƒ©˙˜∆µ≈-@#$%^&*()∑´®†¥¨∂ƒ©˙˜∆µ≈

d.all.raw2 %>% filter(treatment %in% c("baseline") )

(d.all %>% filter(treatment %in% c("baseline") & time.h > 0))$SI_13C_A %>% plot(main = "baseline delta 13C")
(d.all %>% filter(treatment %in% c("baseline") & time.h > 0))$SI_13C_A %>% hist(breaks = 200, main = "baseline delta 13C")
(d.all %>% filter(treatment %in% c("baseline") & time.h > 0))$SI_13C_A %>% mean()
(d.all %>% filter(treatment %in% c("baseline") & time.h > 0))$SI_13C_A %>% sd()

delta13C.control = (d.all %>% filter(treatment %in% c("control") & time.h > 0))$SI_13C_A %>% mean()



# remove outliers of tail vein injection
outlier.tailVein = c(282, # WT, glycerol, too low signal
                     247, # WT, glycerol, too high signal
                     238, 239,  # WT, alanine, too high signal,
                     184, 183, # WT, alanine, way low baseline to start with, with plateu signal start to drop after 2 h
                     286, # WT, glycerol, tanking at plateu due to signal lower than blank average
                     138, 143, 174, 176) # WT, C16:0, too high recovery
d.all = d.all %>% filter(! cage %in% outlier.tailVein) %>% 
  filter(! (phenotype == "HFD" & SI_13C_A < -40)) # remove a single downward spike in HFD glycerol i.v.injection at 3.5 h



# plot
plt.kinetics.delta.13C.tailVein.inj0 =  
  (iii <- d.all %>%
     # temporarily add baseline to phenotype, such that ambient air - phenotype both assigned with color
     mutate(phenotype  = ifelse(treatment == "baseline", "baseline", phenotype)) %>% 
     
     filter(inj == "tail.vein_anesthesia") %>% 
     filter(time.h >= -.5) ) %>% 
  filter(treatment %in% c("13C", "baseline")) %>% 
  filter(phenotype!= "baseline") %>% 
  ggplot( aes(x = time.h, y = SI_13C_A, color = phenotype)) + 
  geom_line(aes(group = round.cage), alpha = .6) +
  geom_hline(yintercept = delta13C.control, color = "firebrick") +
  
  # # control
  geom_line(data = iii %>% filter(treatment == "control"),
            aes(group = round.cage), color = "black", show.legend = F, alpha = 1)  +
  
  # baseline:  non-glucose, with deep color
  geom_line(data = iii %>% filter(phenotype == "baseline" & tracer != "Glucose"),
            aes(group = round.cage), show.legend = F, alpha = .7)  +
  # baseline: glucose, more transparent to show baseline centering at average level
  geom_line(data = iii %>% filter(phenotype == "baseline" & tracer == "Glucose"),
            aes(group = round.cage), show.legend = F, alpha = .3)  +
  
  facet_wrap(~tracer, nrow = 2, scales = "free") +
  theme.myClassic + 
  theme(legend.position = c(.92, .2), 
        legend.key.size = unit(.8, "cm"))  +
  guides(color = guide_legend(override.aes = list(linewidth = 3, alpha = 1))) +
  labs(y = "Delta 13C (per mill)", x = "Time (h)") +
  scale_color_manual(
    values = c("baseline" = "aquamarine4", "HFD" = "chocolate1", 
               "ob/ob" = "deepskyblue3", "WT" = "grey30"),
    limits = c("WT", "HFD", "ob/ob", "baseline"),
    labels = function(x){str_replace(x, "WT", "control")}) 
plt.kinetics.delta.13C.tailVein.inj0


# Add 13C atom injection dose annotation
d.dose.tailVeinInj =  iii %>% group_by(tracer) %>% 
  summarise(inj.umol.13C.atoms.min = min(inj.umol.13C.atoms, na.rm = T) %>% round(0),
            inj.umol.13C.atoms.max = max(inj.umol.13C.atoms, na.rm = T) %>% round(0),
            
            # inj dose x-y position, for facet-free scale
            maxY = max(SI_13C_A ), 
            maxX = integrate.cutoff.time.h %>% unique()) %>% 
  mutate(C13.atom.dose.umol = str_c(inj.umol.13C.atoms.min, " ~ ", inj.umol.13C.atoms.max ))

plt.kinetics.delta.13C.tailVein.inj <-  
  plt.kinetics.delta.13C.tailVein.inj0 + 
  #
  geom_text(data = d.dose.tailVeinInj,
            aes(x = maxX * .6, 
                # y = 200 * .8,  
                y = Inf,
                label = C13.atom.dose.umol), 
            inherit.aes = F , size = 5,
            vjust = 6) +
  #
  geom_text(data = d.dose.tailVeinInj,
            aes(x = maxX * .6, 
                #y = 200 * .65
                y = Inf), 
            inherit.aes = F,  
            label = "µmol 13C-atom \ni.v. injected", size = 4,
            color = "snow4",
            vjust = 4) +
  scale_x_continuous(breaks = seq(0, 6, 1)) +
  scale_y_continuous(n.breaks = 5) 

plt.kinetics.delta.13C.tailVein.inj

ggsave(filename = "delta 13C bolus.pdf", 
       plot = plt.kinetics.delta.13C.tailVein.inj, 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 7, width = 14)





# Check consistency of baseline and control delta 13C
my_viridis_palette <- viridis_pal(option = "D")(n = 6) # You can adjust 'n' to choose the number of colors you need.

plt.baseline.delta13C = d.all %>% filter(treatment %in% c("baseline")) %>% 
  filter(! is.na(tracer)) %>% 
  ggplot(aes(x = time.h, y = SI_13C_A, fill = as.numeric(date))) + 
  # geom_line(aes(group = round), alpha = .3, show.legend = F) +
  geom_point(size = 2, show.legend = T, alpha = .6, shape = 21) +
  scale_fill_viridis_c(labels = pretty(d.all$date), option = "B") +
  facet_wrap(~tracer,  scales = "free_x", nrow = 2) +
  labs(y = "delta 13C (per mill)",  x = "time (h)") +
  geom_hline(yintercept = (d.all %>% filter(treatment %in% c("baseline")))$SI_13C_A %>% median(), 
             size = .5) +
  coord_cartesian(xlim = c(0, 6), clip = "off") +  
  labs(title = "Baseline (ambient air), May 2022 ~ Dec 2023") +
  theme.myClassic +
  theme(panel.spacing = unit(10, units = "pt"),
        legend.position = "bottom",
        panel.border = element_rect(color = "black", fill = NA),
        legend.key.width = unit(100, units = "pt")) 
# plt.baseline.delta13C


plt.baseline.delta13C.densityCompiled = d.all %>%
  filter(treatment %in% c("baseline")) %>% 
  filter(! is.na(tracer)) %>% 
  ggplot(aes(x = time.h, y = SI_13C_A)) + 
  geom_bin2d(bins = 60) +
  scale_fill_continuous(type = "viridis") +
  labs(y = "delta 13C (per mill)",  x = "time (h)") +
  geom_hline(yintercept = (d.all %>% filter(treatment %in% c("baseline")))$SI_13C_A %>% median(), 
             size = 1, color = "black") +
  coord_cartesian(xlim = c(0, 6)) +
  labs(title = "All baseline compiled") +
  theme.myClassic +
  theme(legend.position = "bottom") +
  guides(fill = guide_colorbar(barwidth = unit(200, units = "pt")))
# plt.baseline.delta13C.densityCompiled

plt.baseline.delta13C.consistency2022 = plot_grid(
  plt.baseline.delta13C, 
  ggplot() + theme_void(), # white space
  plt.baseline.delta13C.densityCompiled, 
  rel_widths = c(6, .15,  3), nrow = 1, align = "h")

plt.baseline.delta13C.consistency2022

ggsave(filename = "13CO2 ambient air (baseline).pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       width = 14, height = 8)


# Calculate 13CO2 fraction based on delta 13C value
d.all = d.all %>% mutate(C13.C12.ratio = (SI_13C_A/1000 + 1 ) * VPDB, # 13C/12C ratio
                         C13.fraction = C13.C12.ratio / (1 + C13.C12.ratio)) %>% # 13CO2 fraction in total CO2
  select(-C13.C12.ratio) %>% select(-SI_13C_A) 




# Calculate 13C-natural abundance in influx air and diet & animal carboon pool, and CO2 ppm in influx air, averaged over entire experiment time
r.in = (d.all %>% filter(treatment == "baseline"))$C13.fraction %>% median() # 13CO2 fraction in total air CO2 
ppm.in = (d.all %>% filter(treatment == "baseline"))$SI_CO2_A %>% median() # total CO2 ppm in influx air
r.c.out = (d.all %>% filter(treatment == "control" ))$C13.fraction %>% median() #  13CO2 fraction in total CO2 of outflux air of control mice
ppm.c.out = (d.all %>% filter(treatment == "control" ))$SI_CO2_A %>% median() #  total CO2 ppm in outflux air of control mice
E2 = (ppm.c.out * r.c.out - ppm.in * r.in) / (ppm.c.out - ppm.in) # natural 13C abundance in diet and animal endogenous carboon pool
E1 = 1 # labeled glucose purity as 100%








# Calculate P1: the exhaled flux of carbons derived from the tracer
d.all.treated = d.all %>% filter(treatment != "baseline")
d.all.treated = d.all.treated %>%  # convert variable names to the notation in the Word File
  rename(ppm.t.out = SI_CO2_A,  # total CO2 ppm in the outflux air of 13C-treated mice
         r.t.out = C13.fraction) # 13CO2 fraction in the total CO2 in outflux air of 13C-treated mice)

d.all.treated = d.all.treated %>% 
  mutate(P1.mL.min = ( ppm.t.out * (r.t.out - E2) - ppm.in * (r.in - E2) ) * 2000 * 10^-6 / (E1 - E2)  ) %>%  # 13CO2 exhaled rate, mL/min
  mutate(P1.umol.min = P1.mL.min/1000 / 22.4 * 10^6) # 13CO2 exhaled rate, umol / min 




# To make trend smoother, add mouse endogenous 13CO2 signal at time.h = 0
# (this time point is usually missed as the cage air is not measured immedietly after mouse is put in cage)
d.time.h.zero = d.all.treated %>% filter(treatment == "13C" & inj == "tail.vein_anesthesia")

for (cage.i in d.time.h.zero$cage %>% unique()) {
  d.row.timeZero = (d.time.h.zero %>% filter(cage == cage.i & time.h < 10/60 & time.h > 0))[1, ]
  
  # Update time
  seconds.delta = (d.row.timeZero$time.h - 0) * 3600
  d.row.timeZero$Time =  d.row.timeZero$Time - seconds.delta
  
  # Update time.h
  d.row.timeZero$time.h = 0
  
  # Force P1 to zero
  d.row.timeZero$P1.mL.min = 0
  d.row.timeZero$P1.umol.min = 0
  
  # add this new row to the original dataset
  d.all.treated = rbind(d.all.treated, d.row.timeZero)
}






# Plot kinetic curves; Note the injection dose was different!
plt.kinetics.13C.flx.tailVein = d.all.treated %>% 
  filter(inj == "tail.vein_anesthesia") %>%
  filter(time.h >= -20/60) %>% 
  filter(treatment != "control") %>% 
  ggplot(aes(x = time.h, y = P1.umol.min, color = phenotype)) + 
  geom_hline(yintercept = 0) +
  # geom_line(data = d.all.treated %>% filter(treatment == "control"), color = "grey", aes(group = cage)) +
  geom_point(alpha = .7, size = 1.4) +
  geom_line(aes(group = cage), alpha = .7, size = .5) +
  # geom_smooth(data = d.all.treated %>% filter(time.h >= -.04 & treatment == "13C") ,
  #             aes(group = phenotype), se = F, span = .1, size = 1) +
  # geom_point(data = d.all.treated %>% filter(treatment == "control"), color = "grey", size = .1) +
  # facet_wrap(~cage, nrow = 2, scales = "free") +
  facet_wrap(~tracer, nrow = 2, scales = "free") +
  labs(y = "tracer-derived 13CO2 exhaled (umol/min/animal)") +
  scale_color_manual(values = c("HFD" = "chocolate1", "ob/ob" = "deepskyblue3", "WT" = "grey30")) 


plt.kinetics.13C.flx.tailVein # Note the injection dose was different, leading to different peak height!






# Normalize P1.umol.min based on a fixed 13C-atom injection number
normalization.inj.umol.13C.atoms = 10
d.all.treated = d.all.treated %>% 
  mutate(P1.umol.min.normalized = ifelse(P1.umol.min>=0, P1.umol.min * normalization.inj.umol.13C.atoms / inj.umol.13C.atoms,P1.umol.min))


# tail vein bolus injection, only treated mice data
d.all.treated.bolus.injection <- 
  d.all.treated %>% 
  filter(inj == "tail.vein_anesthesia") %>%
  filter(time.h >= -20/60) %>% 
  filter(treatment != "control")

func.plot.flux.bolus.injection <- function(mydata){
  mydata %>% 
    ggplot(aes(x = time.h, y = P1.umol.min.normalized, color = phenotype)) + 
    geom_point(alpha = .8, size = .7)  +
    geom_line(aes(group = cage), alpha = .5, linewidth = .2) +
    facet_wrap(~tracer, nrow = 2, scales = "free_x") +
    theme.myClassic +
    theme(legend.position = c(.9, .15),
          legend.key.size = unit(1, "cm")) +
    guides(col = guide_legend(override.aes = list(linewidth = 5, alpha = 1))) +
    labs(y = "13CO2 exhaled (nmol/min/animal), 
       10 µmol 13C-atom administered\n", x = "time (h)") +
    # geom_smooth(data = dd %>% filter(time.h >= 0), aes(group = phenotype),
    #             method = "loess", se = F, span = .1, size = 1) + # larger span number gives smoother curve fitting
    scale_color_manual(
      values = c( "chocolate1",  "deepskyblue3", "grey30"),
      labels = function(x){str_replace(x, "WT", "control")}) +
    scale_y_continuous(expand = expansion(mult = c(0, .1)),
                       # convert umol to nmol/min C
                       labels = function(x){x*1000})
}


# show bolus kinetics of all 3 phenotypes
plt.kinetics.13C.flx.tailVein.normalized.3Phenotype <- 
  func.plot.flux.bolus.injection(mydata = d.all.treated.bolus.injection)

plt.kinetics.13C.flx.tailVein.normalized.3Phenotype

ggsave(filename = "13CO2 output flux normalized_bolus_3phenotypes.pdf", 
       plot = plt.kinetics.13C.flx.tailVein.normalized.3Phenotype, 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 6, width = 12)

# show only WT 
plt.kinetics.13C.flx.tailVein.normalized.WT <- 
  func.plot.flux.bolus.injection(
    mydata = d.all.treated.bolus.injection %>% 
      filter(phenotype == "WT" # & tracer == "Glucose"
      )) +
  scale_color_manual(values = c("black")) +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "tomato")

plt.kinetics.13C.flx.tailVein.normalized.WT

ggsave(filename = "13CO2 output flux normalized_bolus_WT.pdf", 
       plot = plt.kinetics.13C.flx.tailVein.normalized.WT, 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 6, width = 12)


# show only WT and glucose 
plt.kinetics.13C.flx.tailVein.normalized.WT.glucose <- 
  func.plot.flux.bolus.injection(
    mydata = d.all.treated.bolus.injection %>% 
      filter(phenotype == "WT" & tracer == "Glucose"
      )) +
  scale_color_manual(values = c("black")) +
  theme(legend.position = "none", 
        axis.text = element_text(size = 17)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "tomato") +
  scale_y_continuous(breaks = seq(0, .12, .03),
                     expand = expansion(mult = c(0, .1)),
                     labels = function(x) {x * 1000})

plt.kinetics.13C.flx.tailVein.normalized.WT.glucose




# Calculate recovery
func.timesection = function(x){ c(diff(x), NA) }
d.all.treated = d.all.treated %>% 
  group_by(cage) %>% arrange(time.h) %>%  # this line is critical: check time section per cage, with time arranged in order
  mutate(timeSection.min = func.timesection(time.h) * 60) %>% 
  mutate(P1.umol.Cumulated.perTimeSection = P1.umol.min * timeSection.min)



CO2.13.bolus.recovered = d.all.treated %>% 
  filter(inj %in% c("tail.vein_anesthesia", "IP")) %>%
  mutate(cage = as.character(cage)) %>% 
  group_by(cage, round.cage) %>% 
  filter(time.h >= 0) %>% 
  summarise(P1.umols.integrated = sum(P1.umol.Cumulated.perTimeSection, na.rm = T)) %>%  # integrating all 13CO2 signal
  mutate(cage = as.double(cage)) %>% 
  left_join(d.mouseID, by = "cage") %>% 
  filter(treatment == "13C") %>% 
  mutate(recovery = P1.umols.integrated / inj.umol.13C.atoms)




# plot recovery
plt.recovery.bolus = CO2.13.bolus.recovered %>% 
  mutate(phenotype = factor(phenotype, levels = ordered.phenotype)) %>% 
  ggplot(aes(x = phenotype, y = recovery, color = cohort )) + 
  # geom_violin(aes(group = phenotype), color = "black") +
  stat_summary(geom = "crossbar", fun = "mean", aes(group = phenotype), color = "black") +
  geom_beeswarm(show.legend = T, cex = 4, size = 2) +
  # geom_text(aes(label = cage), show.legend = F) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, .1), 
                     expand = expansion(mult = c(0, 0.1)),
                     name = "Bolus recovery") + 
  scale_x_discrete(labels = function(x){str_replace(x, "WT", "control")}) +
  theme.myClassic +
  theme(axis.title.x = element_blank(), 
        panel.spacing = unit(10, units = "pt"),
        axis.text.x = element_text(angle = 50, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA)) +
  facet_wrap(~tracer, scales = "free_x", nrow = 2)

plt.recovery.bolus



# Summarise recovery
CO2.13.bolus.recovered.summary = CO2.13.bolus.recovered %>% 
  filter(inj == "tail.vein_anesthesia") %>%
  group_by(phenotype, tracer) %>% 
  summarise(recovery.mean = mean(recovery, na.rm = T),
            recovery.SD = sd(recovery, na.rm = T),
            n.rep = n_distinct(round.cage),
            recovery.SEM = recovery.SD / sqrt(n.rep) ) %>% 
  mutate(method = "bolus.inj")





# plot: recovery vs. glycemia
CO2.13.bolus.recovered %>%
  filter(inj == "tail.vein_anesthesia") %>%
  filter(phenotype %in% c("db/db", "ob/ob")) %>% 
  filter(tracer %in% c("Glucose", "Lactate")) %>%
  filter(! is.na(`glycemia.mg/dL-basal`)) %>% 
  ggplot(aes(x = `glycemia.mg/dL-basal`, y = recovery,  color = round.cage)) + 
  # ggplot(aes(x = recovery, y = `glycemia.mg/dL-post.lamp.heat`, color = round.cage)) + 
  geom_point(size = 4, cex = 2, show.legend = F) + 
  # geom_text_repel(aes(label = mouseID), show.legend = F) +
  geom_text_repel(aes(label = round.cage), show.legend = F) +
  # facet_wrap(~phenotype) +
  facet_grid( tracer ~ phenotype, scales = "free") +
  geom_smooth(aes(group = 1),
              method = "lm", se = F, show.legend = F, color = "black", size = .2) +
  theme_bw()



# plot: recovery vs. body weight
CO2.13.bolus.recovered %>% 
  filter(inj == "tail.vein_anesthesia") %>%
  filter(phenotype %in% c("db/db", "ob/ob")) %>% 
  filter(tracer %in% c("Glucose", "Lactate")) %>% 
  ggplot(aes(x = BW, y = recovery, color = round.cage)) + 
  geom_point(size=2, cex = 2, show.legend = F) + 
  # geom_text_repel(aes(label = round.cage), show.legend = F) +
  facet_wrap(~phenotype) +
  facet_grid( tracer ~ phenotype, scales = "free") +
  geom_smooth(method = "lm", se = F,  show.legend = F, color = "black", size = .2) +
  theme_bw()




#<>-@#$%^&*()™£¢∞§¶•#<>-@#$%^&*()™£¢∞§¶•#<>-@#$%^&*()™£¢∞§¶•#<>-@#$%^&*()™£¢∞§¶•

# calculate cumulated recovery over time

# Define function: sum all numbers (before and include) the current number along each element of a vector
func.cumulated.priorNumbers = function(x){
  xxx = vector()
  for (i in 1:length(x)){
    xxx = append(xxx, x[1:i]  %>% sum())
  }
  return(xxx)
}


d.all.treated.recovery.timeSections = d.all.treated %>% 
  mutate(cage = as.character(cage)) %>% 
  group_by(cage) %>% arrange(Time) %>% 
  
  # For the cumulated recovery over time, force few negative P1 values to zero to avoid confusing "decreasing" cumulated recovery over time
  # mutate(P1.umol.Cumulated.perTimeSection = ifelse(P1.umol.Cumulated.perTimeSection <0, 0, P1.umol.Cumulated.perTimeSection)) %>% 
  
  # sum the "P1.umol.Cumulated.perTimeSection" (each inividual time window) of time sections prior to the current section to calculate recovery curve over time
  mutate(P1.umol.Cumulated.until.current.TimeSection = func.cumulated.priorNumbers(P1.umol.Cumulated.perTimeSection)) %>% 
  # calculate recovery over time course
  mutate(recovery.timeSection = P1.umol.Cumulated.until.current.TimeSection / inj.umol.13C.atoms) %>% 
  
  # calculate recovery kinetics relatives to each mouse's own maximum (metabolic clearance rate)
  left_join(CO2.13.bolus.recovered %>% select(cage, P1.umols.integrated ) %>% mutate(cage = as.character(cage)),
            by = "cage") %>% 
  mutate(recovery.timeSection.normalized = P1.umol.Cumulated.until.current.TimeSection / P1.umols.integrated)


d.all.treated.recovery.timeSections.clean <- 
  d.all.treated.recovery.timeSections %>% 
  filter(inj == "tail.vein_anesthesia") %>%
  filter(phenotype %in% c("WT", "ob/ob", "HFD")) %>% 
  filter(treatment == "13C") %>% 
  filter(time.h >= -10/60) 

# show each individual mouse
plt.cumulatedRecovery.eachMouse =
  fun.plot.cumulatedRecovery.bolus <- function(
    mydata = d.all.treated.recovery.timeSections.clean){
    mydata %>% 
      ggplot(aes(x = time.h, y = recovery.timeSection, color = phenotype )) +
      geom_point(size = .1, alpha = .8) + 
      geom_line(aes(group = cage), size = .2, alpha = .8) +
      facet_wrap(~tracer, scales = "free_x", nrow = 2) +
      scale_color_manual(
        values = c("WT" = "grey30", "HFD" = "chocolate1",  "ob/ob" = "deepskyblue3"),
        limits = c("WT", "HFD", "ob/ob"),
        labels = function(x){str_replace(x, "WT", "control")}) +
      labs(y = "Accumulated recovery in exhaled CO2\n", x = "time (h)") +
      scale_y_continuous(breaks = seq(0, 1, .2), 
                         expand = expansion(mult = c(0, 0))) +
      scale_x_continuous(breaks = seq(0, 6, 1)) +
      coord_cartesian(ylim = c(0, 1)) +
      guides(color = guide_legend(override.aes = list(size = 5))) +
      theme.myClassic +
      theme(legend.position = c(.9, .2))
    # +
    #   geom_text(data = mydata %>% filter(time.h> 2.9), 
    #             aes(label = cage))
  }

# show all 3 phenotypes
plt.cumulatedRecovery.eachMouse.3pheno <- 
  fun.plot.cumulatedRecovery.bolus()
plt.cumulatedRecovery.eachMouse.3pheno

ggsave(filename = "13CO2 accumulated recovery_bolus.pdf", 
       plot = plt.cumulatedRecovery.eachMouse.3pheno, 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 6, width = 12)


# show only WT
plt.cumulatedRecovery.eachMouse.WT <- 
  fun.plot.cumulatedRecovery.bolus(
    mydata = d.all.treated.recovery.timeSections.clean %>% 
      filter(phenotype == "WT")) +
  scale_color_manual(values = "black") +
  theme(legend.position = "none") 
plt.cumulatedRecovery.eachMouse.WT

ggsave(filename = "13CO2 accumulated recovery_bolus_WT.pdf", 
       plot = plt.cumulatedRecovery.eachMouse.WT, 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 5, width = 10)


# show only WT with glucose injection
plt.cumulatedRecovery.eachMouse.WT.glucose <- 
  fun.plot.cumulatedRecovery.bolus(
    mydata = d.all.treated.recovery.timeSections.clean %>% 
      filter(phenotype == "WT" & tracer == "Glucose")) +
  scale_color_manual(values = "black") +
  theme(legend.position = "none",
        axis.text = element_text(size = 17)) +
  labs(y = "recovery\n") +
  scale_x_continuous(breaks = c(0, 2, 4, 6))
plt.cumulatedRecovery.eachMouse.WT.glucose



# Combine the two WT plots together
plot_grid(plt.kinetics.13C.flx.tailVein.normalized.WT.glucose ,
          ggplot() + theme_void(),
          plt.cumulatedRecovery.eachMouse.WT.glucose + labs(y = NULL),
          rel_widths = c(2, .1, 1.8), nrow = 1)

ggsave(filename = "13CO2 WT_glucose.pdf", 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 4, width = 8.3)



# relative to each mouse's own CO2 max-recovery
plt.cumulatedRecovery.eachMouse.normalizedToMax = 
  d.all.treated.recovery.timeSections %>% 
  filter(inj == "tail.vein_anesthesia") %>%
  filter(treatment == "13C") %>% 
  filter(time.h >= -30/60) %>% 
  ggplot(aes(x = time.h, y = recovery.timeSection.normalized, 
             color = phenotype
  )) +
  geom_point(size = .5, alpha = .5) + 
  geom_hline(yintercept = 1) +
  geom_line(aes(group = cage), size = .5, alpha = .5) + 
  # scale_color_viridis(option = "C") +
  labs(y = "recovery of oxidizable carbons\n") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2), expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(breaks = seq(0, 6, 1), expand = expansion(mult = c(0, 0))) +
  facet_wrap(~tracer, scales = "free_x", nrow = 2)   +
  scale_color_manual(values = color.phenotype) +
  theme_bw()

plt.cumulatedRecovery.eachMouse.normalizedToMax 



# save project
save.image(file = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/2_core_13CO2_bolus_injection.RData")


