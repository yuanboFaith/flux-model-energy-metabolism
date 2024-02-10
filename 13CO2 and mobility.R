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



# Global theme setup ========
theme_set(theme_bw() +
            theme(title = element_text(face = "bold", colour = "black", size = 14),
                  plot.title = element_text(size = 14, hjust = .5),
                  plot.subtitle = element_text(size = 10, hjust = .5), 
                  axis.text = element_text(colour = "black", size = 14),
                  strip.background = element_blank(),
                  panel.grid = element_blank(),
                  legend.text = element_text(size = 16),
                  strip.text = element_text(face = "bold", size = 16),
                  panel.spacing = unit(7, "mm"),
                  axis.title.y = element_text(margin = margin(r = 14, "pt")),
                  axis.title.x = element_text(margin = margin(t = 10, "pt"))))

theme.myClassic <-  theme_classic() + 
  theme(strip.background = element_blank(),
        axis.text = element_text(colour = "black", size=13),
        axis.title = element_text(colour = "black", size = 13, face = "bold"),
        strip.text = element_text(colour = "black", size = 14, face = "bold"),
        legend.text = element_text(colour = "black", size = 12),
        plot.title = element_text(size = 14, hjust = .5, face = "bold"),
        legend.title = element_blank(),
        panel.spacing = unit(10, "mm"),
        axis.title.y = element_text(margin = margin(r = 14, "pt")),
        axis.title.x = element_text(margin = margin(t = 10, "pt")))


ordered.phenotype = factor(c("WT",  "HFD", "ob/ob", "db/db", "WT-KD"), ordered = T)
color.phenotype = c("grey30", "chocolate1", "deepskyblue3", "firebrick3",  "darkgreen")
color.phenotype %>% show_col()
names(color.phenotype) = ordered.phenotype

ordered.Compound = c("Bicarbonate", "Glucose", "Lactate",  "Glycerol","Alanine", "Glutamine", "3-HB", "Acetate", "C16:0", "C18:1", "C18:2", "Valine")

# Import data and experimental setup  (NEED MANUAL UPDATE)
pathID = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/13CO2-MouseID.xlsx"
d.mouseID = read_excel(pathID, sheet = "mouseID")
d.integrate.cutOff = read_excel(pathID, sheet = "integrate.cutOff")

# 1-6: inject by catheter
# path1 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2022-May-13C glucose bolus r1.xlsx" # 2022-05-16 experiment
# path2 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2022-May-13C glucose bolus r2.xlsx" # 2022-05-17 experiment
# path3 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2022-May-13C glucose bolus r3.xlsx" # 2022-05-18 experiment
path4 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2022-May-13C glucose bolus r4_control.xlsx" # 2022-05-19 experiment
# path5 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2022-May-13C glucose bolus r5.xlsx" # 2022-05-31 experiment
# path6 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2022-May-13C glucose bolus r6.xlsx" # 2022-06-01 experiment

# starting from round 7, tail vein injection
path7 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2022-May-13C glucose bolus r7.xlsx" # 2022-07-12 experiment
path8 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2022-May-13C glucose bolus r8.xlsx" # 2022-07-13 experiment
path9 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2022-May-13C glucose bolus r9.xlsx" # 2022-07-14 experiment
path10 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2022-May-13C glucose bolus r10.xlsx" # 2022-07-15 experiment
path11 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2022-May-13C glucose bolus r11.xlsx" # 2022-07-21 experiment
path12 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2022-May-13C glucose bolus r12.xlsx" # 2022-07-22 experiment
path13 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2022-May-13C glucose bolus r13.xlsx" # 2022-07-25 experiment
path14 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2022-May-13C glucose bolus r14.xlsx" # 2022-07-26 experiment
path15 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2022-May-13C glucose bolus r15.xlsx" # 2022-07-27 experiment
path89 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2023-May_glucose_13C_r89.xlsx"

# round 90 data was not recorded correctly
path91 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2023-May_glucose_13C_r91.xlsx"
path91.2 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2023-May_glucose_13C_r91-2.xlsx"

path92 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2023-May_glucose_13C_r92.xlsx"
path93 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2023-May_glucose_13C_r93.xlsx"
path94 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2023-May_glucose_13C_r94.xlsx"
path95 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2023-May_glucose_13C_r95.xlsx"
path96 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2023-May_glucose_13C_r96.xlsx"
path97 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2023-May_glucose_13C_r97.xlsx"
path98 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2023-May_glucose_13C_r98.xlsx"


path99 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2023-May_glucose_13C_r99.xlsx"
path100 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2023-May_glucose_13C_r100.xlsx"
path101 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2023-May_glucose_13C_r101.xlsx"
path102 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2023-May_glucose_13C_r102.xlsx"
path103 = "/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/2023-May_glucose_13C_r103.xlsx"

# d.expeData.1 = read_excel(path1)
# d.expeData.2 = read_excel(path2)
# d.expeData.3 = read_excel(path3)
d.expeData.4.control = read_excel(path4)
# d.expeData.5 = read_excel(path5)
# d.expeData.6 = read_excel(path6)
d.expeData.7 = read_excel(path7)
d.expeData.8 = read_excel(path8)
d.expeData.9 = read_excel(path9)
d.expeData.10 = read_excel(path10)
d.expeData.11 = read_excel(path11)
d.expeData.12 = read_excel(path12)
d.expeData.13 = read_excel(path13)
d.expeData.14 = read_excel(path14)
d.expeData.15 = read_excel(path15)

d.expeData.89 = read_excel(path89)

# round #90 was not recorded correctly
d.expeData.91 = read_excel(path91)
d.expeData.91.2 = read_excel(path91.2)
d.expeData.91 <- rbind(d.expeData.91, d.expeData.91.2)

d.expeData.92 = read_excel(path92)
d.expeData.93 = read_excel(path93)
d.expeData.94 = read_excel(path94)
d.expeData.95 = read_excel(path95)
d.expeData.96 = read_excel(path96)
d.expeData.97 = read_excel(path97)
d.expeData.98 = read_excel(path98)

# new clamp - the definite one
d.expeData.99 = read_excel(path99)
d.expeData.100 = read_excel(path100)
d.expeData.101 = read_excel(path101)
d.expeData.102 = read_excel(path102)
d.expeData.103 = read_excel(path103)

# Other parameters, cited globally in project
VPDB = 0.01123720



# Define function: select only equilibrated points, and remove spike point (outliers)
func.equilibrate.outLierRemoved = function(myData){
  #myData = d.expeData.8
  d.expeData = myData %>% 
    rename(Time = DateTime, cage = RespCage_A) %>% 
    mutate(cage = as.double(cage),
           time.h = as.numeric(Time - min(Time)) / 3600  ) # time since the start of data recording
  
  
  # Collect only fully equilibrated data points 
  ## Clean data; this clean-up section can be safely re-run for dataset update after removal # trial and error
  index.rowsCollected = c() # Collect number of rows
  SectionSartMarker = 0
  cages = d.expeData$cage
  for(i in 1: (nrow(d.expeData) - 1 ) ){
    if (cages[i] != cages[i+1]){
      
      #index.rowsCollected = append(index.rowsCollected, c((i-1):(i)) )
      index.rowsCollected = append(index.rowsCollected, c((i)) )  # MANUALLY adjust # of points to include at then end of each dwell time
      # Collect the last # of points during each dwell time. e.g., i-0 : i to collect the last data point at the end of each dwell time
      # It is realized that occasionally a dwell time period does not contain expected number of data points due to instrument error. 
      # Collecting the last XX points for each defined dwell time is robust to this error. The last XX points are visually examined to reach stable state in the above timeline offset procedure
      
      SectionSartMarker = i+1
    }
  }
  d.CO2.equilibrated = d.expeData[index.rowsCollected, ] %>% as_tibble()
  
  
  # 3)
  # Remove outlier spike points 
  # Set up smoothing spline for reach cage 13C readout
  mdl.spline.13C.time = d.CO2.equilibrated %>% nest(-cage) %>%
    # make model and predict
    mutate(mdl = map(data, ~smooth.spline(x = .$time.h, y = .$SI_13C_A)),
           fitted = map2(.x = mdl, .y = data,  ~ predict(.x, x = .y$time.h)$y )) # predict output: $x, time line; $y, 13C fitted value
  
  # Unnest the data
  d.mdl.spline.13C.time.tidy = tibble()
  
  for (i in 1:nrow(mdl.spline.13C.time)){ # loop through each cage and unnest the data
    d.i = mdl.spline.13C.time[i, ]
    output.i = tibble(cage = d.i$cage, # cage number
                      time.h = d.i$data [[1]] [["time.h"]],
                      SI_13C_A = (d.i$data) [[1]] [["SI_13C_A"]], # actual 13C value
                      SI_13C_A.fitted = (d.i$fitted) [[1]]
    )
    d.mdl.spline.13C.time.tidy = rbind(
      d.mdl.spline.13C.time.tidy, output.i) # combine unnested data from each cage
  }
  
  # calculate residual
  d.mdl.spline.13C.time.tidy = d.mdl.spline.13C.time.tidy %>% 
    mutate(residual = SI_13C_A.fitted - SI_13C_A)
  
  # Mark outliers with 3-std away from the overall trend
  d.RMSE = d.mdl.spline.13C.time.tidy %>%
    group_by(cage) %>% 
    summarise(n.rows = length(cage),
              # RMSE = ( sum(residual^2) / n.rows ) %>% sqrt()  ) %>%
              RMSE = ( sum( abs(residual) ) / n.rows )) %>%     
    # use mean of absolute residual, not root mean squared error, to reduce sensitivity to outlier
    # this leads to smaller magnitude of estimated residual, allowing for more delicate control of outlier removal 
    # mark outliers that are 6 RMSE away from the fitted value
    mutate(error.tolerance.cutoff = 6 * RMSE)
  
  d.mdl.spline.13C.time.tidy = d.mdl.spline.13C.time.tidy %>% 
    left_join(d.RMSE, by = "cage") %>%
    mutate(outlier = abs(residual) > error.tolerance.cutoff  )
  
  
  d.CO2.equilibrated.outlierRemoved = d.CO2.equilibrated %>% 
    left_join(d.mdl.spline.13C.time.tidy %>% select(cage, time.h, SI_13C_A, outlier),
              by = c("cage", "time.h", "SI_13C_A")) %>% 
    filter(outlier == F) %>% 
    select(-outlier) %>% 
    arrange(time.h) 
  
  # mark cage #0 as baseline
  # mutate(is.baseline = ifelse(cage == 0, T, F))
  
  # Check
  # d.CO2.equilibrated.outlierRemoved %>% 
  #   ggplot(aes(x = time.h, y = SI_13C_A, color = as.character(cage))) + geom_line() 
  
  return(d.CO2.equilibrated) # output equilibrated, no spike-removal function executed
  # return(d.CO2.equilibrated.outlierRemoved) # output equilibrated and spike-removed data
}

# cage number accumulated to be unique and in numerical order; needs to be consistent with Excel mouse ID cage number!!
d4 = func.equilibrate.outLierRemoved(myData = d.expeData.4.control) %>% mutate(cage = (cage + 3 * 9),  round = 4)
d8 = func.equilibrate.outLierRemoved(myData = d.expeData.8) %>% mutate(cage = (cage + 7 * 9),  round = 8)
d9 = func.equilibrate.outLierRemoved(myData = d.expeData.9) %>% mutate(cage = (cage + 8 * 9),  round = 9)
d10 = func.equilibrate.outLierRemoved(myData = d.expeData.10) %>% mutate(cage = (cage + 9 * 9),  round = 10)
d11 = func.equilibrate.outLierRemoved(myData = d.expeData.11) %>% mutate(cage = (cage + 10 * 9),  round = 11)
d12 = func.equilibrate.outLierRemoved(myData = d.expeData.12) %>% mutate(cage = (cage + 11 * 9),  round = 12)
d13 = func.equilibrate.outLierRemoved(myData = d.expeData.13) %>% mutate(cage = (cage + 12 * 9),  round = 13)
d14 = func.equilibrate.outLierRemoved(myData = d.expeData.14) %>% mutate(cage = (cage + 13 * 9),  round = 14)
d15 = func.equilibrate.outLierRemoved(myData = d.expeData.15) %>% mutate(cage = (cage + 14 * 9),  round = 15)
d89 = func.equilibrate.outLierRemoved(myData = d.expeData.89) %>% mutate(cage = (cage + 88 * 9),  round = 89)

# round 90 data was not recorded correctly
d91 = func.equilibrate.outLierRemoved(myData = d.expeData.91) %>% mutate(cage = (cage + 90 * 9),  round = 91)
d92 = func.equilibrate.outLierRemoved(myData = d.expeData.92) %>% mutate(cage = (cage + 91 * 9),  round = 92)
d93 = func.equilibrate.outLierRemoved(myData = d.expeData.93) %>% mutate(cage = (cage + 92 * 9),  round = 93)
d94 = func.equilibrate.outLierRemoved(myData = d.expeData.94) %>% mutate(cage = (cage + 93 * 9),  round = 94)
d95 = func.equilibrate.outLierRemoved(myData = d.expeData.95) %>% mutate(cage = (cage + 94 * 9),  round = 95)
d96 = func.equilibrate.outLierRemoved(myData = d.expeData.96) %>% mutate(cage = (cage + 95 * 9),  round = 96)
d97 = func.equilibrate.outLierRemoved(myData = d.expeData.97) %>% mutate(cage = (cage + 96 * 9),  round = 97)
d98 = func.equilibrate.outLierRemoved(myData = d.expeData.98) %>% mutate(cage = (cage + 97 * 9),  round = 98)

# new clamp - the definite one
d99 = func.equilibrate.outLierRemoved(myData = d.expeData.99) %>% mutate(cage = (cage + 98 * 9),  round = 99)
d100 = func.equilibrate.outLierRemoved(myData = d.expeData.100) %>% mutate(cage = (cage + 99 * 9),  round = 100)
d101 = func.equilibrate.outLierRemoved(myData = d.expeData.101) %>% mutate(cage = (cage + 100 * 9),  round = 101)
d102 = func.equilibrate.outLierRemoved(myData = d.expeData.102) %>% mutate(cage = (cage + 101 * 9),  round = 102)
d103 = func.equilibrate.outLierRemoved(myData = d.expeData.103) %>% mutate(cage = (cage + 102 * 9),  round = 103)

# Merge all files
ds = list(d4, d8, d9, d10, 
          d11, d12, d13, d14, d15, 
          # insulin vs. glucose fox
          d89, d91, d92, d93, d94, d95, d96, d97, d98,  # round 90 data was not recorded correctly
          
          # new glucose clamp - the definite one
          d99, d100, d101, d102, d103
) 


d.all = tibble()
for(i in 1: length(ds)){ d.all = rbind(d.all, ds[[i]] ) }




# Combine with mouse ID dataset
d.all = d.all %>% 
  left_join(d.mouseID, by = c("cage", "round")) %>%
  mutate(round.cage = str_c(round, "_", cage)) %>% 
  # combine with integration cutoff time
  left_join(d.integrate.cutOff, by = "tracer") %>% 
  mutate(tracer = factor(tracer, levels = ordered.Compound, ordered = T))


# define time = 0 upon tracer administration
d.all = d.all %>% mutate(inj.start = mdy_hms(inj.start)) %>% 
  # Now the time.h, starting from the point of tracer-injection / infusion 
  group_by(round.cage) %>% 
  mutate(time.h = difftime(Time, inj.start, units = "h") %>% as.double()) %>% 
  filter(! (time.h >=  integrate.cutoff.time.h & inj == "tail.vein_anesthesia"))


# correct for water vapor pressure
d.all <- d.all %>% mutate(O2_A = O2_A * 101.325 / (BP_A - WVP_A))

delta13C.control = (d.all %>% filter(treatment %in% c("control") & time.h > 0))$SI_13C_A %>% mean()


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
  mutate(P1.umol.min = P1.mL.min/1000 / 22.4 * 10^6) %>%  # 13CO2 exhaled rate, umol / min 
  
  # total pure CO2 output rate from animal
  mutate(P_total.mL.min = (ppm.t.out - ppm.in) * 2000 / 10^6,
         P_total.umol.min = P_total.mL.min/1000 / 22.4 * 10^6)




# calculate RER

d.O2.baseline <- d.all %>% 
  # filter(round %in% c(99:103) & time.h > 2) %>% 
  select(round, O2_A, BP_A, WVP_A,  mouseID, time.h, treatment, treatment2) %>% 
  filter(treatment == "baseline") 

d.O2.baseline.mean <- d.O2.baseline %>% 
  group_by(round) %>% 
  summarise(O2.baseline.mean = mean(O2_A))

d.O2.baseline %>% ggplot(aes(time.h, O2_A)) + geom_point()

d.O2.mice <- d.all %>% 
  filter(round %in% c(99:103)) %>% 
  filter(treatment != "baseline")



# O2 consumption, subtracting from air separately for each round
d.O2.consume <- d.O2.mice %>% 
  select(round, O2_A, mouseID, time.h, treatment) %>% 
  left_join(d.O2.baseline.mean) %>% 
  mutate(O2.consume.mL.min = (O2.baseline.mean - O2_A)/100 * 2000,
         O2.consume.umol.min = O2.consume.mL.min/1000/22.4 * 10^6 ) %>% 
  select(-c(O2.baseline.mean, O2_A, contains("mL")))

d.clamp <- d.all.treated %>% 
  filter(round %in% c(99:103)) %>% 
  select(time.h, mouseID, P_total.umol.min, P1.umol.min, inj.umol.13C.atoms, treatment2) %>% 
  
  # join O2 consumption
  left_join(d.O2.consume) %>% 
  
  # infusion failure
  filter(round.cage != "99_887") %>% 
  
  mutate(RER = P_total.umol.min / O2.consume.umol.min,
         f_ox = P1.umol.min / inj.umol.13C.atoms,
         
         # glucose carbons contribution to total CO2 production
         # .74 is an empirical correction factor, so that the estiamted Fcirc match glucose infusion rate 
         f_glu.CO2 = (10 * RER - 7) / 3 / RER   * .74, 
         f_FA.CO2 = 1 - f_glu.CO2) %>% 
  
  # Fcirc estimate
  mutate(Fcirc.glu.umol.min = P_total.umol.min * f_glu.CO2 / f_ox,
         Fcirc.FA.umol.min = P_total.umol.min * f_FA.CO2 / .4) %>%  # assume 40% fox for fatty acids
  
  # total carbon oxidation flux
  mutate(F.glu.CO2.umol.min = Fcirc.glu.umol.min * f_ox,
         F.FA.CO2.umol.min = Fcirc.FA.umol.min * .4)

d.clamp


d.clamp.post2h <- d.clamp %>% filter(time.h > 2.5)


colnames(d.clamp.post2h)
d.clamp.post2h.averaged <- d.clamp.post2h %>% 
  group_by(round.cage, treatment2, ) %>% 
  summarise(across(c("RER", "f_ox", "P_total.umol.min", 
                     "f_glu.CO2", "f_FA.CO2", "F.FA.CO2.umol.min", "F.glu.CO2.umol.min",
                     "O2.consume.umol.min",
                     "Fcirc.glu.umol.min", "Fcirc.FA.umol.min"),
                   mean))


pp3 <- d.clamp.post2h.averaged %>% 
  ggplot(aes(treatment2, y = P_total.umol.min, fill = treatment2)) +
  stat_summary(geom = "bar") +
  stat_summary(geom = "text", aes(label = after_stat(y) %>% round(1)), size = 8) +
  geom_beeswarm(cex = 8) +
  scale_x_discrete(limits = c("control", "clamp")) +
  coord_cartesian(ylim = c(0, 69)) +
  theme(legend.position = "none")

pp3




# based on clamp, total glucose oxidation flux = infusion flux * fox
# which is 250 nmol/min/g * 30g / 1000 * 6 carbons * .82 fox = 37
# 37/50 = .74 correction factor

# combine with end point glycemia
d.glycemia.endpoint <- d.all %>% filter(round %in% 99:103) %>% 
  select(Glycemia_end, round.cage) %>% distinct() 

d.clamp.post2h <- d.clamp.post2h %>% left_join(d.glycemia.endpoint)



# add mobility data
d.speed.99 <- read.csv("/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/pedSpeed/pedSpeed_r99.csv")
d.speed.100 <- read.csv("/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/pedSpeed/pedSpeed_r100.csv")
d.speed.101 <- read.csv("/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/pedSpeed/pedSpeed_r101.csv")
d.speed.102 <- read.csv("/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/pedSpeed/pedSpeed_r102.csv")
d.speed.103 <- read.csv("/Users/boyuan/Desktop/Harvard/Research/db db mice/13CO2 analyzer/pedSpeed/pedSpeed_r103.csv")

func.tidy.pedSpeed <- function(myData, whichRound){
  # myData <- d.speed.103; whichRound = 103
  myData %>% as_tibble() %>% 
    mutate(round = whichRound) %>% 
    left_join(d.all %>% ungroup() %>%  select(round, inj.start, round) %>% distinct()) %>% 
    mutate(DateTime = as.character(DateTime),
           inj.start = as.character(inj.start)) %>% 
    mutate(time.h = difftime(DateTime, inj.start, units = "hours") %>% as.double() ) %>% 
    select(-c(DateTime, inj.start, ElSeconds)) %>% 
    pivot_longer(-c(round, time.h), names_to = "metCage", values_to = "pedSpeed") %>% 
    mutate(metCage = str_remove(metCage, "PedSpeed_"))
}

d.speed.tidy.99 <- func.tidy.pedSpeed(myData = d.speed.99, whichRound = 99)
d.speed.tidy.100 <- func.tidy.pedSpeed(myData = d.speed.100, whichRound = 100) 
d.speed.tidy.101 <- func.tidy.pedSpeed(myData = d.speed.101, whichRound = 101) 
d.speed.tidy.102 <- func.tidy.pedSpeed(myData = d.speed.102, whichRound = 102) 
d.speed.tidy.103 <- func.tidy.pedSpeed(myData = d.speed.103, whichRound = 103) 

d.speed <- rbind(d.speed.tidy.99, d.speed.tidy.100) %>% rbind(d.speed.tidy.101) %>% 
  rbind(d.speed.tidy.102) %>% rbind(d.speed.tidy.103) %>% 
  mutate(round.cage = str_c(round, "_", ((round-1) * 9 + as.double(metCage)) )) %>% 
  select(-c(round, metCage)) 





x1 <- d.speed$time.h %>% unique() %>% round(1)
x2 <- d.clamp$time.h %>% unique() %>% round(1)
unique(x1) %in% unique(x2)


d.speed.timeInterval <- d.speed %>% 
  mutate(time.h.interval = cut_width(time.h, width = .1, boundary = 0)) %>% select(-time.h)

d.clamp.speed <-  d.clamp %>% # filter(time.h > 2) %>% 
  filter(round %in% c(99:103)) %>% 
  mutate(time.h.interval = cut_width(time.h, width = .1, boundary = 0)) %>% 
  
  # join with speed data
  left_join(d.speed.timeInterval,
            by = c("time.h.interval", "round.cage"))


p.speed.13CO2 <- d.clamp.speed %>% 
  filter(treatment2 == "control") %>% 
  ggplot(aes(time.h, P1.umol.min / max(P1.umol.min, na.rm = T))) + 
  annotate(geom = "rect",
           xmin = 2.5, xmax = Inf, 
           ymin = -Inf, ymax = Inf, 
           alpha = .1, fill = "green") +
  geom_line(linewidth = .7) +
  geom_line(aes(y = pedSpeed / max(pedSpeed, na.rm = T)), 
            linewidth = .7, color = "red") +
  facet_wrap(~round.cage) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "green") +
  labs(y = "Relative intensity", x = "time (h)") +
  theme(strip.text = element_blank(),
        axis.title.y = element_text(margin = margin(r = 14, "pt")),
        axis.title.x = element_text(margin = margin(t = 10, "pt"))) +
  coord_cartesian(expand = 0)

p.speed.13CO2

ggsave(filename = "13CO2 output vs. activity.pdf", 
       plot = p.speed.13CO2, 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 6, width = 10)


d.clamp.speed2 <- d.clamp.speed %>% 
  # mutate(P1.umol.min.norm = P1.umol.min / max(P1.umol.min, na.rm = T),
  #        pedSpeed.norm = pedSpeed / max(pedSpeed, na.rm = T)) %>% 
  filter(treatment2 == "control") %>% 
  filter(time.h >= 2.5) %>% 
  filter(! round.cage %in% c("99_883", "99_884", "99_885"))

p.P1.activity <- d.clamp.speed2 %>% 
  ggplot(aes(x = P1.umol.min,
             y = pedSpeed)) + 
  geom_point(alpha = .4, size = 4) +
  labs(x = "13CO2 exhaled (µmol/min)",
       y = "Pedestrian speed")

p.P1.activity

library("quantreg")
fm <- rq(data = d.clamp.speed2,
         pedSpeed ~ P1.umol.min)

fm$coefficients

fm %>% summary()

p.P1.activity <- p.P1.activity +
  geom_abline(slope = fm$coefficients[2],
              intercept = fm$coefficients[1]) +
  theme.myClassic
p.P1.activity

ggsave(filename = "13CO2 output vs. activity scatterplot.pdf", 
       plot = p.P1.activity, 
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       device = "pdf", height = 6, width = 4)


d.clamp.speed %>% 
  filter(treatment2 == "control") %>% 
  ggplot(aes(time.h, O2.consume.umol.min / max(O2.consume.umol.min))) + 
  geom_line() +
  geom_line(aes(y = pedSpeed / max(pedSpeed, na.rm = T)), color = "red") +
  facet_wrap(~round.cage) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "green") +
  labs(y = "O2 consumption rate (umol /min)")




