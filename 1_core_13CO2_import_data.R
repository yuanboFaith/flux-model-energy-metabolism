rm(list = ls())

library(plyr)
library(xlsx)
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
theme_set(
  theme_bw() +
    theme(title = element_text(face = "bold", colour = "black"),
          plot.title = element_text(hjust = .5),
          plot.subtitle = element_text(hjust = .5), 
          axis.text = element_text(colour = "black"),
          strip.background = element_blank(),
          panel.grid = element_blank(),
          legend.text = element_text(),
          strip.text = element_text(face = "bold"),
          panel.spacing = unit(7, "mm")))

theme.myClassic <-  
  theme_classic(base_size = 16) +  
  theme(strip.background = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.title = element_text(colour = "black",face = "bold"),
        strip.text = element_text(colour = "black", face = "bold"),
        legend.text = element_text(colour = "black"),
        plot.title = element_text(hjust = .5, face = "bold"),
        legend.title = element_blank(),
        panel.spacing = unit(10, "mm"))


ordered.phenotype = factor(c("WT",  "HFD", "ob/ob"), ordered = T)
color.phenotype = c("grey30", "chocolate1", "deepskyblue3")
color.phenotype %>% show_col()
names(color.phenotype) = ordered.phenotype

ordered.Compound = c(
  "Bicarbonate", "Glucose", "Lactate",  "Glycerol","Alanine", "Glutamine", 
  "3-HB", "Acetate", "C16:0", "C18:1", "C18:2", "Valine", "Methionine", "Tryptophan")

# Import data and experimental setup  (NEED MANUAL UPDATE)
pathID = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2_mouseID_paper.xlsx"
d.mouseID = read_excel(pathID, sheet = "mouseID")
d.integrate.cutOff = read_excel(pathID, sheet = "integrate.cutOff")

# 1-6: inject by catheter
path4 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r4_control.xlsx" # 2022-05-19 experiment
path9 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r9.xlsx" # 2022-07-14 experiment
path10 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r10.xlsx" # 2022-07-15 experiment
path11 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r11.xlsx" # 2022-07-21 experiment
path12 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r12.xlsx" # 2022-07-22 experiment
path13 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r13.xlsx" # 2022-07-25 experiment
path14 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r14.xlsx" # 2022-07-26 experiment
path15 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r15.xlsx" # 2022-07-27 experiment
path16 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r16.xlsx" # 2022-07-28 experiment
path17 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r17.xlsx" # 2022-07-29 experiment
path18 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r18.xlsx" # 2022-07-30 experiment
path19 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r19.xlsx" # 2022-07-31 experiment
path20 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r20.xlsx" # 2022-08-1 experiment
path21 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r21.xlsx" # 2022-08-2 experiment
path22 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r22.xlsx" # 2022-08-3 experiment
path23 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r23.xlsx" # 2022-08-4 experiment
path24 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r24.xlsx" # 2022-08-5 experiment
path25 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r25.xlsx" # 2022-08-6 experiment
path26 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r26.xlsx" # 2022-08-7 experiment
path27 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r27.xlsx" # 2022-08-8 experiment
path28 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r28.xlsx" # 2022-08-9 experiment
path29 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r29.xlsx" # 2022-08-10 experiment
path30 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r30.xlsx" # 2022-08-11 experiment
path31 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r31.xlsx" # 2022-08-12 experiment
path32 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r32.xlsx" # 2022-08-13 experiment
path33 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r33.xlsx" # 2022-08-14 experiment
path34 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r34.xlsx" # 2022-08-15 experiment
path35 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r35.xlsx" # 2022-08-17 experiment
path36 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r36.xlsx" # 2022-08-18 experiment
path37 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r37.xlsx" # 2022-08-19 experiment
path38 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r38.xlsx" # 2022-08-20 experiment
path39 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r39.xlsx" # 2022-08-21 experiment
path40 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r40.xlsx" # 2022-08-22 experiment
path41 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r41.xlsx" # 2022-08-23 experiment
path42 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r42.xlsx" # 2022-08-25 experiment
path43 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r43.xlsx" # 2022-08-26 experiment
path44 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r44.xlsx" # 2022-08-28 experiment
path46 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r46.xlsx" # 2022-08-30 experiment
path47 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r47.xlsx" # 2022-08-31 experiment
path48 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r48.xlsx" # 2022-09-2 experiment
path49 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r49.xlsx" # 2022-09-4 experiment
path50 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r50.xlsx" # 2022-09-8 experiment
path51 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r51.xlsx" # 2022-09-9 experiment
path52 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r52.xlsx" # 2022-09-12 experiment

path53 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r53.xlsx" # 2022-09-14 experiment
path53_2 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r53_2.xlsx" 

path54 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r54.xlsx" # 2022-09-15 experiment
path55 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r55.xlsx" # 2022-09-16 experiment
path56 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r56.xlsx" # 2022-09-19 experiment
path57 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r57.xlsx" # 2022-09-20 experiment
path58 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r58.xlsx" # 2022-09-21 experiment
path59 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r59.xlsx" # 2022-09-23 experiment
path60 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r60.xlsx" # 2022-09-27 experiment
path61 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r61.xlsx" # 2022-10-16 experiment
path62 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r62.xlsx" # 2022-10-17 experiment
path63 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r63.xlsx" 
path63_2 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r63_2.xlsx" 
path64 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r64.xlsx" 
path65 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r65.xlsx" 
path66 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r66.xlsx" 
path67 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r67.xlsx" 
path68 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r68.xlsx" 
path69 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r69.xlsx" 
path70 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r70.xlsx" 
path73 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r73.xlsx" 
path74 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r74.xlsx" 
path75 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r75.xlsx" 
path76 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r76.xlsx" 
path77 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r77.xlsx" 
path79 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r79.xlsx" 
path80 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r80.xlsx" 
path82 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r82.xlsx" 
path83 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r83.xlsx" 
path84 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r84.xlsx" 
path85 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r85.xlsx" 
path86 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r86.xlsx" 
path87 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r87.xlsx" 

path99 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r99.xlsx" 
path100 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r100.xlsx" 
path101 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r101.xlsx" 
path102 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r102.xlsx" 
path103 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r103.xlsx" 

path104 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r104.xlsx" 
path105 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r105.xlsx" 
path106 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r106.xlsx" 
path107 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r107.xlsx" 
path108 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r108.xlsx" 
path109 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r109.xlsx" 
path110 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r110.xlsx" 

path111 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r111.xlsx" 
path112 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r112.xlsx" 
path113 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r113.xlsx" 

path114 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r114.xlsx" 
path115 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r115.xlsx" 
path116 = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_13CO2/r116.xlsx" 


d.data_13CO2.4.control = read_excel(path4)
d.data_13CO2.9 = read_excel(path9)
d.data_13CO2.10 = read_excel(path10)
d.data_13CO2.12 = read_excel(path12)
d.data_13CO2.13 = read_excel(path13)
d.data_13CO2.14 = read_excel(path14)
d.data_13CO2.16 = read_excel(path16)
d.data_13CO2.17 = read_excel(path17)
d.data_13CO2.18 = read_excel(path18)
d.data_13CO2.19 = read_excel(path19)
d.data_13CO2.20 = read_excel(path20)
d.data_13CO2.21 = read_excel(path21)
d.data_13CO2.22 = read_excel(path22)
d.data_13CO2.23 = read_excel(path23)
d.data_13CO2.24 = read_excel(path24)
d.data_13CO2.25 = read_excel(path25)
d.data_13CO2.26 = read_excel(path26)
d.data_13CO2.27 = read_excel(path27)
d.data_13CO2.28 = read_excel(path28)
d.data_13CO2.29 = read_excel(path29)
d.data_13CO2.30 = read_excel(path30)
d.data_13CO2.31 = read_excel(path31)
d.data_13CO2.32 = read_excel(path32)
d.data_13CO2.33 = read_excel(path33)
d.data_13CO2.34 = read_excel(path34)
d.data_13CO2.35 = read_excel(path35)
d.data_13CO2.36 = read_excel(path36)
d.data_13CO2.37 = read_excel(path37)
d.data_13CO2.38 = read_excel(path38)
d.data_13CO2.39 = read_excel(path39)
d.data_13CO2.40 = read_excel(path40)
d.data_13CO2.41 = read_excel(path41)
d.data_13CO2.42 = read_excel(path42)
d.data_13CO2.43 = read_excel(path43)
d.data_13CO2.44 = read_excel(path44) 
d.data_13CO2.46 = read_excel(path46)
d.data_13CO2.47 = read_excel(path47)
d.data_13CO2.48 = read_excel(path48)
d.data_13CO2.49 = read_excel(path49)
d.data_13CO2.50 = read_excel(path50)
d.data_13CO2.51 = read_excel(path51)
d.data_13CO2.52 = read_excel(path52)
d.data_13CO2.53 = read_excel(path53)
d.data_13CO2.53_2 = read_excel(path53_2)
d.data_13CO2.54 = read_excel(path54)
d.data_13CO2.55 = read_excel(path55)
d.data_13CO2.56 = read_excel(path56)
d.data_13CO2.57 = read_excel(path57)
d.data_13CO2.58 = read_excel(path58)
d.data_13CO2.59 = read_excel(path59)
d.data_13CO2.60 = read_excel(path60)
d.data_13CO2.61 = read_excel(path61)
d.data_13CO2.62 = read_excel(path62)
d.data_13CO2.63 = read_excel(path63)
d.data_13CO2.63_2 = read_excel(path63_2)
d.data_13CO2.64 = read_excel(path64)
d.data_13CO2.65 = read_excel(path65)
d.data_13CO2.66 = read_excel(path66)
d.data_13CO2.67 = read_excel(path67)
d.data_13CO2.68 = read_excel(path68)
d.data_13CO2.69 = read_excel(path69)
d.data_13CO2.70 = read_excel(path70)
d.data_13CO2.73 = read_excel(path73)
d.data_13CO2.74 = read_excel(path74)
d.data_13CO2.75 = read_excel(path75)
d.data_13CO2.76 = read_excel(path76)
d.data_13CO2.77 = read_excel(path77)
d.data_13CO2.79 = read_excel(path79)
d.data_13CO2.80 = read_excel(path80)
d.data_13CO2.82 = read_excel(path82)
d.data_13CO2.83 = read_excel(path83)
d.data_13CO2.84 = read_excel(path84)
d.data_13CO2.85 = read_excel(path85)
d.data_13CO2.86 = read_excel(path86)
d.data_13CO2.87 = read_excel(path87)

d.data_13CO2.99 = read_excel(path99)
d.data_13CO2.100 = read_excel(path100)
d.data_13CO2.101 = read_excel(path101)
d.data_13CO2.102 = read_excel(path102)
d.data_13CO2.103 = read_excel(path103)

d.data_13CO2.104 = read_excel(path104)
d.data_13CO2.105 = read_excel(path105)
d.data_13CO2.106 = read_excel(path106)
d.data_13CO2.107 = read_excel(path107)
d.data_13CO2.108 = read_excel(path108)
d.data_13CO2.109 = read_excel(path109)
d.data_13CO2.110 = read_excel(path110)
d.data_13CO2.111 = read_excel(path111)
d.data_13CO2.112 = read_excel(path112)
d.data_13CO2.113 = read_excel(path113)

d.data_13CO2.114 = read_excel(path114)
d.data_13CO2.115 = read_excel(path115)
d.data_13CO2.116 = read_excel(path116)

# Other parameters, cited globally in project
VPDB = 0.01123720



# Define function: select only equilibrated points, and remove spike point (outliers)
func.equilibrate.outLierRemoved = function(myData){
  #myData = d.data_13CO2.8
  d.data_13CO2 = myData %>% 
    rename(Time = DateTime, cage = RespCage_A) %>% 
    mutate(cage = as.double(cage),
           time.h = as.numeric(Time - min(Time)) / 3600  ) # time since the start of data recording
  
  
  # 1)
  # For each individual cage, during early 1 hour (within which time tracer would have been administered to all mice), if the delta 13C reading is 3 or 4 standard deviation above background, this indicates the start of tracer administration
  # once the tracer is administered, the delta 13C signal would rise way above ambient air baseline within one minute
  
  # define baseline threshold
  # background.13C.mean = (d.data_13CO2 %>% filter(cage %in% c(0)))$SI_13C_A %>% mean()
  # background.13C.sd = (d.data_13CO2 %>% filter(cage %in% c(0)))$SI_13C_A %>% sd()
  # delta.13C.threshold.startTime = background.13C.mean + background.13C.sd * (0)
  # 
  # # remove signals of the same level of baseline (within severl SD range) in the first 1 hour
  # d.data_13CO2 = d.data_13CO2 %>% 
  #   group_by(cage) %>% # each cage has different starting time point, as the tracer was administered one mouse after another
  #   mutate(higher.baseline = (SI_13C_A < delta.13C.threshold.startTime) & (time.h < 1) & (cage !=0)  ) %>% 
  #   filter(higher.baseline == F) %>% 
  #   select(-higher.baseline)
  
  # 2)
  # Collect only fully equilibrated data points 
  ## Clean data; this clean-up section can be safely re-run for dataset update after removal # trial and error
  index.rowsCollected = c() # Collect number of rows
  SectionSartMarker = 0
  cages = d.data_13CO2$cage
  for(i in 1: (nrow(d.data_13CO2) - 1 ) ){
    if (cages[i] != cages[i+1]){
      
      #index.rowsCollected = append(index.rowsCollected, c((i-1):(i)) )
      index.rowsCollected = append(index.rowsCollected, c((i)) )  # MANUALLY adjust # of points to include at then end of each dwell time
      # Collect the last # of points during each dwell time. e.g., i-0 : i to collect the last data point at the end of each dwell time
      # It is realized that occasionally a dwell time period does not contain expected number of data points due to instrument error. 
      # Collecting the last XX points for each defined dwell time is robust to this error. The last XX points are visually examined to reach stable state in the above timeline offset procedure
      
      SectionSartMarker = i+1
    }
  }
  d.CO2.equilibrated = d.data_13CO2[index.rowsCollected, ] %>% as_tibble()
  
  
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
d4 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.4.control) %>% mutate(cage = (cage + 3 * 9),  round = 4)
d9 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.9) %>% mutate(cage = (cage + 8 * 9),  round = 9)
d10 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.10) %>% mutate(cage = (cage + 9 * 9),  round = 10)
d12 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.12) %>% mutate(cage = (cage + 11 * 9),  round = 12)
d13 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.13) %>% mutate(cage = (cage + 12 * 9),  round = 13)
d14 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.14) %>% mutate(cage = (cage + 13 * 9),  round = 14)
d16 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.16) %>% mutate(cage = (cage + 15 * 9),  round = 16)
d17 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.17) %>% mutate(cage = (cage + 16 * 9),  round = 17)
d18 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.18) %>% mutate(cage = (cage + 17 * 9),  round = 18)
d19 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.19) %>% mutate(cage = (cage + 18 * 9),  round = 19)
d20 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.20) %>% mutate(cage = (cage + 19 * 9),  round = 20)
d21 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.21) %>% mutate(cage = (cage + 20 * 9),  round = 21)
d22 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.22) %>% mutate(cage = (cage + 21 * 9),  round = 22)
d23 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.23) %>% mutate(cage = (cage + 22 * 9),  round = 23)
d24 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.24) %>% mutate(cage = (cage + 23 * 9),  round = 24)
d25 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.25) %>% mutate(cage = (cage + 24 * 9),  round = 25)
d26 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.26) %>% mutate(cage = (cage + 25 * 9),  round = 26)
d27 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.27) %>% mutate(cage = (cage + 26 * 9),  round = 27)
d28 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.28) %>% mutate(cage = (cage + 27 * 9),  round = 28)
d29 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.29) %>% mutate(cage = (cage + 28 * 9),  round = 29)
d30 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.30) %>% mutate(cage = (cage + 29 * 9),  round = 30)
d31 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.31) %>% mutate(cage = (cage + 30 * 9),  round = 31)
d32 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.32) %>% mutate(cage = (cage + 31 * 9),  round = 32)
d33 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.33) %>% mutate(cage = (cage + 32 * 9),  round = 33)
d34 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.34) %>% mutate(cage = (cage + 33 * 9),  round = 34)
d35 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.35) %>% mutate(cage = (cage + 34 * 9),  round = 35)
d36 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.36) %>% mutate(cage = (cage + 35 * 9),  round = 36)
d37 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.37) %>% mutate(cage = (cage + 36 * 9),  round = 37)
d38 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.38) %>% mutate(cage = (cage + 37 * 9),  round = 38)
d39 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.39) %>% mutate(cage = (cage + 38 * 9),  round = 39)
d40 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.40) %>% mutate(cage = (cage + 39 * 9),  round = 40)
d41 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.41) %>% mutate(cage = (cage + 40 * 9),  round = 41)
d42 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.42) %>% mutate(cage = (cage + 41 * 9),  round = 42)
d43 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.43) %>% mutate(cage = (cage + 42 * 9),  round = 43)
d44 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.44) %>% mutate(cage = (cage + 43 * 9),  round = 44)
d46 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.46) %>% mutate(cage = (cage + 45 * 9),  round = 46)
d47 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.47) %>% mutate(cage = (cage + 46 * 9),  round = 47)
d48 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.48) %>% mutate(cage = (cage + 47 * 9),  round = 48)
d49 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.49) %>% mutate(cage = (cage + 48 * 9),  round = 49)
d50 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.50) %>% mutate(cage = (cage + 49 * 9),  round = 50)
d51 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.51) %>% mutate(cage = (cage + 50 * 9),  round = 51)
d52 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.52) %>% mutate(cage = (cage + 51 * 9),  round = 52)
d53 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.53) %>% mutate(cage = (cage + 52 * 9),  round = 53)
d53_2 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.53_2) %>% mutate(cage = (cage + 52 * 9),  round = 53)
d54 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.54) %>% mutate(cage = (cage + 53 * 9),  round = 54)
d55 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.55) %>% mutate(cage = (cage + 54 * 9),  round = 55)
d56 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.56) %>% mutate(cage = (cage + 55 * 9),  round = 56)
d57 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.57) %>% mutate(cage = (cage + 56 * 9),  round = 57)
d58 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.58) %>% mutate(cage = (cage + 57 * 9),  round = 58)
d59 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.59) %>% mutate(cage = (cage + 58 * 9),  round = 59)
d60 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.60) %>% mutate(cage = (cage + 59 * 9),  round = 60)
d61 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.61) %>% mutate(cage = (cage + 60 * 9),  round = 61)
d62 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.62) %>% mutate(cage = (cage + 61 * 9),  round = 62)
d63 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.63) %>% mutate(cage = (cage + 62 * 9),  round = 63)
d63_2 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.63_2) %>% mutate(cage = (cage + 62 * 9),  round = 63)
d64 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.64) %>% mutate(cage = (cage + 63 * 9),  round = 64)
d65 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.65) %>% mutate(cage = (cage + 64 * 9),  round = 65)
d66 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.66) %>% mutate(cage = (cage + 65 * 9),  round = 66)
d67 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.67) %>% mutate(cage = (cage + 66 * 9),  round = 67)
d68 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.68) %>% mutate(cage = (cage + 67 * 9),  round = 68)
d69 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.69) %>% mutate(cage = (cage + 68 * 9),  round = 69)
d70 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.70) %>% mutate(cage = (cage + 69 * 9),  round = 70)
d73 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.73) %>% mutate(cage = (cage + 72 * 9),  round = 73)
d74 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.74) %>% mutate(cage = (cage + 73 * 9),  round = 74)
d75 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.75) %>% mutate(cage = (cage + 74 * 9),  round = 75)
d76 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.76) %>% mutate(cage = (cage + 75 * 9),  round = 76)
d77 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.77) %>% mutate(cage = (cage + 76 * 9),  round = 77)
d79 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.79) %>% mutate(cage = (cage + 78 * 9),  round = 79)
d80 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.80) %>% mutate(cage = (cage + 79 * 9),  round = 80)
d82 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.82) %>% mutate(cage = (cage + 81 * 9),  round = 82)
d83 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.83) %>% mutate(cage = (cage + 82 * 9),  round = 83)
d84 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.84) %>% mutate(cage = (cage + 83 * 9),  round = 84)
d85 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.85) %>% mutate(cage = (cage + 84 * 9),  round = 85)
d86 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.86) %>% mutate(cage = (cage + 85 * 9),  round = 86)
d87 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.87) %>% mutate(cage = (cage + 86 * 9),  round = 87)

d99 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.99) %>% mutate(cage = (cage + 98 * 9),  round = 99)
d100 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.100) %>% mutate(cage = (cage + 99 * 9),  round = 100)
d101 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.101) %>% mutate(cage = (cage + 100 * 9),  round = 101)
d102 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.102) %>% mutate(cage = (cage + 101 * 9),  round = 102)
d103 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.103) %>% mutate(cage = (cage + 102 * 9),  round = 103)

d104 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.104) %>% mutate(cage = (cage + 103 * 9),  round = 104)
d105 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.105) %>% mutate(cage = (cage + 104 * 9),  round = 105)
d106 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.106) %>% mutate(cage = (cage + 105 * 9),  round = 106)
d107 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.107) %>% mutate(cage = (cage + 106 * 9),  round = 107)
d108 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.108) %>% mutate(cage = (cage + 107 * 9),  round = 108)
d109 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.109) %>% mutate(cage = (cage + 108 * 9),  round = 109)
d110 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.110) %>% mutate(cage = (cage + 109 * 9),  round = 110)
d111 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.111) %>% mutate(cage = (cage + 110 * 9),  round = 111)
d112 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.112) %>% mutate(cage = (cage + 111 * 9),  round = 112)
d113 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.113) %>% mutate(cage = (cage + 112 * 9),  round = 113)

d114 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.114) %>% mutate(cage = (cage + 113 * 9),  round = 114)
d115 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.115) %>% mutate(cage = (cage + 114 * 9),  round = 115)
d116 = func.equilibrate.outLierRemoved(myData = d.data_13CO2.116) %>% mutate(cage = (cage + 115 * 9),  round = 116)



# Merge all files
ds = list(d4, d9, d10, 
          d12, d13, d14, d16, d17, d18, d19, d20, 
          d21, d22, d23, d24, d25, d26, d27, d28, d29, d30, 
          d31, d32, d32, d33, d34, d35, d36, d37, d38, d39, d40,
          d41, d43, d44, d46, d47, d48, d49, d50,  
          d51, d52, d53, d53_2, d54, d55, d56, d57, d58, d59, d60, 
          d61, d62, d63, d63_2, d64, d65, d66, d67, d68, d69, d70,
          d73, d74, d75, d76, d77, d79, d80,
          d82, d83, d84, d85, d86, d87,
          d99, d100, d101, d102, d103,
          d104, d105, d106, d107, d108, 
          # Xmas efforts
          d109, d110, d111, d112, d113, d114, d115, d116) 


d.all.raw = tibble()
for(i in 1: length(ds)){ d.all.raw = rbind(d.all.raw, ds[[i]] ) }



# Combine with mouse ID dataset
d.all.raw2 = d.all.raw %>% 
  left_join(d.mouseID, by = c("cage", "round")) %>%
  mutate(round.cage = str_c(round, "_", cage)) %>% 
  # combine with integration cutoff time
  left_join(d.integrate.cutOff, by = "tracer") %>% 
  mutate(tracer = factor(tracer, levels = ordered.Compound, ordered = T))



# save project
save.image(file = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/1_core_13CO2_import_data.RData")

