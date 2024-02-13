rm(list = ls())

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



# Load my mini R package vivoFlux
setwd("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/vivoFlux")
devtools::load_all()
set.seed(2012)


theme.mybw <- theme_bw() + 
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 13, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = .5, face = "bold", size = 14),
        axis.title = element_text(size = 12, face = "bold"),
        panel.spacing = unit(10, "mm"))

theme.myClassic <-  theme_classic() +  
  theme(strip.background = element_blank(),
        axis.text = element_text(colour = "black", size=13),
        axis.title = element_text(colour = "black", size = 13, face = "bold"),
        strip.text = element_text(colour = "black", size = 12, face = "bold"),
        legend.text = element_text(colour = "black", size = 12),
        plot.title = element_text(size = 12, hjust = .5, face = "bold"),
        legend.title = element_blank(),
        panel.spacing = unit(10, "mm"))

ordered.phenotype = factor(c("WT", "HFD", "ob/ob",  "db/db"), ordered = T)
color.phenotype = c("snow4", "chocolate1", "deepskyblue3",  "firebrick2")
color.phenotype %>% show_col()
names(color.phenotype) = ordered.phenotype


ordered.Compound = c(
  "Glucose", "Lactate",  "Glycerol","Alanine", "Glutamine", "3-HB", "Acetate", 
  "C16:0", "C18:1", "C18:2", "Valine",
  "storage", "CO2", "non-Ox sink")

# ordered.Compound.abbreviated = ordered.Compound %>% str_replace("3-HB", "3-HB")

ordered.gluconeogenic.Set = ordered.Compound[1:5]

color.Compounds = c(pal_nejm()(7), 
                    "seashell", "papayawhip",  "peachpuff", "purple3",
                    "grey", "azure",  "slategrey")

color.Compounds[1] <- "tomato"
color.Compounds[2] <- "turquoise3"
color.Compounds[3] <- "khaki1"
color.Compounds[4] <- "olivedrab4" # "lightseagreen"
color.Compounds[5] <- "mistyrose"

color.Compounds %>% show_col()
# pal_nejm()(7) %>% show_col()

# c("#0072B5FF", "steelblue") %>% show_col()
# c("#0072B5FF", "steelblue") %>% show_col()

names(color.Compounds) = ordered.Compound







# Mouse infusion overview
ph = c("WT", "ob/ob", "HFD", "db/db")
names(ph) = c("L", "O", "H", "d")




# Basic mouse and sample documentation
path.infusion.surgery = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_serum labeling mouse ID.xlsx"
d.mouse = read_excel(path.infusion.surgery, sheet = "mouse data")
d.infusion_rounds = read_excel(path = path.infusion.surgery, sheet = "infusion rounds",
                               # skip = 7, 
                               # range = "A8:O492",
                               range = "A8:F525")
d.standardInfusionParameters = read_excel(path.infusion.surgery, sheet = "infusion parameters")


# Labeling data =======
path.labeling = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/data_serum labeling.xlsm"

d.LCMS_sampleInfo = read_excel(path.labeling, sheet = "sample_info") 
# Check there is no duplicated rows
if (d.LCMS_sampleInfo %>% duplicated() %>% sum() >0) stop("There is duplicated rows in the mouse_ID dataset.")

d = read_excel(path.labeling, sheet = "a_h")
d = d %>% select( - contains("d-")) %>% # remove infusion round d; likely cross contamination of artery blood by tracer residual on the button; ca. 80% labeling in art. blood
  select( - contains( "a-tail-d1")) %>% # mouse jumped out of cage and fell from table to floor during infusion; ca. only 5% M+6 glucose labeling, aberrantly low; 
  select(-c(Q1,Q2, Q3, Q4)) # remove quality control samples

# background subtraction & natural abundance correction
d.backgroundSubtracted = d %>% flx.subtract_background(blankSamples = c("blank1", "blank2", "blank3", "blank4", "blank_2ndSeq"))
list.corrected = d.backgroundSubtracted %>% flx.correct_natural_abundance()

# combine infusion and sample ID data
func.combine_sample_ID = function(dataset.tidy){
  d.tidy = dataset.tidy %>% 
    mutate(mouse_ID = str_extract(sample, pattern = WRD %R% one_or_more(DIGIT) %R% END), # L/O/d, with mouse number within phenotype
           infusion_round = str_extract(sample, START %R% one_or_more(WRD) )) %>%
    left_join(d.LCMS_sampleInfo, by = c("sample", "infusion_round", "mouse_ID")) %>% 
    left_join(d.infusion_rounds, by = c("infusion_round", "mouse_ID")) 
  return(d.tidy)
} 

d.normalized.tidy.1 = func.combine_sample_ID(
  dataset.tidy = list.corrected$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )

d.corrected.tidy.1 = func.combine_sample_ID(
  dataset = list.corrected$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  # the intensity will be normalized separately for each separate MS.run 
  mutate(MS.run = "1_a_h")





# Sheet 2: RU analysis in Feb 2022
list.corrected.2 = read_excel(path.labeling, sheet = "RU_x-aa") %>% 
  select(-c("z-tail-O1","aa-art-O3", "aa-tail-O3")) %>%  # glucose labeling only 3~4%; signifiant mice age by time of infusion, button leaking in O1 in Dec 2021 - Feb 2022 glucose clamp study; 
  flx.subtract_background(blankSamples =  c("Blank1",	"Blank3",	"Blank4")) %>% 
  natural_abundance_correction(resolution = 70000)

d.normalized.tidy.2 = func.combine_sample_ID(
  dataset.tidy = list.corrected.2$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )

d.corrected.tidy.2 = func.combine_sample_ID(
  dataset = list.corrected.2$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  mutate(MS.run = "2_RU_x-aa")





# Sheet 3: RU analysis in March 2022, analyzed by Rutgers with resolution of 70,000
list.corrected.3 = read_excel(path.labeling, sheet = "RU_i-n") %>% 
  flx.subtract_background(blankSamples =  paste0("Blank", 3:9)) %>% 
  natural_abundance_correction(resolution = 70000)

d.normalized.tidy.3 = func.combine_sample_ID(
  dataset.tidy = list.corrected.3$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )

d.corrected.tidy.3 = func.combine_sample_ID(
  dataset = list.corrected.3$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  mutate(MS.run = "3_RU_i-n")




# Sheet 4: clamp data, basal & clamp
d4 = read_excel(path.labeling, sheet = "clamp_af-ah") %>% select(!contains("-2"))  # remove -2 (being the enrichment at the end of the clamp)

# rename file names to be consistent with the general naming rule
d4.1 = d4[, 1:3] # first three columns
d4.2 = d4[, 4:ncol(d4)]
d4.rounds = colnames(d4.2) %>% str_extract(pattern = START %R% one_or_more(WRD))
d4.mouse =  colnames(d4.2) %>% str_extract(pattern = "-" %R% one_or_more(WRD) %R% one_or_more(DGT) ) %>% str_remove("-")
d4.2.filenames = d4.rounds %>% str_c("-tail-", d4.mouse)
d4.2.filenames[28:32] = paste0("blank", 1:5)
colnames(d4.2) = d4.2.filenames
d4 = cbind(d4.1, d4.2) %>% as_tibble()

list.corrected.4 = d4 %>% 
  flx.subtract_background(blankSamples =  paste0("blank", 1:5)) %>% 
  natural_abundance_correction(resolution = 12000)

d.normalized.tidy.4 = func.combine_sample_ID(
  dataset.tidy = list.corrected.4$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )

d.corrected.tidy.4 = func.combine_sample_ID(
  dataset = list.corrected.4$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  mutate(MS.run = "4_clamp_af-ah")




# Sheet 5: Glycerol-3P
list.corrected.5 = read_excel(path.labeling, sheet = "G3P_o-p") %>% 
  
  # remove outlier samples
  select(-`p-tail-d5`) %>% 
  
  flx.subtract_background(blankSamples =  paste0("blank_glycerol3P_", 1:3)) %>% 
  natural_abundance_correction(resolution = 12000)

d.normalized.tidy.5 = func.combine_sample_ID(
  dataset.tidy = list.corrected.5$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )
d.corrected.tidy.5 = func.combine_sample_ID(
  dataset = list.corrected.5$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  mutate(MS.run = "5_G3P_o-p")

# remove artery sample due to contamination from lock solution; so much higher labeling in venous blood than artery blood due to dilution by lock solution-glycerol
d.normalized.tidy.5 = d.normalized.tidy.5 %>% filter(sample %>% str_detect("tail"))
d.corrected.tidy.5 = d.corrected.tidy.5 %>% filter(sample %>% str_detect("tail"))




# Sheet 6: r to w
list.corrected.6 = read_excel(path.labeling, sheet = "o-p-r-s-t-u-v-w") %>% 
  flx.subtract_background(blankSamples =  paste0("blank_", 6:9)) %>% 
  natural_abundance_correction(resolution = 12000)

d.normalized.tidy.6 = func.combine_sample_ID(
  dataset.tidy = list.corrected.6$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )

d.corrected.tidy.6 = func.combine_sample_ID(
  dataset = list.corrected.6$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  mutate(MS.run = "6_o-p-r-s-t-u-v-w")




# Sheet 7: glycerol 3-P, rounds # a to n;  q to z
list.corrected.7 = read_excel(path.labeling, sheet = "G3P_a-n_q-z-aa-ae") %>% 
  flx.subtract_background(blankSamples =  c("blank2", "blank4")) %>% 
  natural_abundance_correction(resolution = 12000)

d.normalized.tidy.7 = func.combine_sample_ID(
  dataset.tidy = list.corrected.7$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )

d.corrected.tidy.7 = func.combine_sample_ID(
  dataset = list.corrected.7$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  mutate(MS.run = "7_G3P_a-n_q-z-aa-ae")


# Sheet 8: gluconeogenic, ab, ac, ad
list.corrected.8 = read_excel(path.labeling, sheet = "ab_ac_ad") %>% 
  flx.subtract_background(blankSamples =  paste0("blank", 1:6) ) %>% 
  natural_abundance_correction(resolution = 12000)

d.normalized.tidy.8 = func.combine_sample_ID(
  dataset.tidy = list.corrected.8$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )

d.corrected.tidy.8 = func.combine_sample_ID(
  dataset = list.corrected.8$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  mutate(MS.run = "8_ab_ac_ad")



# Sheet 9:  glycerol 3-P & FA & HILIC metabolites, rounds # am to ao
list.corrected.9 = read_excel(path.labeling, sheet = "G3P_HILIC_am-ao") %>% 
  flx.subtract_background(blankSamples =  paste0("blk", 1:3) ) %>% 
  natural_abundance_correction(resolution = 12000)

d.normalized.tidy.9 = func.combine_sample_ID(
  dataset.tidy = list.corrected.9$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )

d.corrected.tidy.9 = func.combine_sample_ID(
  dataset = list.corrected.9$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  mutate(MS.run = "9_G3P_HILIC_am-ao")




# Sheet 10:  glycerol 3-P, rounds # ar to bd
list.corrected.10 = read_excel(path.labeling, sheet = "G3P_ar_bd") %>% 
  natural_abundance_correction(resolution = 12000)

d.normalized.tidy.10 = func.combine_sample_ID(
  dataset.tidy = list.corrected.10$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )

d.corrected.tidy.10 = func.combine_sample_ID(
  dataset = list.corrected.10$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  mutate(MS.run = "10_G3P_ar_bd")


# Sheet 11:  # ar to bd, early elution 
list.corrected.11 = read_excel(path.labeling, sheet = "ar_bd_EarlyElution") %>% 
  natural_abundance_correction(resolution = 12000)

d.normalized.tidy.11 = func.combine_sample_ID(
  dataset.tidy = list.corrected.11$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )

d.corrected.tidy.11 = func.combine_sample_ID(
  dataset = list.corrected.11$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  mutate(MS.run = "11_ar_bd_EarlyElution")



# Sheet 12:  # ar to bd, late elution 
list.corrected.12 = read_excel(path.labeling, sheet = "ar_bd_LateElution") %>% 
  natural_abundance_correction(resolution = 12000)

d.normalized.tidy.12 = func.combine_sample_ID(
  dataset.tidy = list.corrected.12$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )

d.corrected.tidy.12 = func.combine_sample_ID(
  dataset = list.corrected.12$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  mutate(MS.run = "12_ar_bd_LateElution")


# Sheet 13:  # be to bf
list.corrected.13 = read_excel(path.labeling, sheet = "be-bf_valine") %>% 
  flx.subtract_background(blankSamples =  c("blk2", "blk4") ) %>% 
  natural_abundance_correction(resolution = 12000)

d.normalized.tidy.13 = func.combine_sample_ID(
  dataset.tidy = list.corrected.13$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )

d.corrected.tidy.13 = func.combine_sample_ID(
  dataset = list.corrected.13$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  mutate(MS.run = "13_be-bf_valine")


# Sheet 14:  # bg to bk
list.corrected.14 = read_excel(path.labeling, sheet = "bg-bk") %>% 
  flx.subtract_background(blankSamples =  c("blk_3", "blk_4") ) %>% 
  natural_abundance_correction(resolution = 12000)

d.normalized.tidy.14 = func.combine_sample_ID(
  dataset.tidy = list.corrected.14$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )

d.corrected.tidy.14 = func.combine_sample_ID(
  dataset = list.corrected.14$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  mutate(MS.run = "14_bg-bk")

# ----------
# sheet 15: bl - bp: C18:1 and C18:2 
d15 = read_excel(path.labeling, sheet = "bl-bp")

# remove the second rep, and rename the sample to remove -1 in the sample ID
d15 = d15 %>% select(-contains("-2"))
colnames(d15) = colnames(d15) %>% str_remove(pattern = "-1")

list.corrected.15 = d15 %>% 
  flx.subtract_background(blankSamples =  c("blk_3", "blank_1")) %>% 
  natural_abundance_correction(resolution = 12000)

d.normalized.tidy.15 = func.combine_sample_ID(
  dataset.tidy = list.corrected.15$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )

d.corrected.tidy.15 = func.combine_sample_ID(
  dataset = list.corrected.15$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  mutate(MS.run = "15_bl-bp")


# sheet 16: bq - bs
d16 = read_excel(path.labeling, sheet = "bq-bs")
d16 = d16 %>% 
  # for mice with 3 replicates of blood collection (spaced by 10 min), take the 2nd collection
  select(-contains("-3")) %>% select(-contains("-1"))
colnames(d16) = colnames(d16) %>% str_remove("-2") %>% str_replace("-", "-tail-")

list.corrected.16 = d16 %>% 
  flx.subtract_background(blankSamples =  paste0("buffer_", 1:4)) %>% 
  natural_abundance_correction(resolution = 12000)

d.normalized.tidy.16 = func.combine_sample_ID(
  dataset.tidy = list.corrected.16$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )

d.corrected.tidy.16 = func.combine_sample_ID(
  dataset = list.corrected.16$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  mutate(MS.run = "16_bq-bs")



# 
# sheet 17: G3P: be - bk
d17 = read_excel(path.labeling, sheet = "G3P_be_bk")
list.corrected.17 = d17 %>% flx.subtract_background(blankSamples = "Enz-blank") %>% 
  natural_abundance_correction(resolution = 12000)

d.normalized.tidy.17 = func.combine_sample_ID(
  dataset.tidy = list.corrected.17$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )

d.corrected.tidy.17 = func.combine_sample_ID(
  dataset = list.corrected.17$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  mutate(MS.run = "17_G3P_be_bk")



# sheet 18: bt - bw 
d18 = read_excel(path.labeling, sheet = "bt-bw")
list.corrected.18 = d18 %>% 
  flx.subtract_background(blankSamples = c("blank_3", "blank_4", "blank_5")) %>% 
  natural_abundance_correction(resolution = 12000)

d.normalized.tidy.18 = func.combine_sample_ID(
  dataset.tidy = list.corrected.18$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )

d.corrected.tidy.18 = func.combine_sample_ID(
  dataset = list.corrected.18$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  mutate(MS.run = "18_bt-bw")


# sheet 19: bz-ca
d19 = read_excel(path.labeling, sheet = "bz-ca")
list.corrected.19 = d19 %>% 
  flx.subtract_background(blankSamples = c("blank_2",	"blank_3",	"blank_4")) %>% 
  natural_abundance_correction(resolution = 12000)

d.normalized.tidy.19 = func.combine_sample_ID(
  dataset.tidy = list.corrected.19$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )

d.corrected.tidy.19 = func.combine_sample_ID(
  dataset = list.corrected.19$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  mutate(MS.run = "19_bz-ca")

# sheet 20: ce
d20 = read_excel(path.labeling, sheet = "ce")
list.corrected.20 = d20 %>% 
  flx.subtract_background(blankSamples = c("blk_2")) %>% 
  natural_abundance_correction(resolution = 12000)

d.normalized.tidy.20 = func.combine_sample_ID(
  dataset.tidy = list.corrected.20$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )

d.corrected.tidy.20 = func.combine_sample_ID(
  dataset = list.corrected.20$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  mutate(MS.run = "20_ce")


# sheet 21: cf
d21 = read_excel(path.labeling, sheet = "cf")
list.corrected.21 = d21 %>% 
  natural_abundance_correction(resolution = 12000)

d.normalized.tidy.21 = func.combine_sample_ID(
  dataset.tidy = list.corrected.21$Normalized %>% 
    gather(-c(Compound, C_Label), key = sample, value = enrichment) )

d.corrected.tidy.21 = func.combine_sample_ID(
  dataset = list.corrected.21$Corrected %>% 
    gather(-c(Compound, C_Label), key = sample, value = intensity)) %>% 
  mutate(MS.run = "21_cf")




# Compile all imported data
func.combine.tidies = function(inputList, cumulator) {
  for (i in 1:length(inputList)) {
    cumulator = rbind(cumulator, inputList[[i]])
  }
  return(cumulator)
}  


cumulator = NULL

d.normalized.tidy = func.combine.tidies(inputList = list(
  d.normalized.tidy.1,  d.normalized.tidy.2,  d.normalized.tidy.3,  d.normalized.tidy.4, 
  d.normalized.tidy.5,  d.normalized.tidy.6,  d.normalized.tidy.7,  d.normalized.tidy.8, 
  d.normalized.tidy.9,  d.normalized.tidy.10, d.normalized.tidy.11, d.normalized.tidy.12, 
  d.normalized.tidy.13, d.normalized.tidy.14, d.normalized.tidy.15, d.normalized.tidy.16,
  d.normalized.tidy.17, d.normalized.tidy.18, d.normalized.tidy.19, d.normalized.tidy.20,
  d.normalized.tidy.21), 
  cumulator = cumulator)

cumulator = NULL
d.corrected.tidy = tibble()

d.corrected.tidy = func.combine.tidies(inputList = list(
  d.corrected.tidy.1,  d.corrected.tidy.2,  d.corrected.tidy.3,  d.corrected.tidy.4, 
  d.corrected.tidy.5,  d.corrected.tidy.6,  d.corrected.tidy.7,  d.corrected.tidy.8, 
  d.corrected.tidy.9,  d.corrected.tidy.10, d.corrected.tidy.11, d.corrected.tidy.12,
  d.corrected.tidy.13, d.corrected.tidy.14, d.corrected.tidy.15, d.corrected.tidy.16,
  d.corrected.tidy.17, d.corrected.tidy.18, d.corrected.tidy.19, d.corrected.tidy.20,
  d.corrected.tidy.21),
  cumulator = d.corrected.tidy)

d.corrected.tidy$blood %>% unique()
# check duplication
d.normalized.tidy %>% duplicated() %>% sum(); d.corrected.tidy %>% duplicated() %>% sum()




# Remove outliers from the compiled dataset; select only relevant metabolites 
samples.toRemove = c(
  "am-tail-O30", "am-tail-O32", "ap-tail-L30", "ap-tail-L31",
  "ac-tail-d3", "ac-tail-d2", "ah-tail-O18", "b-tail-O3",
  "p-tail-d8", # tracer glycerol, glycerol labeling higher than usual, but no labeling in glucose 
  "an-tail-O32", # tracer C16:0, extremely low labeling 
  "O-art-O5", # tracer glycerol, there is no labeling in lactate and glucose in the artery blood
  "O-art-d7", # tracer glycerol, little labeling in alanine
  "am-tail-O16", # tracer glucose, higher infusion rate to label TG-glycerol, not ordinary infusion experiment
  "e-art-d3", "e-art-O1", "e-tail-d3", "e-tail-d2", "e-tail-O1", # tracer lactate, abnormally low labeling in glucose
  "ad-tail-O3", "ad-art-O3", "e-tail-O4", # tracer lactate: too much high labeling in tail labeling and too low Fcirc
  "bc-tail-L64", # too low labeling, WT, palmitate infusion
  "ar-tail-O59", "r-tail-O4", "s-tail-O7", "s-art-O7", # two low labeling, ob/ob, 3-HB infusion
  "be-tail-L70", # too low labeling, WT, valine infusion
  "bd-tail-L66", # too high enrichment at 60%, WT, valine infusion
  "bk-tail-L64", # Valine, WT, extremely low labeling
  # HFD, remove double cathetered mice, as the enrichment is consistently higher in serum than other genotypes and serum has unusual bright yellow corlor. 
  "bg-art-H32", "bg-art-H33", "bh-art-H36", "bh-art-H37",
  "bg-tail-H31",  "bg-tail-H32", "bg-tail-H33", "bh-tail-H35", "bh-tail-H36", "bh-tail-H37",
  "bw-tail-O105", # can draw blood, but canNOT push saline through button!
  "ad-art-L8", "ad-tail-L8" # way too high labeling of lactate
)


# fatty acid outliers 
outlier.C16.0 <- c("bs-tail-L93", "bs-tail-L94", "bs-tail-L92", 
                   "br-tail-H53", "bq-tail-H51", "bf-tail-H19", "br-tail-H56", "br-tail-H57",
                   "br-tail-H51", "bq-tail-H56",  "bf-tail-H8", "bf-tail-H3")

# outlier.C18.1 <- c("bn-tail-H53", "bp-tail-H51", "bm-tail-H14", "br-tail-O84",
#                    "br-tail-O86", "br-tail-O83", "bp-tail-O85",
#                    
#                    "au-tail-L66", "bo-tail-L93", "aw-tail-L63", "bm-tail-L82", "bm-tail-H7",
#                    "bl-tail-H52", "au-tail-O61", "au-tail-O53",
#                    
#                    "aw-tail-O56", "br-tail-O85", "aw-tail-O60", "bi-tail-O82")
# 
# outlier.C18.2 <- c("ay-tail-L64", "ax-tail-L66", "bk-tail-L82", "bk-tail-H53",
#                    "bk-tail-O82", "ax-tail-O61", "ax-tail-O57", "ay-tail-O60")

outlier.fattyacids <- c(outlier.C16.0)

samples.toRemove <- c(samples.toRemove, outlier.fattyacids)

func.RemoveOutlierSample = function(dataset, outliers){
  dataset %>% filter(! sample %in% outliers ) %>% filter(Compound %in% ordered.Compound) %>%   return()
}

d.corrected.tidy = func.RemoveOutlierSample(dataset = d.corrected.tidy, outliers = samples.toRemove)
d.normalized.tidy = func.RemoveOutlierSample(dataset = d.normalized.tidy, outliers = samples.toRemove)



# keep C18:2 data only in rounds bt, bu, and bw, where lock solution was used, and blood can be drawn from the cathether
d.corrected.tidy <- d.corrected.tidy %>% 
  filter(!(infused.tracer == "C18:2" & !infusion_round %in% c("bt", "bu", "bw") ))

d.normalized.tidy <- d.normalized.tidy %>% 
  filter(!(infused.tracer == "C18:2" & !infusion_round %in% c("bt", "bu", "bw") ))


# keep C18:1 K (the free package) data only in rounds ce, cf, where lock solution was used, and blood can be drawn from the cathether 
d.corrected.tidy <- d.corrected.tidy %>% 
  filter(!(infused.tracer == "C18:1" & !infusion_round %in% c("ce", "cf") ))

d.normalized.tidy <- d.normalized.tidy %>% 
  filter(!(infused.tracer == "C18:1" & !infusion_round %in% c("ce", "cf") ))



# Remove the early experiments of C16:0 infusion where tracer enrichment is abnormally low
d.corrected.tidy = d.corrected.tidy %>% 
  filter(! ( (infused.tracer == "C16:0") & (infusion_round %in% c("q", "t")) ))

d.normalized.tidy = d.normalized.tidy %>% 
  filter(! ( (infused.tracer == "C16:0") & (infusion_round %in% c("q", "t")) ))



save.image(file = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/4_core_serum_labeling_import data.RData")

