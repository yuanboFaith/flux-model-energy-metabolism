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

load("/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/7_core_production_flux.RData")

options(dplyr.print_max = 4)

# Extrapolate from valine OX and NOS to total protein utilization

# import amino acid basics
path.aa <- "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/amino acid basics.xlsx"
d.aa <- read_excel(path.aa, sheet = "AA.basics") 
d.aa$C.number <- d.aa$C.number %>% as.double()

# define function to count amino acid frequency in a protein
func.AA.count <- function(aminoAcid.sequence, whichProtein){
  # split
  seq.splitted = str_split(aminoAcid.sequence, pattern = "")[[1]]
  # count
  seq.count = table(seq.splitted) / length(seq.splitted) * 100
  # organize into a data frame
  t <- tibble(Abbreviation_2 = names(seq.count),
              percent = c(seq.count)) 
  colnames(t)[2] <- whichProtein
  
  return(t)
}

# mouse myosin 1: Q5SX40 · MYH1_MOUSE : https://www.uniprot.org/uniprotkb/Q5SX40/entry#sequences
s1 <- "MSSDAEMAVFGEAAPYLRKSEKERIEAQNKPFDAKSSVFVVDAKESFVKATVQSREGGKVTAKTEGGTTVTVKDDQVYPMNPPKYDKIEDMAMMTHLHEPAVLYNLKERYAAWMIYTYSGLFCVTVNPYKWLPVYNAEVVAAYRGKKRQEAPPHIFSISDNAYQFMLTDRENQSILITGESGAGKTVNTKRVIQYFATIAVTGEKKKEEATSGKMQGTLEDQIISANPLLEAFGNAKTVRNDNSSRFGKFIRIHFGTTGKLASADIETYLLEKSRVTFQLKAERSYHIFYQIMSNKKPDLIEMLLITTNPYDYAFVSQGEITVPSIDDQEELMATDSAIDILGFTSDERVSIYKLTGAVMHYGNMKFKQKQREEQAEPDGTEVADKAAYLQNLNSADLLKALCYPRVKVGNEYVTKGQTVQQVYNSVGALAKAVYEKMFLWMVTRINQQLDTKQPRQYFIGVLDIAGFEIFDFNSLEQLCINFTNEKLQQFFNHHMFVLEQEEYKKEGIEWEFIDFGMDLAACIELIEKPMGIFSILEEECMFPKATDTSFKNKLYEQHLGKSNNFQKPKPAKGKVEAHFSLVHYAGTVDYNIAGWLDKNKDPLNETVVGLYQKSSMKTLAYLFSGAAAAAEAESGGGGGKKGAKKKGSSFQTVSALFRENLNKLMTNLRSTHPHFVRCIIPNETKTPGAMEHELVLHQLRCNGVLEGIRICRKGFPSRILYADFKQRYKVLNASAIPEGQFIDSKKASEKLLGSIDIDHTQYKFGHTKVFFKAGLLGLLEEMRDDKLAQLITRTQAMCRGYLARVEYQKMVERRESIFCIQYNVRAFMNVKHWPWMKLYFKIKPLLKSAETEKEMANMKEEFEKAKENLAKAEAKRKELEEKMVALMQEKNDLQLQVQSEADSLADAEERCDQLIKTKIQLEAKIKEVTERAEDEEEINAELTAKKRKLEDECSELKKDIDDLELTLAKVEKEKHATENKVKNLTEEMAGLDETIAKLTKEKKALQEAHQQTLDDLQAEEDKVNTLTKAKIKLEQQVDDLEGSLEQEKKIRMDLERAKRKLEGDLKLAQESTMDVENDKQQLDEKLKKKEFEMSNLQSKIEDEQALGMQLQKKIKELQARIEELEEEIEAERASRAKAEKQRSDLSRELEEISERLEEAGGATSAQIEMNKKREAEFQKMRRDLEEATLQHEATAATLRKKHADSVAELGEQIDNLQRVKQKLEKEKSEMKMEIDDLASNMEVISKSKGNLEKMCRTLEDQVSELKTKEEEQQRLINELTAQRGRLQTESGEYSRQLDEKDSLVSQLSRGKQAFTQQIEELKRQLEEEIKAKSALAHALQSSRHDCDLLREQYEEEQEAKAELQRAMSKANSEVAQWRTKYETDAIQRTEELEEAKKKLAQRLQDAEEHVEAVNAKCASLEKTKQRLQNEVEDLMIDVERTNAACAALDKKQRNFDKILAEWKQKYEETHAELEASQKESRSLSTELFKIKNAYEESLDHLETLKRENKNLQQEISDLTEQIAEGGKRIHELEKIKKQIEQEKSELQAALEEAEASLEHEEGKILRIQLELNQVKSEIDRKIAEKDEEIDQLKRNHIRVVESMQSTLDAEIRSRNDAIRLKKKMEGDLNEMEIQLNHSNRMAAEALRNYRNTQGILKDTQLHLDDALRGQEDLKEQLAMVERRANLLQAEIEELRATLEQTERSRKIAEQELLDASERVQLLHTQNTSLINTKKKLETDISQIQGEMEDIVQEARNAEEKAKKAITDAAMMAEELKKEQDTSAHLERMKKNLEQTVKDLQHRLDEAEQLALKGGKKQIQKLEARVRELEGEVENEQKRNVEAIKGLRKHERRVKELTYQTEEDRKNVLRLQDLVDKLQSKVKAYKRQAEEAEEQSNVNLAKFRKIQHELEEAEERADIAESQVNKLRVKSREVHTKIISEE"
s1 <- func.AA.count(aminoAcid.sequence = s1, whichProtein = "myosin I")


# mouse actin in sckeletal muscle : P68134 · ACTS_MOUSE : https://www.uniprot.org/uniprotkb/P68134/entry#sequences
s2 <- "MCDEDETTALVCDNGSGLVKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIITNWDDMEKIWHHTFYNELRVAPEEHPTLLTEAPLNPKANREKMTQIMFETFNVPAMYVAIQAVLSLYASGRTTGIVLDSGDGVTHNVPIYEGYALPHAIMRLDLAGRDLTDYLMKILTERGYSFVTTAEREIVRDIKEKLCYVALDFENEMATAASSSSLEKSYELPDGQVITIGNERFRCPETLFQPSFIGMESAGIHETTYNSIMKCDIDIRKDLYANNVMSGGTTMYPGIADRMQKEITALAPSTMKIKIIAPPERKYSVWIGGSILASLSTFQQMWITKQEYDEAGPSIVHRKCF"
s2 <- func.AA.count(aminoAcid.sequence = s2, whichProtein = "actin")

# collagen : Q63870 · CO7A1_MOUSE : https://www.uniprot.org/uniprotkb/Q63870/entry#sequences
s3 <- "MRLRLLVAALCAAEILMGAPEVWAQPRDRVTCTRLYAADIVFLLDGSSSIGRSNFREVRGFLEGLVLPFSGAASAQGVRFATVQYSDDPQTEFGLDTLGSGSDTIRAIRELSYKGGNTRTGAALHHVSDRVFLPRLTRPGVPKVCILITDGKSQDLVDTAAQKLKGQGVKLFAVGIKNADPEELKRVASQPTSDFFFFVNDFSILRTLLPLISRRVCTTAGGVPVTLPSDDTPSGPRDLVLSEPSSQSLRVQWTAASGPVTGYKVQYTPLTGLGQPLPSERQEVNIPAGETSTRLQGLRPLTDYQVTVVALYANSIGEAVSGTARTTAKEGLELSLQNITSHSLLVAWRRVPGANGYRVTWRDLSGGPTQQQDLSPGQGSVFLDHLEPGTDYEVTVSALFGHSVGPAASLTARTASSVEQTLHPIILSPTSILLSWNLVPEARGYRLEWRRESGLETPQKVELPPDVTRHQLDGLQPGTEYRLTLYTLLEGREVATPATVVPTGLEQLVSPVMNLQAIELPGQRVRVSWNPVPGATEYRFTVRTTQGVERTLLLPGSQTTFDLDDVRAGLSYTVRVSARVGAQEGDASILTIHRDPEAPLVVPGLRVVASDATRIRVAWGLVPGASGFRISWRTGSGPESSRTLTPDSTVTDILGLQPSTSYQVAVSALRGREEGPPVVIVARTDPLGPVRRVHLTQAGSSSVSITWTGVPGATGYRVSWHSGHGPEKSLLVSGDATVAEIDGLEPDTEYIVRVRTHVAGVDGAPASVVVRTAPEPVGSVSKLQILNASSDVLRVTWVGVPGATSYKLAWGRSEGGPMKHRILPGNKESAEIRDLEGGVSYSVRVTALVGDREGAPVSIVITTPPATPALLETLQVVQSGEHSLRLRWEPVPGAPGFRLHWQPEGGQEQSLTLGPESNSYNLVGLEPATKYQVWLTVLGQTGEGPPRKVTAYTEPSHIPSTELRVVDTSIDSVTLTWTPVSGASSYILSWRPLRGTGQEVPRAPQTLPGTSSSHRVTGLEPGISYVFSLTPIQSGVRGSEISVTQTPACSHGPVDVVFLLHATRDNAHNAEAVRRVLERLVSALGPLGPQAAQVGLLTYSHRPSPLFPLNSSHDLGIILRKIRDIPYVDPSGNNLGTAVTTAHRYLLASNAPGRRQQVPGVMVLLVDEPLRGDILSPIREAQTSGLKVMALSLVGADPEQLRRLAPGTDPIQNFFAVDNGPGLDRAVSDLAVALCQAAVTIEPQTGPCAVHCPKGQKGEPGVTGLQGQAGPPGPPGLPGRTGAPGPQGPPGSTQAKGERGFPGPEGPPGSPGLPGVPGSPGIKGSTGRPGPRGEQGERGPQGPKGEPGEPGQITGGGGPGFPGKKGDPGPSGPPGSRGPVGDPGPRGPPGLPGISVKGDKGDRGERGPPGPGIGASEQGDPGLPGLPGSPGPQGPAGRPGEKGEKGDCEDGGPGLPGQPGPPGEPGLRGAPGMTGPKGDRGLTGTPGEPGVKGERGHPGPVGPQGLPGAAGHPGVEGPEGPPGPTGRRGEKGEPGRPGDPAVGPGGAGAKGEKGEAGLPGPRGASGSKGEQGAPGLALPGDPGPKGDPGDRGPIGLTGRAGPTGDSGPPGEKGEPGRPGSPGPVGPRGRDGEAGEKGDEGIPGEPGLPGKAGERGLRGAPGPRGPVGEKGDQGDPGEDGRNGSPGSSGPKGDRGEPGPPGPPGRLVDAGIESRDKGEPGQEGPRGPKGDPGPPGVSGERGIDGLRGPPGPQGDPGVRGPAGDKGDRGPPGLDGRSGLDGKPGAPGPPGLHGASGKAGDPGRDGLPGLRGEHGPPGPPGPPGVPGKAGDDGKPGLNGKNGDPGDPGEDGRKGEKGDSGAPGREGPDGPKGERGAPGNPGLQGPPGLPGQVGPPGQGFPGVPGITGPKGDRGETGSKGEQGLPGERGLRGEPGSLPNAERLLETAGIKVSALREIVDTWDESSGSFLPVPERRPGPKGDPGDRGPPGKEGLIGFPGERGLKGERGDPGPQGPPGLALGERGPPGPPGLAGEPGKPGIPGLPGRAGGSGEAGRPGERGERGEKGERGDQGRDGLPGLPGPPGPPGPKVAIEEPGPGLAREQGPPGLKGAKGEPGSDGDPGPKGDRGVPGIKGDVGEPGKRGHDGNPGLPGERGVAGPEGKPGLQGPRGTPGPVGSHGDPGPPGAPGLAGPAGPQGPSGLKGEPGETGPPGRGLPGPVGAVGLPGPPGPSGLVGPQGSPGLPGQVGETGKPGPPGRDGSSGKDGDRGSPGVPGSPGLPGPVGPKGEPGPVGAPGQVVVGPPGAKGEKGAPGDLAGALLGEPGAKGDRGLPGPRGEKGEAGRAGGPGDPGEDGQKGAPGLKGLKGEPGIGVQGPPGPSGPPGMKGDLGPPGAPGAPGVVGFPGQTGPRGETGQPGPVGERGLAGPPGREGAPGPLGPPGPPGSAGAPGASGLKGDKGDPGAGLPGPRGERGEPGVRGEDGHPGQEGPRGLVGPPGSRGEQGEKGAAGAAGLKGDKGDSAVIEGPPGPRGAKGDMGERGPRGIDGDKGPRGESGNPGDKGSKGEPGDKGSAGSIGVRGLTGPKGEPGAAGIPGEPGAPGKDGIPGFRGDKGDIGFMGPRGLKGEKGIKGTCGRDGERGDKGEAGFPGRPGLAGKKGDMGEPGLPGQSGAPGKEGLIGPKGDRGFDGQSGPKGDQGEKGERGPPGVGGFPGPRGNDGSSGPPGPPGGVGPKGPEGLQGQKGERGPPGESVVGAPGAPGTPGERGEQGRPGPAGPRGEKGEAALTEDDIRDFVRQEMSQHCACQGQFIASGSRPLPGYAADTAGSQLHHVPVLRVSHVEEEGQVPPEDDDDFSEYSVYSVEDYQEPEVPWDGEAEIKGWDQRGSDLCSLPLDEGSCTAYTLRWYHRAVPGGTACHPFVYGGCGGNANRFGTREACERRCPPQGVHSQKTGAA"
s3 <- func.AA.count(aminoAcid.sequence = s3, whichProtein = "collagen")

s <- left_join(s1, s2) %>% left_join(s3)

# combine with amino acid basics dataset
s.all <- s %>% left_join(d.aa, by = "Abbreviation_2") %>% 
  select(aminoAcids, C.number, Abbreviation_1, `myosin I`, actin, collagen)

# tidy and calculate carbon percent
s.all.carbon.pct <- s.all %>% 
  pivot_longer(-c(aminoAcids, C.number, Abbreviation_1), names_to = "protein", values_to = "pct") %>% 
  group_by(protein) %>% 
  mutate(pct_C.atom = C.number * pct / sum(C.number * pct) * 100) %>% 
  # put both molecule and carbon-based fraction into a single column
  pivot_longer(contains("pct"), names_to = "base", values_to = "pct") %>% 
  # mark BCAA
  mutate(whichAA = ifelse(aminoAcids %in% (BCAA <- c("Leucine", "Isoleucine", "Valine")), aminoAcids, "others")) %>% 
  # protein in order
  mutate(protein = factor(protein, levels = c("myosin I", "actin", "collagen" )))

s.all.carbon.pct


# summary, grouping non-BCAA together
s.all.carbon.pct.summary <- s.all.carbon.pct %>% 
  group_by(protein, whichAA, base) %>% 
  summarise(pct = sum(pct)) %>% 
  mutate(whichAA = factor(whichAA, levels = c("Valine", "Leucine", "Isoleucine", "others") %>% rev()))


# plot
func.plt.AA <- function(mydata){
  mydata %>% 
    ggplot(aes(protein, y = pct, fill = whichAA)) +
    facet_wrap(~base) +
    # column border for BCAA only
    geom_col(alpha = .7, color = "black") +
    # label with percentage
    geom_text(aes(label = paste(round(pct, 1), "%")),
              position = position_stack(vjust = .5),
              hjust = .5) +
    theme.myClassic +
    scale_y_continuous(expand = expansion(mult = c(0, 0)),
                       breaks = seq(0, 100, 20)) +
    scale_x_discrete(expand = expansion(add = .8)) +
    labs(x = NULL, y = "amino acid %") +
    scale_fill_manual(values = c("grey90", "orange", "steelblue1", "tomato"))
}

plt.aa.composition <- s.all.carbon.pct.summary %>% 
  filter(protein != "collagen") %>% 
  filter(base == "pct") %>% 
  func.plt.AA() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2)) +
  facet_wrap(~"")

plt.aa.composition

ggsave(filename = "amino acid distribution.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 4.5, width = 3)




# ATOM-based correction factor
favorite <- s.all.carbon.pct.summary %>%
  filter(protein != "collagen" & base == "pct_C.atom" & whichAA == "Valine")
correct.factor.aminoAcids <- 100 / mean(favorite$pct)

# MOLECULE-based corretion factor
favorite2 <- s.all.carbon.pct.summary %>%
  filter(protein != "collagen" & base == "pct" & whichAA == "Valine")
correct.factor.aminoAcids.molecule <- 100 / mean(favorite2$pct)

# plot glutamine and glutamate
s.all %>% 
  filter(aminoAcids %in% c("Glutamate", "Glutamine"))




# - <>-# - <>-# - <>-# - <>-# - <>-# - <>-# - <>-# - <>-# - <>-
# Extrapolate from C16:0, C18:1, C18:2 to total NEFA utilization

d.FA.content <- read_excel(
  "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/fatty acids content in serum.xlsm",
  skip = 1, range = "A2:AA26")


# average content
d.FA.content.summary <- d.FA.content %>% 
  select(-formula) %>% 
  pivot_longer(-c(compound, C_number), names_to = "phenotype", values_to = "conc.") %>% 
  mutate(phenotype = str_remove(phenotype, "_" %R% DIGIT) %>% factor(levels = ordered.phenotype),
         C_number = as.double(C_number)) %>% 
  group_by(compound,C_number, phenotype) %>% 
  summarise(conc.mean = mean(conc.)) 

d.FA.content.summary

# visualize TIC intensity (likely molecule basis)
d.FA.content.summary %>% ungroup() %>% 
  mutate(compound = fct_reorder(compound, conc.mean, .fun = mean)) %>% 
  ggplot(aes(x = phenotype, y = conc.mean, fill = compound)) +
  geom_col(color = "black") +
  geom_text(aes(label = compound), 
            position = position_stack(vjust = .5)) +
  scale_fill_viridis_d(option = "B") +
  theme.myClassic +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  scale_x_discrete(expand = expansion(add = .8)) +
  theme(legend.position = "none") +
  labs(y = "TIC intensity")



# calculate relative percent on an atom basis
d.FA.carbon_pct <- d.FA.content.summary %>% 
  filter(C_number <= 20) %>% 
  group_by(phenotype) %>% 
  mutate(C.pct = (conc.mean * C_number) / sum(conc.mean * C_number) * 100,
         top3 = ifelse(compound %in% c("C16:0", "C18:1", "C18:2"), compound, "others"))

# average carbon number in a FA (for later reesterification calculation)
d.FA.length.mean <- d.FA.carbon_pct %>% 
  group_by(phenotype) %>% 
  summarise(C.number_average = sum(C.pct * C_number) / sum(C.pct))
# on average, there are 17.8 atoms in a NEFA. This number will be used to calculate the reesterification rate
average.C.of.FA <- d.FA.length.mean$C.number_average %>% mean()



# group all - non-top3 NEFA
d.FA.carbon_pct.summary <- d.FA.carbon_pct %>% 
  group_by(phenotype, top3) %>% 
  summarise(C.pct = sum(C.pct)) %>% 
  mutate(top3 = factor(top3, levels = c("C16:0", "C18:1", "C18:2", "others") ))

# plot
func.plt.FA <- function(mydata, colwidth = .9){
  mydata %>% 
    ggplot(aes(x = phenotype, y = C.pct, fill = top3)) +
    geom_col(color = "black", position = position_stack(reverse = T), width =  colwidth) +
    geom_text(aes(label = round(C.pct, 0) %>% paste("%")),
              position = position_stack(vjust = .5, reverse = T),
              size = 5, color = "black") +
    labs(y = "carbon % in NEFA", x = NULL) +
    scale_fill_manual(values = c("grey95", "pink", "turquoise", "yellow") %>% rev()) +
    guides(fill = guide_legend(reverse = T)) +
    theme.myClassic + 
    scale_y_continuous(expand = expansion(mult = c(0, 0)),
                       breaks = seq(0, 100, 20)) 
}

# WT
plt.FA.composition.WT <- d.FA.carbon_pct.summary %>% 
  filter(phenotype == "WT") %>% func.plt.FA(colwidth = .5) +
  guides(fill = guide_legend(nrow = 4)) +
  coord_flip(clip = "off") +
  theme(plot.margin = margin(rep(10, 4)),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 

plt.FA.composition.WT

ggsave(filename = "fatty acid distribution WT.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 3, width = 5.5)



# all phenotypess
d.FA.carbon_pct.summary %>% func.plt.FA() +
  scale_x_discrete(labels = function(x){str_replace(x, "WT", "control")})

ggsave(filename = "fatty acid distribution 3 pheno.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 3.3, width = 4.2)

# calculate correction factor to convert sum of top3 NEFA to total NEFA
d.FA.correction.factor <- d.FA.carbon_pct.summary %>% 
  filter(top3 == "others") %>% 
  mutate(correction.factor = 100/(100 - C.pct)) %>% 
  select(phenotype, correction.factor)

d.FA.correction.factor

# create a data frame of corretion factor for amino acids and fat
d.correct.factors <- d.FA.correction.factor %>% 
  rename(fat = correction.factor) %>% 
  mutate(protein = correct.factor.aminoAcids,
         carbs = 1,
         KB = 1) %>%
  pivot_longer(-phenotype, names_to = "category", values_to = "factor")

d.correct.factors



# - <>-# - <>-# - <>-# - <>-# - <>-# - <>-# - <>-# - <>-# - <>-# - <>-

# extrapolated total OX and NOS

category <- c(
  "Glucose" = "carbs", 
  "Lactate" = "carbs",
  "Glycerol" = "carbs",
  "3-HB" = "KB",
  "Valine" = "protein",
  "Glutamine" = "protein",
  "Alanine" = "protein",
  "C16:0" = "fat",
  "C18:1" = "fat",
  "C18:2" = "fat")

ordered.category <- c("carbs", "fat", "protein", "KB")

# mark each nutrient its associated category
d.flx.CO2.nonOx.sink.categorized <- d.flx.CO2.nonOx.sink %>% 
  mutate(Compounds = as.character(Compounds), 
         category = category[Compounds]) %>% 
  # remove alanine and glutamine, as they're included in total protein
  filter(! Compounds %in% c("Alanine", "Glutamine"))

# sum the flux in each category of nutrients
d.flx.CO2.nonOx.sink.categorized.summary <- 
  d.flx.CO2.nonOx.sink.categorized %>% 
  group_by(phenotype, destiny, category) %>% 
  summarise(nmol.min.animal.mean = sum(nmol.min.animal.mean),
            nmol.min.animal.SEM = sum(nmol.min.animal.SEM^2) %>% sqrt()) %>% 
  # combine with the correction factor data frame
  left_join(d.correct.factors, by = c("phenotype", "category")) %>% 
  # calculate total flux with extrapolation
  mutate(categoryTotal.nmol.min.animal.mean = nmol.min.animal.mean * factor,
         categoryTotal.nmol.min.animal.SEM = nmol.min.animal.SEM * factor) %>% 
  # remove fluxes before applying the correction factor
  select(-c(nmol.min.animal.mean, nmol.min.animal.SEM, factor)) %>% 
  # category in order
  mutate(category = factor(category, levels = ordered.category %>% rev())) %>% 
  # mark error position
  arrange(desc(category)) %>% 
  mutate(error.y = cumsum(categoryTotal.nmol.min.animal.mean)) %>% 
  # combine with total CO2, calculate % of total CO2
  left_join(d.totalCO2) %>% 
  mutate(pct = ifelse(destiny == "CO2",
                      categoryTotal.nmol.min.animal.mean/1000/ totalCO2.umol.min * 100,
                      NA))

# plot

func.plt.category.OX.NOS <- function(mydata = d.flx.CO2.nonOx.sink.categorized.summary){
  # mydata = d.flx.CO2.nonOx.sink.categorized.summary
  mydata %>% 
    ggplot(aes(x = phenotype, y = categoryTotal.nmol.min.animal.mean, fill = category)) +
    geom_col(color = "black", alpha = .7) +
    facet_wrap(~destiny, nrow = 1) +
    geom_errorbar(aes(ymin = error.y - categoryTotal.nmol.min.animal.SEM,
                      ymax = error.y), width = .3) +
    scale_y_continuous(expand = expansion(mult = c(0, .1)),
                       breaks = seq(0, 160000, 20000),
                       labels = function(x){x/1000} ) + # to umol /min/animal)
    scale_x_discrete(
      expand = expansion(add = 1),
      labels = function(x){str_replace(x, "WT", "control")}) +
    
    theme.myClassic +
    theme(axis.title.x = element_blank(),
          panel.spacing = unit(0.1, "pt"))  +
    labs(y = "µmol C-atoms / min / animal\n") +
    guides(fill = guide_legend(reverse = F)) +
    scale_fill_manual(values = c("pink", "steelblue1", "yellow", "firebrick1"))
}

# plot WT
plt.category.OX.NOS.WT <-  func.plt.category.OX.NOS(
  mydata = d.flx.CO2.nonOx.sink.categorized.summary %>% filter(phenotype == "WT")) +
  scale_y_continuous(expand = expansion(mult = c(0, .01)),
                     breaks = seq(0, 160000, 10000),
                     labels = function(x){x/1000} ) +  # to umol /min/animal)
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  coord_cartesian(ylim = c(0, 70*1000))

plt.category.OX.NOS.WT  

# Combine with total CO2, nutrient-wise CO2/NOS, and category CO2/NOS
plot_grid(
  # total CO2
  plt.CO2.total.WT + 
    scale_x_discrete(expand = expansion(add = 1)) +
    theme(axis.text.x = element_blank(),
          panel.spacing = unit(-20, "pt"),
          axis.text = element_text(size = 15)), 
  
  # gap
  ggg <- ggplot() + theme_void(),
  
  # nutrient-wise contribution
  plt.flx.CO2.nonOx.sink.animal.WT.2.theme <- 
    plt.flx.CO2.nonOx.sink.animal.WT.2 + 
    scale_x_discrete(expand = expansion(add = 1)) +
    theme(
      axis.title.x = element_blank(), 
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text = element_text(size = 15),
      legend.text = element_text(size = 14),
      panel.spacing = unit(-20, "pt")),
  
  # gap
  ggg,
  
  # category contribution
  plt.category.OX.NOS.WT.theme <- 
    plt.category.OX.NOS.WT + 
    theme(panel.spacing = unit(-20, "pt"),
          axis.text = element_text(size = 15),
          legend.text = element_text(size = 14)), 
  
  nrow = 1, rel_widths = c(.8, .1, 2, .1, 2)
)

ggsave(filename = "OX NOS WT.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 3.7, width = 11.5)



# -

# Plot per gram body weight # -# -# -# -# -# -# -# -# -# -# -# -# -# -
# define function to remove the right side y-axis
f.removeRightYaxis <- function(p){
  p + theme(axis.text.y.right = element_blank(),
            axis.title.y.right = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y.right = element_blank())
}



# update y axis into per gram body weight using average body weight mass 29 g for WT

# <1> CO2 output
# created at 6_core_consumption_flux.R file
plt.CO2.total.WT.BW # normalized by body weight of each individual animal

plt.CO2.total.WT.BW.theme.gBW <- plt.CO2.total.WT.BW + 
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(axis.text.x = element_blank(),
        panel.spacing = unit(-20, "pt"),
        axis.text = element_text(size = 15)) +
  coord_cartesian(ylim = c(0, 2.5))

# plt.CO2.total.WT.BW.theme.gBW

# <2> 10 nutrients Ox and Sink
plt.flx.CO2.nonOx.sink.animal.WT.2.theme.gBW <- 
  f.removeRightYaxis(
    plt.flx.CO2.nonOx.sink.animal.WT.2.theme +
      scale_y_continuous(
        position = "right", name = "nmol/min/animal",
        expand = expansion(mult = c(0, .1)),
        sec.axis = sec_axis(trans = ~./29, 
                            name = "nmol C atom / min / g BW",
                            breaks = seq(0, 2400, 500)
        )
      ) +
      coord_cartesian(ylim = c(0, 2200*30))
  )

# <3> per class of nutrients Ox and Sink
plt.category.OX.NOS.WT.theme.gBW <- f.removeRightYaxis(
  plt.category.OX.NOS.WT.theme +
    scale_y_continuous(
      position = "right", name = "nmol/min/animal",
      expand = expansion(mult = c(0, .1)),
      sec.axis = sec_axis(trans = ~./29, 
                          name = "nmol C atom / min / g BW",
                          breaks = seq(0, 2500, 500)
      )
    ) + 
    coord_cartesian(ylim = c(0, 2200*30))
)

plot_grid(
  plt.CO2.total.WT.BW.theme.gBW, 
  ggg,
  plt.flx.CO2.nonOx.sink.animal.WT.2.theme.gBW, ggg,
  plt.category.OX.NOS.WT.theme.gBW,
  nrow = 1, rel_widths = c(.8, .15, 2, .15, 2)
)


ggsave(filename = "OX NOS WT gBW.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 3.7, width = 12)





# all 3 phenotype
plt.category.OX.NOS.WT.3pheno <- func.plt.category.OX.NOS() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) 
plt.category.OX.NOS.WT.3pheno

# Combine with total CO2, nutrient-wise CO2/NOS, and category CO2/NOS
plt.stacked.bars.OX.NOS.3pheno.ALL <- plot_grid(
  # 1. total CO2
  plt.CO2.total.3pheno +
    facet_wrap(~1, scales = "free_x") + 
    coord_cartesian(ylim = c(0, 140), clip = "off")+
    scale_y_continuous(breaks = seq(0, 140, 20),
                       expand = expansion(mult = c(0, .08)),
                       name = "umol C/min/animal",
                       position = "right",
                       sec.axis = sec_axis(transform = ~ . / 28,
                                           breaks = seq(0, 5, .5),
                                           labels = function(x) x * 1000,
                                           name = "nmol C / min / g lean mass\n"))+
    theme(axis.text.x = element_text(angle = 45),
          axis.text = element_text(size = 15),
          # remove the right (original main y axis elements) on unit of umol/min/animal
          axis.title.y.right = element_blank(),
          axis.text.y.right = element_blank(),
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank()), 
  #gap
  ggplot() + theme_void(),
  
  # 2. nutrient-wise contribution
  plt.flx.CO2.nonOx.sink.3pheno + 
    theme(# strip.text = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank(),
      axis.text = element_text(size = 15),
      # remove the rght (original. main y axis elements) on unit of umol/min/animal
      axis.title.y.right = element_blank(),
      axis.text.y.right = element_blank(),
      axis.line.y.right = element_blank(),
      axis.ticks.y.right = element_blank()) +
    scale_y_continuous(limits = c(0, 100000),
                       breaks = seq(0, 140000, 20000),
                       labels = function(x){x/1000},
                       expand = expansion(mult = c(0, .08)),
                       name = "umol C/min/animal",
                       position = "right",
                       sec.axis = sec_axis(transform = ~ . / 28,
                                           breaks = seq(0, 5000, 500),
                                           name = "nmol C / min / g lean mass\n")) +
    coord_cartesian(ylim = c(0, 140*1000)),
  
  #gap
  ggplot() + theme_void(),
  
  # 3. category-wise contribution
  plt.category.OX.NOS.WT.3pheno +
    # theme(strip.text = element_blank()) +
    scale_y_continuous(expand = expansion(mult = c(0, .08)),
                       labels = function(x){x/1000},
                       breaks = seq(0, 140*1000, 20*1000),
                       position = "right",
                       sec.axis = sec_axis(transform = ~ . / 28,
                                           breaks = seq(0, 5000, 500),
                                           name = "nmol C / min / g lean mass\n")) +
    coord_cartesian(ylim = c(0, 140*1000)) +
    theme(axis.text = element_text(size = 15),
          # remove the right (original main y axis elements) on unit of umol/min/animal
          axis.title.y.right = element_blank(),
          axis.text.y.right = element_blank(),
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank()),
  
  nrow = 1, rel_widths = c(1.13, .1, 2.7, .06, 2.7)
  
)

plt.stacked.bars.OX.NOS.3pheno.ALL

ggsave(filename = "OX NOS 3 pheno per g lean mass.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 5, width = 13.5)




# mark with numbers

# calcualte total NOS flux and pct of each category
d.cats.stats <- d.flx.CO2.nonOx.sink.categorized.summary %>% 
  filter(destiny == "non-Ox sink") %>% 
  group_by(destiny, phenotype, category) %>% 
  summarise(categoryTotal.nmol.min.animal.mean = sum(categoryTotal.nmol.min.animal.mean)) %>% 
  group_by(phenotype) %>% 
  mutate(pct = categoryTotal.nmol.min.animal.mean/ sum(categoryTotal.nmol.min.animal.mean))


plt.category.OX.NOS.WT.3pheno +
  # label with absolute fluxes 
  geom_text(aes(label = round(categoryTotal.nmol.min.animal.mean/1000, 1)),
            position = position_stack(vjust = .5)) +
  # label with % of the total CO2
  geom_text(
    data = d.flx.CO2.nonOx.sink.categorized.summary %>% filter(destiny == "CO2"),
    aes(label = ifelse(!is.na(pct), round(pct, 1) %>% paste("%"), NA), # NOS has NA values, not label them
        x = rep((1:3)+.4, 4)),
    position = position_stack(vjust = .5),
    color = "blue") +
  # label with total % explained for OX
  geom_text(data = d.flx.CO2.nonOx.sink.categorized.summary %>% group_by(phenotype, destiny) %>% 
              summarise(pct.total = sum(pct)) %>% filter(destiny == "CO2"),
            aes(x = 1:3, 
                y = c(60, 67, 69) * 1000, 
                label = round(pct.total, 1) %>% paste("%\nexplained")),
            inherit.aes = F, fontface = "bold", 
            color = "cyan3") +
  # label % for NOS
  geom_text(data = d.cats.stats, 
            aes(label = (pct*100) %>% round() %>% paste("%"),
                x = (c(1:3)+.4) %>% rep(each = 4)),
            position = position_stack(vjust = .5), color = "red") +
  # label TOTAL height of the bar
  geom_text(# calculate sum of all categories 
    data = d.flx.CO2.nonOx.sink.categorized.summary %>% group_by(phenotype, destiny) %>% 
      summarise(grandAll = sum(categoryTotal.nmol.min.animal.mean)), 
    aes(label = paste("sum:", (grandAll/1000) %>% round()) ,
        x = phenotype,
        y = grandAll + 5*1000),
    inherit.aes = F, color = "green3", fontface = "bold",size = 5) 


ggsave(filename = "OX NOS 3 pheno numbers.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 8, width = 12)





# -<>-# -<>-# -<>-# -<>-# -<>-# -<>-# -<>-# -<>-
# NEFA reesterification, continue using carbon atom based flux

d.ER <- d.flx.CO2.nonOx.sink.categorized.summary %>% 
  filter(category == "fat") %>% 
  # nmol carbons / min / animal bound to CO2 and extracellular reesterification 
  mutate(ER = categoryTotal.nmol.min.animal.mean, 
         ER.SEM = categoryTotal.nmol.min.animal.SEM) %>%    
  ungroup() %>% select(phenotype, ER, ER.SEM, destiny)

x1 <- d.ER %>% select(-ER.SEM) %>% spread(key = destiny, value = ER) %>% rename(ER = `non-Ox sink`)
x2 <- d.ER %>% select(-ER) %>% spread(key = destiny, value = ER.SEM) %>% rename(ER.SEM = `non-Ox sink`, CO2.SEM = CO2)

d.ER  <-  left_join(x1, x2, by = "phenotype")

# total lipolysis rate, 3 x (Ra of glycerol) (here Ra is the rate of appearance, or total production rate P)
d.lipolysis <- d.directFlux_interConversion %>% 
  filter(targetCompound == "Glycerol" & sources == "storage") %>% ungroup() %>% 
  
  # unit for glycerol, the flux is in molecules (as exception)
  mutate(Ra.glycerol.molec = nmol.min.animal.mean/3, 
         Ra.glycerol.molec.SEM = nmol.min.animal.SEM / 3 ) %>% 
  
  select(phenotype, Ra.glycerol.molec, Ra.glycerol.molec.SEM) %>% 
  
  # total NEFA release rate (carbon basis), 3 x Ra of glycerol * carbon number per fatty acid
  mutate(total.Lys = Ra.glycerol.molec * 3 * average.C.of.FA,
         total.Lys.SEM = Ra.glycerol.molec.SEM * 3 * average.C.of.FA)

# calculate intracellular reesterification as total lipolysis - circulatory production
d.reesterified <- d.ER %>% left_join(d.lipolysis, by = "phenotype") %>% 
  mutate(IR = total.Lys - CO2 - ER,
         IR.SEM = sum(CO2.SEM^2 + ER.SEM^2) %>% sqrt())

# tidy up 
a1 <- d.reesterified %>% select(phenotype, CO2, ER, IR) %>% 
  pivot_longer(-phenotype, values_to = "nmol.min.animal.mean", names_to = "destiny") 

a2 <- d.reesterified %>% select(phenotype, CO2.SEM, ER.SEM, IR.SEM) %>% 
  pivot_longer(-phenotype, values_to = "nmol.min.animal.SEM", names_to = "destiny") %>% 
  mutate(destiny = str_remove(destiny, ".SEM"))

d.reesterified.tidy <- left_join(a1, a2, by = c("phenotype", "destiny")) %>% 
  mutate(destiny = factor(destiny, levels = c("ER", "CO2", "IR"))) %>% 
  group_by(phenotype) %>% 
  arrange(desc(destiny)) %>% 
  mutate(error.y = cumsum(nmol.min.animal.mean)) %>% 
  # calculate percentage relative to total lipolysis rate
  mutate(fraction = nmol.min.animal.mean / (sum(nmol.min.animal.mean)) )



# visualize
fun.plt.esterification <-  
  func.plt.ester <- function(dataset = d.reesterified.tidy){
    dataset %>% 
      ggplot(aes(x = phenotype, y = nmol.min.animal.mean, fill = destiny)) +
      geom_col(color = "black", alpha = .7) +
      geom_errorbar(aes(ymax = error.y,
                        ymin = error.y - nmol.min.animal.SEM),
                    width = .2) +
      scale_fill_brewer(palette = "Pastel1", labels = c("IR" = "Intra", "ER" =  "Extra", "CO2" = "OX")) +
      scale_x_discrete(
        expand = expansion(add = .8),
        labels = function(x) {str_replace(x, "WT", "control")}) +
      scale_y_continuous(breaks = seq(0, 160*1000, 25*1000),
                         labels = function(x){x/1000},
                         expand = expansion(mult = c(0, 0.05))) +
      theme.myClassic +
      labs(x = NULL, y = "µmol C / min / animal")
  }


fun.plt.esterification(dataset = d.reesterified.tidy %>% 
                         filter(phenotype == "WT")) +
  scale_fill_brewer(palette = "Spectral") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

ggsave(filename = "reesterification WT.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 3.5, width = 3)

# per animal
plt.ester.3pheno <-  fun.plt.esterification(dataset = d.reesterified.tidy) 

ggsave(filename = "reesterification 3 pheno.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 3.5, width = 4)


# label the bars with percentage
plt.ester.3pheno +
  # label with absolute flux
  geom_text(aes(label = round(nmol.min.animal.mean/1000, 1)),  # nmol to µmol
            position = position_stack(vjust = .5)) +
  # label with %
  geom_text(aes(label = round(fraction*100, 1) %>% paste("%"),
                x = rep(c(1:3)+.3, time = 3)), 
            position = position_stack(vjust = .5),
            color = "red", size = 2.5) 

ggsave(filename = "reesterification numbered.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 3.5*1.5, width = 4*1.5)




# -<>-# -<>-# -<>-# -<>-# -<>-# -<>-# -<>-# -<>-# -<>-# -<>-
# Amino acid Fcirc and Sink flux 
d.aa.Fcirc <- read_excel(path = path.aa, sheet = "AA.Fcirc") %>% 
  left_join(d.aa, by = "aminoAcids") %>% 
  # convert to carbon atom basis, whole animal level
  mutate(P.nmol.min.anim = `Fcirc_nmol/g/min` * C.number * 27) %>% 
  select(aminoAcids, contains("Abbre"), status, P.nmol.min.anim)

# amino acid average percent content 
d.aa.pct <- s.all.carbon.pct %>% 
  filter(protein != "Collagen") %>% 
  filter(base == "pct_C.atom") %>% 
  group_by(aminoAcids) %>% 
  summarise(pct = mean(pct))

# combine fatty acid content with total production flux
d.aa.Fcirc <- d.aa.Fcirc %>% 
  left_join(d.aa.pct)


# sink flux
# d.release <- d.directFlux_interConversion %>% 
#   filter(phenotype == "WT" & sources == "storage" & targetCompound %in% c("Glutamine", "Alanine", "Valine")) %>% 
#   ungroup() %>%  select(targetCompound, nmol.min.animal.mean) %>% 
#   rename(store.release.nmol.min = nmol.min.animal.mean, aminoAcids = targetCompound)
#   
# d.aa.Fcirc %>% left_join(d.release)

# direct recycling flux
d.sink.measured <- d.flx.CO2.nonOx.sink %>% 
  filter(phenotype == "WT" & destiny == "non-Ox sink" & 
           Compounds %in% c("Alanine", "Glutamine", "Valine")) %>% ungroup() %>% 
  select(Compounds, nmol.min.animal.mean) %>% 
  rename(aminoAcids = Compounds, sink = nmol.min.animal.mean)

# assume same fox of EAA as valine
d.sink <- d.aa.Fcirc %>% select(aminoAcids, P.nmol.min.anim) %>% 
  filter(! aminoAcids %in% c("Alanine", "Glutamine", "Valine")) %>% 
  mutate(sink = P.nmol.min.anim * (
    1- filter(d.CO2.fox, phenotype == "WT" & infused.tracer == "Valine")$fox.mean)) %>% 
  select(-P.nmol.min.anim) %>% 
  # combine with the measured sink flux of amino acids
  rbind(d.sink.measured)

# combine the sink flux into the bigger dataset
d.aa.Fcirc.sink <- d.sink %>% left_join(d.aa.Fcirc)


# plot content vs Fcirc

# plot sink flux vs content
func.plt.aa <-  function(p){
  p + geom_smooth(data = d.aa.Fcirc.sink %>% filter(! aminoAcids %in% c ("Glutamine", "Alanine")),
                  method = "lm", se = F, color = "black",
                  aes(group = 1)) +
    geom_point(size = 10, shape = 21, alpha = .7) + 
    geom_text(aes(label = Abbreviation_1), color = "black") +
    scale_y_continuous(labels = function(x){x/1000},
                       name = "carbon µmol / min \n") +
    theme.myClassic +
    theme(legend.position = "none") +
    coord_cartesian(clip = "off") +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Pastel1") 
}


# plot Fcirc  
plt.fcirc <- d.aa.Fcirc.sink %>%  
  ggplot(aes(x = pct, y = P.nmol.min.anim, fill = status, color = status)) 
plt.fcirc <- plt.fcirc %>% func.plt.aa()

# sink plot
plt.sink <- d.aa.Fcirc.sink %>% 
  ggplot(aes(x = pct, y = sink, fill = status, color = status)) 
plt.sink <- plt.sink %>%  func.plt.aa()

plot_grid(plt.fcirc + ggtitle("total production"), 
          ggplot() + theme_void(),
          plt.sink + ggtitle("recycling flux"),
          rel_widths = c(2, .5, 2.3), nrow = 1)

ggsave(filename = "amino acid flux vs. content.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 3.8, width = 7.5)



# plot direct storage-releasing flux
d.ala.storageRelease <- d.directFlux_interConversion %>% 
  filter(phenotype == "WT" & sources == "storage" & targetCompound == "Alanine") %>% 
  select(nmol.min.animal.mean)

d.storageRelease <- d.aa.Fcirc.sink %>% filter(status == "Essential" | Abbreviation_1 == "Ala") %>% 
  rename(storage_release.min.animal = P.nmol.min.anim) 

d.storageRelease[d.storageRelease$Abbreviation_1 == "Ala", ]$storage_release.min.animal <- d.ala.storageRelease$nmol.min.animal.mean 

d.storageRelease %>% 
  ggplot(aes(x = pct, y = storage_release.min.animal, fill = status, color = status)) +
  geom_smooth(method = "lm", alpha = .3, fill = "grey", color = "black") +
  geom_point(size = 10) +
  geom_text(aes(label = Abbreviation_1), color = "black", size = 4) +
  scale_y_continuous(labels = function(x){x/28} %>% round(),
                     name = "nmol C / min / g BW \n",
                     limits = c(0, NA)) +
  scale_x_continuous(limits = c(0, NA), name = "C percentage") +
  theme.myClassic +
  theme(legend.position = "none") +
  coord_cartesian(clip = "off") +
  scale_color_brewer(palette = "Pastel1")




# Significance analysis ------------------------
# the mean
d.production.consumption2 <- 
  d.production.consumption %>% 
  ungroup() %>% 
  mutate(to = as.character(to), from = as.character(from)) %>% 
  mutate(
    # lipolysis
    from = ifelse(from == "storage" & to %in% c("C16:0", "C18:1", "C18:2", "Glycerol"), "TAG", from),
    # glycogenolysis
    from = ifelse(from == "storage" & to %in% c("Glucose", "Lactate"), "Glycogen", from),
    # protein breakdown
    from = ifelse(from == "storage" & to %in% c("Glutamine", "Alanine", "Valine"), "Protein", from),
    
    # FA storing
    to = ifelse(from %in% c("C16:0", "C18:1", "C18:2", "Glycerol") & to == "non-Ox sink", "TAG", to),
    # glucose storing
    to = ifelse(from == "Glucose" & to == "non-Ox sink", "Storage", to),
  )

d1 <-  d.production.consumption2 %>% 
  select(phenotype, to, from, nmol.min.animal.mean) %>% 
  mutate(nmol.min.animal.mean = nmol.min.animal.mean + 1) # to avoid zero values

d1 <- d1[!duplicated(d1), ] %>% 
  spread(phenotype, nmol.min.animal.mean)

# the SEM
d2 <-  d.production.consumption2 %>% 
  select(phenotype, to, from, nmol.min.animal.SEM)

d2 <- d2[!duplicated(d2), ] %>% 
  spread(phenotype, nmol.min.animal.SEM) 
colnames(d2)[3:5] <- colnames(d2)[3:5] %>% str_c("_se")
d2

# combine
d.sig <- d1 %>% left_join(d2, by = c("to", "from")) %>% rename(ob = "ob/ob", ob_se = `ob/ob_se`)
d.sig <- d.sig %>% 
  mutate(ratio_HFD_WT = HFD / WT,
         ratio_ob_WT = ob / WT,
         t_HFD_WT = (HFD - WT) / sqrt( WT_se^2 + HFD_se^2 ),
         t_ob_WT =  (ob - WT) / sqrt( WT_se^2 + ob_se^2 )) %>% 
  
  # calculate p value (two tails) based on t value (Inf. approach normal distribution)
  mutate(across(c(t_HFD_WT, t_ob_WT), 
                # calculated accumulated p value, but taking the upper tail, i.e., get 1-accumulated probability
                .fns = ~ pt(abs(.), df = Inf, lower.tail = FALSE, log.p = FALSE) * 2,
                .names = "{.col}_pvalue")) 
d.sig

# tidy up
# flux (mean of all 3 phenotypes)
d.flx <- d.sig %>% select(to, from, WT, HFD, ob) %>% 
  pivot_longer(-c(1:2), names_to = "phenotype", values_to = "flux") %>% 
  group_by(to, from) %>% 
  summarise(flux = mean(flux))
# fold change
d.ratio <- d.sig %>% select(to, from, contains("ratio")) %>% 
  pivot_longer(-c(1:2), names_to = "compare", values_to = "ratio") %>% 
  mutate(compare = str_remove(compare, "ratio_"))
# p value
d.p <- d.sig %>% select(to, from, contains("pvalue")) %>% 
  pivot_longer(-c(1:2), names_to = "compare", values_to = "pvalue") %>% 
  mutate(compare = str_remove(compare, "t_") %>% str_remove("_pvalue"))

# adjusted p value
d.p <- d.p %>% mutate(
  # p.adj = p.adj(pvalue, method = "bonferroni")
  p.adj = p.adjust(pvalue, method = "hochberg")
  
) %>% 
  # to facilitate visualization, 
  # bring the ob/ob Glycogen to Glucose to plotting range: -log10 p value from 32 to 24 
  mutate(p.adj = ifelse(
    compare == "ob_WT" & from == "Glycogen" & to == "Glucose", 10^-24, p.adj))

# combine average flux, fold change, and adjusted p value
d.volcano <- d.p %>% 
  left_join(d.ratio, by = c("to", "from", "compare")) %>% 
  left_join(d.flx, by = c("to", "from")) %>% 
  mutate(neg.log.p.adj = - log10(p.adj),
         neg.log.pvalue = -log10(pvalue)) 

# visualize
p.threshold <- -log10(0.05)


func.plt.volcano <- function(d, whichcompare = "ob_WT", textsize = 3) {
  
  # whichcompare = "ob_WT"
  d <- d.volcano %>% filter(compare == whichcompare)
  
  # significant subset
  d.significant.subset <- d %>% 
    # filter(neg.log.p.adj > 0 ) # %>% # stats significant
    filter(neg.log.p.adj > p.threshold ) %>%  # %>% # stats significant
    filter(flux > 200) # flux > 200
  
  d %>% 
    ggplot(aes(x = log2(ratio), 
               # y = neg.log.pvalue, 
               y = neg.log.p.adj,
               size = flux, color = from)) +
    
    geom_point(alpha = .4) +
    geom_hline(yintercept = p.threshold, linetype = "dashed") +
    geom_vline(xintercept = 0) +
    
    # label significant points
    ggrepel::geom_text_repel(
      data = d.significant.subset,
      aes(label = paste(from, "→", to)),
      size = textsize, max.overlaps = Inf, show.legend = F,
      fontface = "bold") +
    
    scale_size(range = c(1, 10), 
               limits = c(0, NA), 
               breaks = c(100, 500, 2000, 5000, 10000, 20000 ),
               labels = function(x){x/1000},
               name = "µmol/min") +
    scale_x_continuous(breaks = seq(-6, 6, 1)) +
    scale_y_continuous(breaks = seq(0, 30, 5)) +
    guides(
      #color = guide_legend(override.aes = list(size = 4, alpha = .7), reverse = T)
      color = "none"
      
    ) +
    theme_classic(base_size = 17) +
    theme(axis.text = element_text(colour = "black"),
          axis.title = element_text(colour = "black", face = "bold")) +
    labs(y = "-log10 (adjusted p value)", 
         x = paste("log2 (", whichcompare, ")" ))
}

p1 <- func.plt.volcano(whichcompare = "HFD_WT", textsize = 3.5) +
  coord_cartesian(xlim = c(-3., 3)) +
  labs(x = "log2 (HFD vs. WT)")

p2 <- func.plt.volcano(whichcompare = "ob_WT", textsize = 3) +
  coord_cartesian(xlim = c(-3, 5), ylim = c(0, 25)) +
  labs(x = "log2 (ob/ob vs. WT )")

plot_grid(p1 + theme(legend.position = "none"), 
          p2, 
          nrow = 1, rel_widths = c(3, 3.5))

ggsave(filename = "significant analysis.png",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 6, width = 12, device = "png")



## - <>- needed C3 flux for glycerol backbone in TAG resynthesis

d.reesterified.tidy %>% 
  filter(destiny != "CO2") %>% 
  group_by(phenotype) %>% 
  summarise(
    total.esterify.nmol.min.anim = sum(nmol.min.animal.mean)) %>% 
  mutate(glycerol.needed.nmol.min.animal = 
           # divide (carbon # per FA) divide (3 FA per TAG ) x 3 cabon/glycerol
           total.esterify.nmol.min.anim / average.C.of.FA / 3 * 3)

# mass of extracellular reesterification, carbon mass only
mass.TG.ER_carbon <- 42.8 * 12 / (10^6) * 60 * 24 # C mass in 24h

# mass (g) of equivalent TAG formed by extracellular reesterification
mass.TG.ER <- mass.TG.ER_carbon / .86


# Triolein
# C57H104O6
MW <- 12*57 + 10 + 16 *6
(12*57) / MW
# carbons account for 86% mass in a TAG  


# mass (g) of equivalent TAG formed by both Extra- and intracellular reesterification
mass.TAG_IR.ER <- (42.8+26.1) * 12 / (10^6) * 60 * 24 / .86



# mmol of TAG total formed by both ER and IR
mmol.TAG.ER.IR <- mass.TAG_IR.ER / MW * 1000
mmol.TAG.ER.IR

# mmol of equivalent glucose needed to replenish glycerol backbone
mmol.TAG.ER.IR/2

# glucose carbons sink: 2.6 umol carbons / min /animal
2.6 * 60 * 24 / 1000 / 6

# glycerol sink : 2 umol carbons / min /animal
2 * 60 * 24 / 1000 / 3

# glucose production: mmol molecules / day / animal 
glc.produced = 11507 / 6 / (10^6) * 60 * 24
glc.produced

mmol.TAG.ER.IR/2 / glc.produced

18/24



# calculate futile cycles
# dataset containing all flower fluxes
d.from.to <- d.production.consumption %>% ungroup() %>% 
  select(phenotype, from, to, nmol.min.animal.mean, nmol.min.animal.SEM) %>% 
  distinct() %>% 
  arrange(phenotype)

d.from.to

# energy wasted in Cori cycle (lactate to glucose)
d.Cori <- d.from.to %>% 
  filter(from == "Lactate" & to == "Glucose") %>% 
  mutate(cons_ATP.mean = nmol.min.animal.mean / 6 * 6,  # each glucose or two lactates molecules consumes 6 ATP
         cons_ATP.SEM = nmol.min.animal.SEM / 6 * .6) %>% 
  select(phenotype, contains("ATP")) %>% 
  mutate(futile = "Cori")
d.Cori

# glycolysis of circulating glucose
d.glycolysis.circ.glu <- d.from.to %>% 
  filter( (from == "Glucose" & to == "Lactate") |   # glucose → lactate 
            (from == "Glucose" & to == "CO2") )   # glucose → CO2
d.glycolysis.circ.glu

# glycolysis of glycogen storage → lactate
d.glycolysis.GLN.lactate <- d.from.to %>% 
  filter( (from == "storage" & to == "Lactate")) 
d.glycolysis.GLN.lactate


# total glycolysis flux
d.glycolysis <- rbind(d.glycolysis.circ.glu, d.glycolysis.GLN.lactate)
d.glycolysis

# sum up total glycolysis
d.glycolysis.total <- d.glycolysis %>% 
  group_by(phenotype) %>% 
  summarise(nmol.min.animal.mean = sum(nmol.min.animal.mean),
            nmol.min.animal.SEM = sum(nmol.min.animal.SEM^2) %>% sqrt())

d.ATP.glycolysis <-  d.glycolysis.total %>% 
  # per glucose generates 2 NADH and 2 ATP, total 7 ATP equivalent
  mutate(prod_ATP.mean = nmol.min.animal.mean / 6 * 7, 
         prod_ATP.SEM = nmol.min.animal.SEM / 6 * 7) %>% 
  select(-contains("nmol"))

# TCA 

d.ATP.TCA.gluc.lactate <- d.from.to %>% 
  filter( (from == "Glucose" & to == "CO2") |  
            (from == "Lactate" & to == "CO2") ) %>% 
  group_by(phenotype) %>% 
  summarise(nmol.min.animal.mean = sum(nmol.min.animal.mean),
            nmol.min.animal.SEM = sum(nmol.min.animal.SEM^2) %>% sqrt()) %>% 
  # per pyruvate generates 12.5 ATP flowing into the TCA cycle (including the formation of acetyl-CoA)
  mutate(prod_ATP.mean = nmol.min.animal.mean / 3 * 12.5, 
         prod_ATP.SEM = nmol.min.animal.SEM / 3 * 12.5) %>% 
  select(-contains("nmol"))

d.ATP.TCA.gluc.lactate

# sum of TCA and glycolysis of both glucose and lactate
d.ATP.glycolysis.TCA <- (d.ATP.glycolysis %>% mutate(path = "Glycolysis")) %>% 
  rbind(
    d.ATP.TCA.gluc.lactate %>% mutate(path = "TCA")) 

d.ATP.Glc.Lac.total <- d.ATP.glycolysis.TCA %>% 
  group_by(phenotype) %>% 
  summarise(prod_ATP.mean = sum(prod_ATP.mean),
            prod_ATP.SEM = sum(prod_ATP.SEM^2) %>% sqrt())

d.ATP.Glc.Lac.total

# compare the cori cycle with total energy produced
d.Cori %>% left_join(d.ATP.Glc.Lac.total, by = "phenotype") %>% 
  mutate(waste.frc = cons_ATP.mean / prod_ATP.mean)



# save project
save.image(file = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/8_core_flux extrapolation.RData")



# compare needed glycerol backbone for TAG synthesis
d.C3.needed <- d.reesterified.tidy %>% filter(destiny != "CO2") %>% 
  group_by(phenotype) %>% 
  select(-c(error.y, fraction)) %>% 
  summarise(total.reesterification.C.atom = sum(nmol.min.animal.mean),
            total.reesterification.C.atom.SEM = sum(nmol.min.animal.SEM^2) %>% sqrt()) %>% 
  mutate(C3     = total.reesterification.C.atom     / (average.C.of.FA * 3),
         C3.SEM = total.reesterification.C.atom.SEM / (average.C.of.FA * 3), 
         Compounds = factor("Reester"),
         direction = "needed",
         .keep = "unused") %>% 
  relocate(Compounds, .after = phenotype) %>% 
  # add cumulated error position
  mutate(error.y = C3)


# depositing sources of C3 from glucose, glycerol, and lactate
d.C3.source <- d.destiny %>% ungroup() %>% 
  filter(destiny == "non-Ox sink", Compounds %in% c("Glucose", "Lactate", "Glycerol")) %>% 
  select(phenotype, Compounds, nmol.min.animal.mean, nmol.min.animal.SEM) %>% 
  mutate(C3 = nmol.min.animal.mean / 3,
         C3.SEM = nmol.min.animal.SEM / 3) %>% 
  select(-c(contains("nmol"))) %>% 
  mutate(direction = "source") %>% 
  # add cumulated error position
  group_by(phenotype) %>% 
  arrange(phenotype, desc(Compounds)) %>% 
  mutate(error.y = cumsum(C3))


# Combine dataset together
d.C3.all <- d.C3.source %>% rbind(d.C3.needed) %>% 
  mutate(phenotype = str_replace(phenotype, "WT", "control"))

# plot
d.C3.all %>% 
  ggplot(aes(x = direction, y = C3, fill = Compounds)) +
  geom_col(alpha = .6, color = "black") +
  geom_errorbar(aes(ymax = error.y, 
                    ymin = error.y - C3.SEM),
                width = .3) + 
  facet_wrap(~phenotype) +
  scale_fill_viridis_d(option = "H") +
  scale_x_discrete(expand = expansion(add = .8)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)),
                     breaks = seq(0, 3000, 500),
                     labels = function(x){x/1000}) +
  labs(y = "C3 flux ( µmol / min )", x = NULL) +
  theme.myClassic +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y.left = element_text(margin = margin(r = 15)),
        strip.text = element_text(size = 14, face = "bold")) 

ggsave(filename = "glycerol recycle.pdf",
       path = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/R Figures",
       height = 3*1.1, width = 6*1.1)


save.image(file = "/Users/boyuan/Desktop/Harvard/Manuscript/1. fluxomics/raw data/8_core_flux_extrapolation.RData")


