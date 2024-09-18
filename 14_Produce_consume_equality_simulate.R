library(tidyverse)

n.nutrient <- 20 # compounds noted as a big letters
nutrientNames <- LETTERS[1:n.nutrient]

# simulate the n x n labeling matrix

# matrix of labeling, an upper bound of 1 (literally it can be any finite number); 
# in practice, the infusion rate is controlled such that labeling is bound below 20% to minimize metabolic perturbation
m.L <- runif(n.nutrient^2, min = 0, max = 1) %>%  
  matrix(nrow = n.nutrient)
m.L

# optional validation #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# load("/Users/boyuan/Desktop/Harvard/Research/db db mice/Infusion data/Flux analysis script/labeling matrix.RData")
# m.L <- m.L[, -1] %>% as.matrix()
#$$$$$$$$$$$$$$$#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


rownames(m.L) <- nutrientNames
colnames(m.L) <- nutrientNames



# in each row, let the diagonal value (same tracer-tracee pair) be the max;
# in each row, if the diagonal value is not the max, swap the diagonal value with the max value position
for (i in 1:n.nutrient){ 
  if (m.L[i, i] != max(m.L[i, ])) { 
    position <- m.L[i, ] %>% which.max() # column index of the max value
    biggest <- m.L[i, position] # extract max value 
    
    # switch the max and diagonal values
    m.L[i, position] <- m.L[i, i] 
    m.L[i, i] <- biggest
  }
}


# simulate the infusion rate, from 1 to 1000
R <- runif(n.nutrient, min = 1, max = 1000)

# optional validation #$$$$$$$$$$$$$$$
# load("/Users/boyuan/Desktop/Harvard/Research/db db mice/Infusion data/Flux analysis script/infusion_nmolC.min.gBW.RData")
# R <- R[-7]
#$$$$$$$$$$$$$$$#$$$$$$$$$$$$$$$$$$$$$




# simulate the fox values
fox <- runif(n.nutrient, min = 0, max = 1)
# optinal validation. #$$$$$$$$$$$$$$$
# load("/Users/boyuan/Desktop/Harvard/Research/db db mice/Infusion data/Flux analysis script/fox.RData")
# fox <- fox[-7]
#$$$$$$$$$$$$$$$#$$$$$$$$$$$$$$$$$$$$$





# calculate the production flux 
# for each target nutrient, e.g., nutrient D,
# move the Dth row and Dth column to 1st row and column, respectively. 
# the following rows and column correspond to source nutrients

t.produce <- tibble()

for (i in 1:n.nutrient){
  
  m.L.i <- m.L # m.L.i fresh start to be manipulated in each loop
  
  # i <- 3
  TN <- nutrientNames[i] # target nutrient
  
  # TN row as first row
  row.i <- m.L.i[i, ] %>% t() # extract ith row of TN
  rownames(row.i) <- TN # assign row name as TN
  m.L.i <- rbind(row.i, m.L.i[-i, ]) # TN row as first row
  
  # TN column as first column
  col.i <- m.L.i[, i] %>% as.matrix()
  colnames(col.i) <- TN
  m.L.i <- cbind(col.i, m.L.i[, -i])
  
  # create the L.TN.TN - L.i.TN - composed matrix (naming: superscript - subscript)
  # take the first column of m.L.i, and turn it into matrix
  L1 <- m.L.i[, 1] %>% rep(n.nutrient) %>% matrix(nrow = n.nutrient, byrow = F)
  
  # create the classical pairwise labeling matrix
  # excluding the first column of m.L.i, and augmented with last column of 0
  col.zeros <- rep(0, n.nutrient) %>% as.matrix()
  colnames(col.zeros) <- "SP" # SP for storage pool
  L2 <- m.L.i[, -1] %>% cbind(col.zeros)
  
  # define M as the left-side matrix result
  M <- (L1 - L2) / (1 - L2)
  
  # create vector v on the right side of equation
  R.i <- R[i] # infusion rate of TN tracer
  L.TN.TN <- L1[1, 1] # labeling of TN under infusion of TN tracer
  v <- c(R.i * (1-L.TN.TN), rep(0, n.nutrient-1))
  
  # solve the equation without constraint
  flx.produce <- solve(M, v) 
  
  # tidy up in the form of tibble
  t.i <- tibble(from =  c(rownames(m.L.i)[-1], paste0("SP.", TN) ),
         to = TN,
         flx = flx.produce)
  
  t.produce <- rbind(t.produce, t.i)
}

# total production flux
t.total.P <- t.produce %>% 
  group_by(to) %>% 
  summarise(total.produce = sum(flx)) %>% 
  rename(nutrient = to)


t.total.P$total.produce * 28/1000 # umol/min/animal
  

#- <> -- <> -- <> -- <> -- <> -- <> -- <> -- <> -- <> -- <> -- <> -- <> -- <> -

# calculate consumption flux

# create the total production matrix using the results calculated above
m.P <- t.total.P$total.produce %>% rep(each = n.nutrient) %>% 
  matrix(byrow = F, nrow = n.nutrient)
m.P
colnames(m.P) <- t.total.P$nutrient

# the all-inclusive labeling matrix is already simulated at the very early step
m.L

# calculate the left side matrix
M <- m.P * m.L / (1 - m.L)

# create the right side vector
R.fox <- R * fox
R.sink <- R * (1 - fox)

# solve the fraction of oxidation theta, and fraction of non-Ox sink lambda
θ <- solve(M, R.fox)
λ <- solve(M, R.sink)

# calculate the flux to CO2 and non-ox sink
t.total.CO2.sink <- t.total.P %>% 
  mutate(θ = θ, 
         λ = λ) %>% 
  mutate(flx.CO2 = total.produce * θ,
         flx.sink = total.produce * λ)
  
t.total.CO2.sink


# tidy up the CO2 and sink fluxes
t.destiny <- t.total.CO2.sink %>% 
  select(nutrient, contains("flx")) %>% 
  pivot_longer(-nutrient, 
               names_to = "to",
               values_to = "flx") %>% 
  rename(from = nutrient) %>% 
  mutate(to = str_remove(to, "flx."))

t.destiny


# integrate the CO2 and sink fluxes with the inter-conversion (production) flux data
t.totalConsume <- t.produce %>% filter(! str_detect(from, "SP")) %>% 
  rbind(t.destiny) %>% 
  arrange(from) %>% 
  group_by(from) %>% 
  summarise(total.consume = sum(flx)) %>% 
  rename(nutrient = from)

t.totalConsume

# Compare total consumption vs. production
t.total.P %>% 
  mutate(totalConsume = t.totalConsume$total.consume) %>% 
  mutate(diff = round(total.produce - totalConsume, 4)) %>% 
  as.data.frame()





