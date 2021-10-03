library(GA)
library(tidyverse)
library(lubridate)
library(rstan)
library(doSNOW)
library(doParallel)
library(abind)
library(Rcpp)


source("Functions.R")

load("Result\\Disease model\\Disease_model_realizations_6.RData")

dat.6 <- read_csv(".\\Data\\Sichuan_pref_yearly.csv") %>% 
  as.data.frame() %>%                                       # filter to data before 2015
  filter(Year < 2015)

# get the historical number of typing
dat.agg <- dat.6 %>%
  mutate(Z = Z_mild_1 + Z_mild_2 + Z_mild_3 + Z_severe_1 + Z_severe_2 + Z_severe_3) %>%
  group_by(Year) %>%
  summarize(Z = sum(Z))

# archetypal designs
GA.sugguestion <- matrix(NA, 7, 20)

dat.agg.alter <- dat.6 %>%
  mutate(Z_mild = Z_mild_1 + Z_mild_2 + Z_mild_3, Z_severe = Z_severe_1 + Z_severe_2 + Z_severe_3, inc = Y/pop*100000, inc_severe = Y_severe/pop*100000) %>%
  group_by(Year, Pref) %>%
  summarize(Y = sum(Y), Y_severe = sum(Y_severe), Z_mild = sum(Z_mild), Z_severe = sum(Z_severe), pop = mean(pop)) %>%
  mutate(inc = Y/pop*100000, inc_severe = Y_severe/pop*100000) %>%
  group_by(Pref) %>%
  summarize(Y = sum(Y), Y_severe = sum(Y_severe), Z_mild = sum(Z_mild), Z_severe = sum(Z_severe), pop = mean(pop), inc = mean(inc), inc_severe = mean(inc_severe)) %>%
  mutate(Z = Z_mild + Z_severe)


GA.sugguestion[1, ] <- convert.inverse(21, dat.agg.alter$Z/sum(dat.agg.alter$Z))   # Existing
GA.sugguestion[2, ] <- convert.inverse(21, rep(1/21, 21))   # Equal
GA.sugguestion[3, ] <- convert.inverse(21, dat.agg.alter$pop/sum(dat.agg.alter$pop)) # PopSize
GA.sugguestion[4, ] <- convert.inverse(21, dat.agg.alter$Y/sum(dat.agg.alter$Y))  # Case
GA.sugguestion[5, ] <- convert.inverse(21, dat.agg.alter$inc/sum(dat.agg.alter$inc))  # IncRate
GA.sugguestion[6, ] <- convert.inverse(21, dat.agg.alter$Y_severe/sum(dat.agg.alter$Y_severe))  # CaseSevere
GA.sugguestion[7, ] <- convert.inverse(21, dat.agg.alter$inc_severe/sum(dat.agg.alter$inc_severe))  # IncRate Severe


# # code for running on Savio
# ncoresPerNode <-as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))
# nodeNames <-strsplit(Sys.getenv("SLURM_NODELIST"), ",")[[1]]
# machines=rep(nodeNames, each = ncoresPerNode)
# repi = length(machines)
# repi
# 
# cl = makeCluster(machines, type = "SOCK", outfile = "log.txt")
# registerDoSNOW(cl)
# clusterExport(cl, varlist = ls(globalenv()))
# clusterEvalQ(cl, library("rstan"))
# clusterEvalQ(cl, library("abind"))

cl <- makeCluster(detectCores() - 1)
registerDoSNOW(cl)
repi = detectCores() - 1


p_severe <- seq(0.01, 0.99, 0.02)

result.all <- expand.grid(p_severe = p_severe, goal = c("all", "severe"), test.limit = "X1", design = c("Existing", "Equal", "PopSize", "Case", "IncRate", "Case_Severe", "IncRate_Severe"))
result.all$ofv <- NA



for(k in 1:nrow(result.all))
{
  print(result.all[k,])
  
  current.design <- array(NA, 21)
  current.design[21] <- result.all[k, "p_severe"]
  
  if(result.all[k,"design"] == "Existing")
  {
    current.design[1:20] <- GA.sugguestion[1, ] 
  }
  if(result.all[k,"design"] == "Equal")
  {
    current.design[1:20] <- GA.sugguestion[3, ] 
  }
  if(result.all[k,"design"] == "PopSize")
  {
    current.design[1:20] <- GA.sugguestion[5, ] 
  }
  if(result.all[k,"design"] == "Case")
  {
    current.design[1:20] <- GA.sugguestion[2, ] 
  }
  if(result.all[k,"design"] == "IncRate")
  {
    current.design[1:20] <- GA.sugguestion[4, ] 
  }
  if(result.all[k,"design"] == "Case_Severe")
  {
    current.design[1:20] <- GA.sugguestion[6, ] 
  }
  if(result.all[k,"design"] == "IncRate_Severe")
  {
    current.design[1:20] <- GA.sugguestion[7, ] 
  }
  
  result.all[k, "ofv"] <- obj.fun.par(current.design, dm.data = dm.realization, test.total = dat.agg$Z, goal = result.all[k,"goal"], repi = repi)
  
  save(result.all, file = "Result\\Optimal_psevere_benchmark_scenario.Rdata")
}

stopCluster(cl)
