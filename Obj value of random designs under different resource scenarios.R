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
repi.obj = detectCores() - 1



random.design <- matrix(runif(300*21), 300, 21)

ofv.result <- matrix(NA, 300, 8)


for(repi in 1:nrow(random.design))
{
  print(paste("repi = ", repi, sep = ""))
  
  ofv.result[repi, 1] <- obj.fun.par(random.design[repi,], dm.data = dm.realization, test.total = dat.agg$Z*0.5,  goal = "all", repi = repi.obj)
  
  ofv.result[repi, 2] <- obj.fun.par(random.design[repi,], dm.data = dm.realization, test.total = dat.agg$Z, goal = "all", repi = repi.obj)
  
  ofv.result[repi, 3] <- obj.fun.par(random.design[repi,], dm.data = dm.realization, test.total = dat.agg$Z*2, goal = "all", repi = repi.obj)
  
  ofv.result[repi, 4] <- obj.fun.par(random.design[repi,], dm.data = dm.realization, test.total = dat.agg$Z*5, goal = "all", repi = repi.obj)
  
  ofv.result[repi, 5] <- obj.fun.par(random.design[repi,], dm.data = dm.realization, test.total = dat.agg$Z*0.5, goal = "severe", repi = repi.obj)
  
  ofv.result[repi, 6] <- obj.fun.par(random.design[repi,], dm.data = dm.realization, test.total = dat.agg$Z, goal = "severe", repi = repi.obj)
  
  ofv.result[repi, 7] <- obj.fun.par(random.design[repi,], dm.data = dm.realization, test.total = dat.agg$Z*2, goal = "severe", repi = repi.obj)
  
  ofv.result[repi, 8] <- obj.fun.par(random.design[repi,], dm.data = dm.realization, test.total = dat.agg$Z*5, goal = "severe", repi = repi.obj)
  
  save(ofv.result , file = "Result\\OBJ_testlimit_cor.Rdata")
}

