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
GA.sugguestion <- matrix(NA, 7, 21)

dat.agg.alter <- dat.6 %>%
  mutate(Z_mild = Z_mild_1 + Z_mild_2 + Z_mild_3, Z_severe = Z_severe_1 + Z_severe_2 + Z_severe_3, inc = Y/pop*100000, inc_severe = Y_severe/pop*100000) %>%
  group_by(Year, Pref) %>%
  summarize(Y = sum(Y), Y_severe = sum(Y_severe), Z_mild = sum(Z_mild), Z_severe = sum(Z_severe), pop = mean(pop)) %>%
  mutate(inc = Y/pop*100000, inc_severe = Y_severe/pop*100000) %>%
  group_by(Pref) %>%
  summarize(Y = sum(Y), Y_severe = sum(Y_severe), Z_mild = sum(Z_mild), Z_severe = sum(Z_severe), pop = mean(pop), inc = mean(inc), inc_severe = mean(inc_severe)) %>%
  mutate(Z = Z_mild + Z_severe)

# Append the data to the shapefile for plotting maps
GA.sugguestion[1, 1:20] <- convert.inverse(21, dat.agg.alter$Z/sum(dat.agg.alter$Z))   # Existing
GA.sugguestion[2, 1:20] <- convert.inverse(21, rep(1/21, 21))   # Equal
GA.sugguestion[3, 1:20] <- convert.inverse(21, dat.agg.alter$pop/sum(dat.agg.alter$pop)) # PopSize
GA.sugguestion[4, 1:20] <- convert.inverse(21, dat.agg.alter$Y/sum(dat.agg.alter$Y))  # Case
GA.sugguestion[5, 1:20] <- convert.inverse(21, dat.agg.alter$inc/sum(dat.agg.alter$inc))  # IncRate
GA.sugguestion[6, 1:20] <- convert.inverse(21, dat.agg.alter$Y_severe/sum(dat.agg.alter$Y_severe))  # CaseSevere
GA.sugguestion[7, 1:20] <- convert.inverse(21, dat.agg.alter$inc_severe/sum(dat.agg.alter$inc_severe))  # IncRate Severe

GA.sugguestion[,21] <- 0.5


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

ga.monitornew <- function (object, digits = getOption("digits"), ...) 
{
  fitness <- na.exclude(object@fitness)
  sumryStat <- c(mean(fitness), max(fitness))
  sumryStat <- format(sumryStat, digits = digits)
  best.des <- round( unlist(object@population[which.max(object@fitness), ]), 5)
#  save(object, file = "GA_intermediate_all_TestX5.Rdata")
  cat(paste("GA | iter =", object@iter, "| Mean =", sumryStat[1], 
            "| Best =", sumryStat[2]), "| Best solution = ", best.des, "| Time = ", as.character(Sys.time()))
  cat("\n")
  flush.console()
}

GA.All.X5 <- ga(type = "real-valued", fitness = obj.fun.par, dm.data = dm.realization, test.total = dat.agg$Z*5, goal = "all", repi = repi, lower = rep(0, 21), upper = rep(1, 21),  keepBest = TRUE, popSize = 50, maxiter = 100, monitor = ga.monitornew, run = 100, pmutation = 0.05, suggestions = GA.sugguestion)

save(GA.All.X5, file = "Result\\Optimization\\GA_all_TestX5.Rdata")

stopCluster(cl)   


