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

dat <- read_csv(".\\Data\\Sichuan_pref_yearly.csv") %>% 
  as.data.frame()

# get the historical number of typing
dat.agg <- dat %>%
  mutate(Z = Z_mild_1 + Z_mild_2 + Z_mild_3 + Z_severe_1 + Z_severe_2 + Z_severe_3) %>%
  group_by(Year) %>%
  summarize(Z = sum(Z))

# archetypal designs
GA.sugguestion <- matrix(NA, 7, 20)

dat.agg.alter <- dat %>%
  filter(Year < 2015) %>%
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
# 
# cl = makeCluster(machines, type = "SOCK", outfile = "log.txt")
# registerDoSNOW(cl)
# clusterExport(cl, varlist = ls(globalenv()))
# clusterEvalQ(cl, library("rstan"))
# clusterEvalQ(cl, library("abind"))

cl <- makeCluster(detectCores() - 1)
registerDoSNOW(cl)

rep.times = 1000




#======================= for all =========================
load("Result\\Optimization\\GA_all_TestX1.Rdata")
optimal.inc <- GA.All.X1@solution 
optimal.inc[21] <- 0.11

optimal.psevere <- read.csv("Result\\Optimal_psevere_for_benchmarking.csv") %>%
  filter(goal == "(A) For all")

Alternative.name <- c("Optimal", "Existing", "Equal", "PopSize", "Case", "IncRate", "CaseSevere","IncRateSevere")

Alternative.design <- list(optimal.inc,
                           c(GA.sugguestion[1,], optimal.psevere$p_severe[optimal.psevere$design == "Existing"]), 
                           c(GA.sugguestion[2,], optimal.psevere$p_severe[optimal.psevere$design == "Equal"]),
                           c(GA.sugguestion[3,], optimal.psevere$p_severe[optimal.psevere$design == "PopSize"]),
                           c(GA.sugguestion[4,], optimal.psevere$p_severe[optimal.psevere$design == "Case"]),
                           c(GA.sugguestion[5,], optimal.psevere$p_severe[optimal.psevere$design == "IncRate"]),
                           c(GA.sugguestion[6,], optimal.psevere$p_severe[optimal.psevere$design == "SevereCase"]),
                           c(GA.sugguestion[7,], optimal.psevere$p_severe[optimal.psevere$design == "SevereIncRate"])
)



#--------- 2009-2014 -----------
load("Result\\Disease model\\Disease_model_realizations_6.RData")
result.all <- data.frame(Type = rep(Alternative.name, each = rep.times), ofv = NA)

for(k in 1:8)
{
  result.all$ofv[result.all$Type == Alternative.name[k]] <- foreach(i = 1:rep.times, .combine = "c",  .packages = c("rstan", "abind")) %dopar% {
    cat(Alternative.name[k], i, "\n")
    obj.fun(Alternative.design[[k]], dm.data = dm.realization, test.total = dat.agg$Z[1:6], goal = "all")
  }
}
save(result.all, file = "Result\\Compare_OFV_of_benchmark_all_09to14.Rdata")




#--------------- 2015 ------------------
load("Result\\Disease model\\Disease_model_realizations_7.RData")
result.all <- data.frame(Type = rep(Alternative.name, each = rep.times), ofv = NA)
for(k in 1:8)
{
  result.all$ofv[result.all$Type == Alternative.name[k]] <- foreach(i = 1:rep.times, .combine = "c",  .packages = c("rstan", "abind")) %dopar% {
    cat(Alternative.name[k], i, "\n")
    inter.result <- obj.fun(Alternative.design[[k]], dm.data = dm.realization, test.total = dat.agg$Z, goal = "all", raw = TRUE)
    mean(abs(inter.result[[1]][,,7] - inter.result[[2]][,,7]))
  }
}
save(result.all, file = "Result\\Compare_OFV_of_benchmark_all_15.Rdata")






#======================= for severe =========================
load("Result\\Optimization\\GA_severe_TestX1.Rdata")
optimal.inc.severe <- GA.severe.X1@solution 
optimal.inc.severe[21] <- 0.65

optimal.psevere <- read.csv("Result\\Optimal_psevere_for_benchmarking.csv") %>%
  filter(goal == "(B) For severe")

Alternative.name <- c("Optimal", "Existing", "Equal", "PopSize", "Case", "IncRate", "CaseSevere","IncRateSevere")

Alternative.design <- list(optimal.inc.severe,
                           c(GA.sugguestion[1,], optimal.psevere$p_severe[optimal.psevere$design == "Existing"]), 
                           c(GA.sugguestion[2,], optimal.psevere$p_severe[optimal.psevere$design == "Equal"]),
                           c(GA.sugguestion[3,], optimal.psevere$p_severe[optimal.psevere$design == "PopSize"]),
                           c(GA.sugguestion[4,], optimal.psevere$p_severe[optimal.psevere$design == "Case"]),
                           c(GA.sugguestion[5,], optimal.psevere$p_severe[optimal.psevere$design == "IncRate"]),
                           c(GA.sugguestion[6,], optimal.psevere$p_severe[optimal.psevere$design == "Case_Severe"]),
                           c(GA.sugguestion[7,], optimal.psevere$p_severe[optimal.psevere$design == "IncRate_Severe"])
)


#--------- 2009-2014 -----------
load("Result\\Disease model\\Disease_model_realizations_6.RData")
result.all <- data.frame(Type = rep(Alternative.name, each = rep.times), ofv = NA)
for(k in 1:8)
{
  result.all$ofv[result.all$Type == Alternative.name[k]] <- foreach(i = 1:rep.times, .combine = "c",  .packages = c("rstan", "abind")) %dopar% {
    cat(Alternative.name[k], i, "\n")
    obj.fun(Alternative.design[[k]], dm.data = dm.realization, test.total = dat.agg$Z[1:6], goal = "severe")
  }
}
save(result.all, file = "Result\\Compare_OFV_of_benchmark_severe_09to14.Rdata")


#-------------- 2015 ---------------
load("Result\\Disease model\\Disease_model_realizations_7.RData")
result.all <- data.frame(Type = rep(Alternative.name, each = rep.times), ofv = NA)
for(k in 1:8)
{
  result.all$ofv[result.all$Type == Alternative.name[k]] <- foreach(i = 1:rep.times, .combine = "c",  .packages = c("rstan", "abind")) %dopar% {
    cat(Alternative.name[k], i, "\n")
    inter.result <- obj.fun(Alternative.design[[k]], dm.data = dm.realization, test.total = dat.agg$Z, goal = "severe", raw = TRUE)
    
    mean(abs(inter.result[[1]][,,7] - inter.result[[2]][,,7]))
  }
}
save(result.all, file = "Result\\Compare_OFV_of_benchmark_severe_15.Rdata")




stopCluster(cl)
