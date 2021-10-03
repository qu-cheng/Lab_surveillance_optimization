#=================== yearly GP model===========================
library(rstan)
library(tidyverse)

# common stan settings
setMKLthreads(1)
options(mc.cores = 4)
rstan_options(auto_write = TRUE)

#######################################################################
#              Input data and data processing
#######################################################################
dat <- read_csv(".\\Data\\Sichuan_pref_yearly.csv") %>% as.data.frame()   # Surveillance data
dat.6 <- dat %>%                                       # filter to data before 2015
  filter(Year < 2015)

D <- read.csv(".\\Data\\D_21.csv")    # matrices reflecting the spatial adjacencies
W <- read.csv(".\\Data\\W_21.csv")    # matrices reflecting the spatial adjacencies

# reshape testing data
T = length(unique(dat$Year))
S = length(unique(dat$Pref))
K = 3
Z_mild_array <- Z_severe_array <- array(NA, dim = c(S, K, T))

for(k in 1:K){
  Z_mild_array[,k,] <- as.array(matrix(dat[,k+5], ncol = T, byrow = TRUE))
  Z_severe_array[,k,] <- as.array(matrix(dat[,k+8], ncol = T, byrow = TRUE))
}


#######################################################################
#                  Smooth testing probablity
#######################################################################
# stan model for smoothing probability of testing mild cases so it will be never 0, this smoothed probability will be used in the disease model 
testmod = stan_model('Stan code\\p_test_model.stan')

#======================= for mild cases ==============================
# use all 7 years of data
tdata = list(
  Ni = S, Nt = T,
  Y = as.array(matrix(dat$Y_mild, ncol = T, byrow = TRUE)),
  Z = as.array(matrix(dat$Z_mild_1 + dat$Z_mild_2 + dat$Z_mild_3, ncol = T, byrow = TRUE))
)
mtest_fit_7 = sampling(testmod, data = tdata)
save(mtest_fit_7, file ="Result\\Disease model\\Mild_testing_prob_smooth_7.RData")

# use the first 6 years of data
tdata = list(
  Ni = S, Nt = 6,
  Y = as.array(matrix(dat$Y_mild, ncol = T, byrow = TRUE)[, 1:6]),
  Z = as.array(matrix(dat$Z_mild_1 + dat$Z_mild_2 + dat$Z_mild_3, ncol = T, byrow = TRUE)[, 1:6])
)
mtest_fit_6 = sampling(testmod, data = tdata)
save(mtest_fit_6, file ="Result\\Disease model\\Mild_testing_prob_smooth_6.RData")



#======================= for severe cases ============================
# with 7 years of data
tdata = list(
  Ni = S, Nt = T,
  Y = as.array(matrix(dat$Y_severe, ncol = T, byrow = TRUE)),
  Z = as.array(matrix(dat$Z_severe_1 + dat$Z_severe_2 + dat$Z_severe_3, ncol = T, byrow = TRUE))
)
stest_fit_7 = sampling(testmod, data = tdata)
save(stest_fit_7, file ="Result\\Disease model\\Severe_testing_prob_smooth_7.RData")

# with 6 years of data
tdata = list(
  Ni = S, Nt = 6,
  Y = as.array(matrix(dat$Y_severe, ncol = T, byrow = TRUE)[, 1:6]),
  Z = as.array(matrix(dat$Z_severe_1 + dat$Z_severe_2 + dat$Z_severe_3, ncol = T, byrow = TRUE)[, 1:6])
)
stest_fit_6 = sampling(testmod, data = tdata)
save(stest_fit_6, file ="Result\\Disease model\\Severe_testing_prob_smooth_6.RData")



#######################################################################
#                  disease model fit
#######################################################################
stan.model <- stan_model('Stan code\\stanmodel.stan')

# all 7 years of data
fdata = list(
  Ni = S, Nk = K, Nt = T, D = D, W = W, 
  N_it = matrix(dat$pop, ncol = T, byrow = T)/1e4,
  Y = as.array(matrix(dat$Y, ncol = T, byrow = TRUE)),
  Y_severe = as.array(matrix(dat$Y_severe, ncol = T, byrow = TRUE)),
  Z_mild = aperm(array(c(dat$Z_mild_1,dat$Z_mild_2,
                         dat$Z_mild_3),dim = c(T,S,K)),
                 c(2,3,1)),
  Z_severe = aperm(array(c(dat$Z_severe_1,dat$Z_severe_2,
                           dat$Z_severe_3),dim = c(T,S,K)),
                   c(2,3,1)),
  p_test_severe = matrix(summary(stest_fit_7,'theta_it')$summary[,6],
                         ncol = T, byrow = T),
  p_test_mild = matrix(summary(mtest_fit_7,'theta_it')$summary[,6],
                       ncol = T, byrow = T)
)

fit <- sampling(stan.model, data = fdata,chains = 4, refresh = 10, iter = 2000)
save(fit, file ="Result\\Disease model\\Disease_model_fit_7.RData")

# 6 years of data
T.new = 6
fdata = list(
  Ni = S, Nk = K, Nt = T.new, D = D, W = W, 
  N_it = matrix(dat.6$pop, ncol = T.new, byrow = TRUE)/1e4,
  Y = as.array(matrix(dat.6$Y, ncol = T.new, byrow = TRUE)),
  Y_severe = as.array(matrix(dat.6$Y_severe, ncol = T.new, byrow = TRUE)),
  Z_mild = aperm(array(c(dat.6$Z_mild_1,dat.6$Z_mild_2,
                         dat.6$Z_mild_3),dim = c(T.new,S,K)),
                 c(2,3,1)),
  Z_severe = aperm(array(c(dat.6$Z_severe_1,dat.6$Z_severe_2,
                           dat.6$Z_severe_3),dim = c(T.new,S,K)),
                   c(2,3,1)),
  p_test_severe = matrix(summary(stest_fit_6,'theta_it')$summary[,6],
                         ncol = T.new, byrow = T),
  p_test_mild = matrix(summary(mtest_fit_6,'theta_it')$summary[,6],
                       ncol = T.new, byrow = T)
)

fit <- sampling(stan.model, data = fdata,chains = 4, refresh = 10, iter = 2000)
save(fit, file ="Result\\Disease model\\Disease_model_fit_6.RData")