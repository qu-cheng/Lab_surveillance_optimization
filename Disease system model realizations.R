library(rstan)

#######################################################################
#              Input data and data processing
#######################################################################
dat <- read_csv(".\\Data\\Sichuan_pref_yearly.csv") %>% as.data.frame()   # Surveillance data
dat.6 <- dat %>%                                       # filter to data before 2015
  filter(Year < 2015)

D <- read.csv(".\\Data\\D_21.csv")    # matrices reflecting the spatial adjacencies
W <- read.csv(".\\Data\\W_21.csv")    # matrices reflecting the spatial adjacencies


#######################################################################
#                    use all 7 years of data
#######################################################################
load("Result\\Disease model\\Disease_model_fit_7.RData")   # output of Disease system model fit.R

T = length(unique(dat$Year))
S = length(unique(dat$Pref))
K = 3


# set.seed(20201002)
sampleID <- 1:4000
dm.realization <- list(Ni = S, Nk = K, Nt = T, D = D, W = W,
                       N_it = matrix(dat$pop, ncol = T, byrow = T)/1e4,
                       Y = as.array(matrix(dat$Y, ncol = T, byrow = TRUE)),
                       Y_severe = as.array(matrix(dat$Y_severe, ncol = T, byrow = TRUE)))


dm.realization$beta0 <- rstan::extract(fit, 'beta0')[[1]]
dm.realization$G_t_L <- rstan::extract(fit,'G_t_L')[[1]][sampleID, ,]
dm.realization$sigma_k <- rstan::extract(fit,'sigma_k')[[1]][sampleID, ]
dm.realization$rho <- rstan::extract(fit,'rho')[[1]][sampleID]
dm.realization$alpha <- rstan::extract(fit,'alpha')[[1]][sampleID]
dm.realization$p_severe_k <- rstan::extract(fit,'p_severe_k')[[1]][sampleID, ]
dm.realization$lambda_ikt <- rstan::extract(fit,'lambda_ikt')[[1]][sampleID, , ,]

save(dm.realization, file ="Result\\Disease model\\Disease_model_realizations_7.RData")



#######################################################################
#                    use first 6 years of data
#######################################################################
load("Result\\Disease model\\Disease_model_fit_6.RData")   # output of Disease system model fit.R

T = length(unique(dat.6$Year))
S = length(unique(dat.6$Pref))
K = 3


# set.seed(20201002)
sampleID <- 1:4000

dm.realization <- list(Ni = S, Nk = K, Nt = T, D = D, W = W, 
                       N_it = matrix(dat.6$pop, ncol = 6, byrow = T)/1e4,
                       Y = as.array(matrix(dat.6$Y, ncol = 6, byrow = TRUE)),
                       Y_severe = as.array(matrix(dat.6$Y_severe, ncol = 6, byrow = TRUE)))

dm.realization$beta0 <- rstan::extract(fit, 'beta0')[[1]]
dm.realization$G_t_L <- rstan::extract(fit,'G_t_L')[[1]][sampleID, ,]
dm.realization$sigma_k <- rstan::extract(fit,'sigma_k')[[1]][sampleID, ]
dm.realization$rho <- rstan::extract(fit,'rho')[[1]][sampleID]
dm.realization$alpha <- rstan::extract(fit,'alpha')[[1]][sampleID]
dm.realization$p_severe_k <- rstan::extract(fit,'p_severe_k')[[1]][sampleID, ]
dm.realization$lambda_ikt <- rstan::extract(fit,'lambda_ikt')[[1]][sampleID, , ,]

save(dm.realization, file ="Result\\Disease model\\Disease_model_realizations_6.RData")
