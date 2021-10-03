pmod = stan_model('Stan code\\stanmodel_pars_as_data.stan', auto_write = TRUE)
options(mc.cores = 1)
rstan::rstan_options(auto_write = TRUE)


#' Converting an n-1 d unconstrained parameters to an n d sum-to-one parameters
#'
#' @param n the number of constrained parameters
#' @param x the input
#'
#' @return  an n dimension array
#' @export
#'
#' @examples
convert.var <- function(n, x)
{
  xout <- array(NA, n)
  
  xout[1] <- 1 - x[1]^(1/(n-1))
  
  for(i in 2:(n-1))
  {
    xout[i] <- (1 - sum(xout[1:(i-1)]))*(1 - x[i]^(1/(n - i)))
  }
  
  xout[n] = 1 - sum(xout[1:(n-1)])
  
  xout
}



#' Converting an n d sum-to-one parameters to an n-1 d unconstrained parameters
#'
#' @param n the number of constrained parameters
#' @param xout the input (n d parameters)
#'
#' @return  an n dimension array
#' @export
#'
#' @examples
convert.inverse <- function(n, xout)
{
  x <- array(NA, n-1)
  
  x[1]  <- (1 - xout[1])^(n-1)
  
  for(i in 2:(n-1))
  {
    x[i] <- (1 - xout[i]/(1 - sum(xout[1:(i-1)])))^(n - i)
  }
  
  x
}







#' Evaluating the objective function
#'
#' @param unconstrained.design   the 20-dimension unconstrained design
#' @param dm.data                the disease model data
#' @param test.total             the total number of tests
#' @param goal                   if the surveillance goal is minimizing MAE in the estimation of incidence rate or composition
#'
#' @return
#' @export
#'
#' @examples
obj.fun <- function(unconstrained.design, dm.data, test.total, goal = c("all", "severe"), raw = FALSE, i = sample(1:4000, 1))
{
  if(!goal %in% c("all", "severe"))
    stop("goal should be either all or severe!")
  
  unconstrained.design <- unlist(unconstrained.design)
  
  # convert the unconstrained design parameter to meet the sum-to-one constraint
  constrained.design <- convert.var(n = 21, x = unconstrained.design[1:20])
  p.severe <- unconstrained.design[21]
  
  pdata = list(Ni = dm.data$Ni, Nk = dm.data$Nk, Nt = dm.data$Nt, D = dm.data$D, W = dm.data$W, N_it = dm.data$N_it, Y = dm.data$Y, Y_severe = dm.data$Y_severe) # format input for the surveillance model
  
  test.spatial <- sapply(test.total, function(x) round(x*constrained.design))  # estimate the number of test allocate to each spatial unit determined by the current design parameter, each row represent one year
  
  test.spatial[test.spatial > pdata$Y] <- pdata$Y[test.spatial > pdata$Y]  # can't higher than the total number of cases
  
  pdata$p_test_severe <- matrix(p.severe, nrow(pdata$Y_severe), ncol(pdata$Y_severe))
  pdata$p_test_severe[pdata$Y_severe*p.severe > test.spatial] <- test.spatial[pdata$Y_severe*p.severe > test.spatial]/pdata$Y_severe[pdata$Y_severe*p.severe > test.spatial]
  
  test_severe <- round(pdata$Y_severe*p.severe)    # assuming all severe cases get tested
  test_severe[test_severe > test.spatial] <- test.spatial[test_severe > test.spatial]
  test_mild <- test.spatial - test_severe
  
  pdata$p_test_mild <- test_mild/(pdata$Y - pdata$Y_severe)
  pdata$p_test_mild[pdata$p_test_mild == 0] <- 1e-8
  
  # generate test results as surveillance system realizations
  # attach parameters to pdata
  pdata$G_t_L = dm.data$G_t_L[i,,]
  pdata$sigma_k = dm.data$sigma_k[i,]
  pdata$rho = dm.data$rho[i]
  pdata$alpha = dm.data$alpha[i]
  pdata$beta0 = dm.data$beta0[i]
  
  lambda_it <- apply(dm.data$lambda_ikt[i,,,], c(1, 3), sum)
  p_ikt <- sweep(dm.data$lambda_ikt[i,,,], c(1,3), lambda_it, "/")
  
  p_severe_k <- dm.data$p_severe_k[i,]
  lambda_ikt_severe <- sweep(dm.data$lambda_ikt[i,,,], 2, dm.data$p_severe_k[i, ], "*")
  lambda_ikt_mild <- sweep(dm.data$lambda_ikt[i,,,], 2, (1 - dm.data$p_severe_k[i, ]), "*")
  
  
  sample.from.mild <- abind::abind(lambda_ikt_mild, test_mild, along = 2) # append pathogen-specific incidence rate with the total test number
  sample.from.severe <- abind::abind(lambda_ikt_severe, test_severe, along = 2)
  
  current.realization.mild <- apply(sample.from.mild, c(1,3), function(x) rmultinom(1, x[4], prob = x[1:3]))
  current.realization.severe <- apply(sample.from.severe, c(1,3), function(x) rmultinom(1, x[4], prob = x[1:3]))
  
  # format data for refit
  pdata$Z_mild <- aperm(current.realization.mild, c(2,1,3))
  pdata$Z_severe <- aperm(current.realization.severe, c(2,1,3))
  
  pfit <- rstan::sampling(pmod, data = pdata, chains = 1, refresh = 0, iter = 500, pars = c("lambda_ikt", "p_severe_k"), warmup = 200, cores = 1)
  
  if(length(pfit@sim) == 0)
  {
    if(!raw)
    {
      ofv = -10
      return(ofv)
    } else
    {
      return(NA)
    }
  } else {
    # calculate objective function value
    lambda_ikt.estimate <- apply(rstan::extract(pfit,'lambda_ikt')[[1]], c(2,3,4), median)
    
    p_severe_k.est <- apply(rstan::extract(pfit, "p_severe_k")[[1]], 2, median)
    
    
    if(!raw)
    {

      
      if(goal == "all")
      {
        ofv <- -mean(abs(dm.data$lambda_ikt[i,,,] - lambda_ikt.estimate))*10   # unit: per 100,000; negative since GA can only maximize
      }

      if(goal == "severe")
      {
        ofv <- -mean(abs(sweep(dm.data$lambda_ikt[i,,,], 2, dm.data$p_severe_k[i,], "*") - sweep(lambda_ikt.estimate, 2, p_severe_k.est, "*")))*10   # unit: per 100,000; negative since GA can only maximize
      }
      
      return(ofv)
    }
    
    if(raw)
    {
      if(goal == "all")
      {
        raw.data <- list(dm.data$lambda_ikt[i,,,]*10, lambda_ikt.estimate*10)
      }
      
      if(goal == "severe")
      {
        raw.data <- list(sweep(dm.data$lambda_ikt[i,,,], 2, dm.data$p_severe_k[i,], "*")*10, sweep(lambda_ikt.estimate, 2, p_severe_k.est, "*")*10)
      }
      return(raw.data)
    }
  }
}








#' Function for evaluating the obj fn of multiple disease realizations in parallel
#'
#' @param design.par       the design parameter
#' @param dm.data          the disease model data
#' @param test.total       the total number of tests
#' @param goal             if the surveillance goal is minimizing MAE in the estimation of incidence rate or composition
#' @param repi             the number of realizations to sample
#'
#' @return
#' @export
#'
#' @examples
obj.fun.par <- function(design.par, dm.data, test.total, goal = c("all", "severe"), repi)
{
    all.obj.value <- foreach(i = 1:repi, .combine = "c",  .packages = c("rstan", "abind")) %dopar% {
      obj.fun(design.par, dm.data, test.total, goal)
    }
    
    mean(all.obj.value, na.rm = TRUE)
}