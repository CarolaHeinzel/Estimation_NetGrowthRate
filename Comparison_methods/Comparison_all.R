library(gtools)
library(pracma) 
library(cloneRate)
library(Rmpfr)
library(reshape2) 
source("Simulation_self.R")

compute_list_properties <- function(liste, n) {
  normal_list = liste
  # ∑_i ∑_j (H_i - H_j)^+
  double_sum <- sum(outer(normal_list, normal_list, function(x, y) pmax(0, x - y)))
  return(double_sum)
}

calc_estimator_adapted <- function(n, h){
  values <- compute_list_properties(h, n)
  return(1/values)
}

# repeat the simulation num_Rep times
# t is the sampling time, r the true net growth rate and n the sampling size.
repeat_simulation <- function(num_Rep, n, r, t){
  estimator_mcmc <- c(num_Rep)
  estimator_old <- c(num_Rep)
  estimator_JS <- c(num_Rep)
  set.seed(r+1)
  for(i in 1:num_Rep){
    print(i)
    a = runif(1, min = r, max = r+1)
    test <- simUltra(a,b = a - r,cloneAge = t,n = n,nTrees = 1,precBits = 1000,addStem = FALSE,nCores = 1)
    
    res <- cloneRate::internalLengths(test[[1]], alpha = 0.05)
    il= res$sumInternalLengths
    h = test[[2]]
    
    # New estimator
    estimator_JS[i] = calc_estimator_adapted(n, h)#[[1]]
    
    # Estimator based on Johnson et al.
    estimator_old[i] = n/il
    
    # Method by Tanja Stadler
    test_mcmc <-birthDeathMCMC(test[[1]],maxGrowthRate = 4,alpha = 0.05,verbose = TRUE,nChains = 4,nCores = 1,chainLength = 2000 )
    estimator_mcmc[i] <- test_mcmc$estimate
  }
  filename <- paste0("CPPs_eva_2_1706_test_",t,"_n_", n, "_r_", r, ".txt")
  est_all = cbind(estimator_JS, estimator_old, estimator_mcmc)
  write.table(est_all, file = filename, row.names = FALSE, col.names = FALSE)
  
  # Johnson, Method by Tanja Stadler, New estimator (without the constant)
  return(list(estimator_old, estimator_mcmc, estimator_JS))
}

est = repeat_simulation(2, 5, 1, 40)
