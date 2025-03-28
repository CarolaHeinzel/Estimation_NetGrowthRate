library(gtools)
library(pracma) 
library(cloneRate)
library(Rmpfr)
library(reshape2) 


# Calculates \hat r_{inverse}, \hat r_{JS, MSE}
calc_estimator_adapted <- function(n, h){
  values <- compute_list_properties(h, n)
  est <-  values
  return(1/est)
  
}

compute_list_properties <- function(liste, n) {
  normal_list = liste
  # ∑_i ∑_j (H_i - H_j)^+
  double_sum <- sum(outer(normal_list, normal_list, function(x, y) pmax(0, x - y)))
  return(double_sum)
}

# repeat the simulation num_Rep times
repeat_simulation_dist <- function(num_Rep, n, r, t){
  inverse_s <- c(num_Rep)

  #set.seed(n+1)
  for(i in 1:num_Rep){
    a = runif(1, min = r, max = r+1)
    test <- simUltra(a,b = a - r,cloneAge = t,n = n,nTrees = 1,precBits = 1000,addStem = FALSE,nCores = 1)
    
    res <- cloneRate::internalLengths(test[[1]], alpha = 0.05)
    il= res$sumInternalLengths
    h = test[[2]]
    # New estimator
    inverse_s[i] = calc_estimator_adapted(n, h)
  }
  # Yubo, Phylofit, New estimator
  return(inverse_s)
}

t_new = repeat_simulation_dist(100, 10, 0.5, 40)
t_new = repeat_simulation_dist(100, 10, 1, 40)
t_new = repeat_simulation_dist(100, 10, 0.5, 100)

#  Plot the distribution

density_values <- density(t_new)

# Plotten der Dichte
plot(density_values, 
     main = "Density n = 10, r = 0.5, T = 100", 
     xlab = "Values", 
     ylab = "Density", 
     col = "blue", 
     lwd = 2)
