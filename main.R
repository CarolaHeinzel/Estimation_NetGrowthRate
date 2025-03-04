# This code calculates the different estimators for the net growth rate and compares them to each other.

library(gtools)
library(pracma) 
library(cloneRate)
library(Rmpfr)
source("new_Estimator.R")
source("Simulation_CPP.R")

# Important to calculate the expected value according to Durrett
ew_durrett <- function(n){
  sum_result <- 0
  for (k in 2:(n - 1)) {
    value <-  n / (k*(k-1))
    sum_result <- sum_result + value 
  }
  return(sum_result)
}

compute_list_properties <- function(liste, n) {
  normal_list = liste
  max_value <- max(normal_list)
  mean_value <- mean(normal_list)
  
  # ∑_i ∑_j (H_i - H_j)^+
  double_sum <- sum(outer(normal_list, normal_list, function(x, y) pmax(0, x - y)))
  
  return(list(Maximum = max_value, Mittelwert = mean_value, d_sum = double_sum))
}

# Calculates \hat r_{new,2}
calc_estimator <- function(n, h){
  values <- compute_list_properties(h, n)
  est <-  1/(n-1)*values$Doppelsumme + values$Maximum - values$Mittelwert 
  ev_th <- ew_durrett(n)
  return(ev_th/est)
  
}

# Calculates \hat r_{inverse}, \hat r_{JS, MSE}
calc_estimator_adapted <- function(n, h){
  if(n==5){
    c_n <- 12*115/144
  }else if(n == 6){
    c_n <- 5*4*163/200
  }else if(n == 7){
    c_n <- 6*5*497/600
  }else if(n == 8){
    c_n <- 7*6*617/735
  }else if(n == 9){
    c_n <- 8*7*13311/15680
  }else if(n == 10){
    c_n <- 9*8*15551/18144
  }
  values <- compute_list_properties(h, n)
  est <-  values$d_sum
  return(list(c_n/est, 1/est))
  
}

# repeat the simulation num_Rep times
repeat_simulation <- function(num_Rep, n, r, t){
  estimator <- numeric(num_Rep)
  estimator_mcmc <- c(num_Rep)
  estimator_new <- c(num_Rep)
  estimator_JS <- c(num_Rep)
  estimator_JS_adapted <- c(num_Rep)
  ew_durrett <- ew_durrett(n)
  ew_CH <- c(num_Rep)
  set.seed(n)
  for(i in 1:num_Rep){
    a = runif(1, min = r, max = r+1)
    test <- simUltra(a,b = a - r,cloneAge = t,n = n,nTrees = 1,precBits = 1000,addStem = FALSE,nCores = 1)
    res <- cloneRate::internalLengths(test[[1]], alpha = 0.05)
    il= res$sumInternalLengths
    h = test[[2]]
    test_mcmc <-birthDeathMCMC(test[[1]],maxGrowthRate = 4,alpha = 0.05,verbose = TRUE,nChains = 4,nCores = 1,chainLength = 2000 )
    estimator_JS[i] = calc_estimator(n, h)
    estimator_JS_adapted[i] = calc_estimator_adapted(n, h)[[1]]
    res_new =  internalLengths(test[[1]], alpha = 0.05)
    ew_new <- ew_durrett_adapted(n, res_new[[2]])
    estimator[i] <-  ew_durrett/il # durrett
    estimator_mcmc[i] <- test_mcmc$estimate
    estimator_new[i] <- ew_new/il
    ew_CH[i] = calc_estimator_adapted(n, h)[[2]]
  }
  # durrett, Pyhloft, \hat r_{new}, \hat r_{JS}, \hat r_{new, 2}, \hat r_{JS, MSE}
  return(list(estimator, estimator_mcmc, estimator_new, estimator_JS, estimator_JS_adapted, ew_CH))
}

# repeat the simulation for different r
# save the estimators in a .txt file
repeat_simulation_T <- function(n_rep, t, r_list, n){
  K = length(r_list)
  results_repeat <- rep(0, K)
  estimator <- numeric(K)
  estimator_old <- c(K)
  estimator_mcmc <- c(K)
  estimator_new <- c(K)
  estimator_JS <- c(K)
  estimator_CH <- c(K)
  i = 1
  for(r in r_list){
    temp <- repeat_simulation(n_rep, n, r, t)
    temp_df <- data.frame(
      list1 = temp[[1]],                    
      list2 = temp[[2]],                   
      list3 = temp[[3]],                    
      list4 = temp[[4]],                    
      list5 = temp[[5]],
      list6 = temp[[6]]
    )
    # in list6, we have to include the prefactor afterwards
    output_file = paste0("table_test_t_",t,"_r_",r,"_n_",n , ".txt")
    write.table(
      temp_df,
      file = output_file,
      sep = "\t",          
      row.names = FALSE,    
      col.names = TRUE,      
    )
    mse1 = mean((temp[[1]]-r)^2)
    mse2 = mean((temp[[2]]-r)^2)
    mse3 = mean((temp[[3]]-r)^2)
    mse4 = mean((temp[[4]]-r)^2)
    mse5 = mean((temp[[5]]-r)^2)
    if(n==5){
      c_CH = 4.769173
    }else if(n==6){
      c_CH = 9.833995
    }else if(n==7){
      c_CH = 17.367272
    }else if(n==8){
      c_CH = 26.940960
    }else if(n==9){
      c_CH = 36.827248
    }else if(n==10){
      c_CH = 51.135133
    }
    mse6 = mean((temp[[6]]*c_CH-r)^2)
    estimator_old[i] = mse1
    estimator[i] = mse2
    estimator_mcmc[i] = mse3
    estimator_new[i] = mse4
    estimator_JS[i] = mse5
    estimator_CH[i] = mse6
    i = i+1
  }
  return(list(estimator_old, estimator, estimator_mcmc, estimator_new, estimator_JS, estimator_CH))
}
# Example Usage
temp = repeat_simulation_T(100, 40, c(0.5, 1, 2), 5) 
