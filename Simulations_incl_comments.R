library(cloneRate)

# Expected Value new
ew_func <- function(n){
  sum_result <- 0
  for (k in 2:(n - 1)) {
    for (i in 1:(n-k-1)){
      value <-   (k+n) / ((k - 1) * k*(n+1)) - 2/(k*(n+1))
      sum_result <- sum_result + value 
    }
    sum_result <- sum_result + 2*n/(k*(n+1))
  }
  return(sum_result)
}
# Expected Value Durrett
ew_durrett <- function(n){
  sum_result <- 0
  for (k in 2:(n - 1)) {
    value <-  n / (k*(k-1))
    sum_result <- sum_result + value 
  }
  return(sum_result)
}

# Simulate the Data
repeat_simulation <- function(num_Rep, n){
  estimator <- numeric(num_Rep)
  estimator_old <- c(num_Rep)
  estimator_durrett <- c(num_Rep)
  ew <- ew_func(n)
  ew_durrett <- ew_durrett(n)
  
  for(i in 1:num_Rep){
    r = 0.5
    a = runif(1, min = 0.5, max = 1.5)
    test <- simUltra(a,b = a - r,cloneAge = 40,n = n,nTrees = 1,precBits = 1000,addStem = FALSE,nCores = 1)
    res <- internalLengths(test, alpha = 0.05)
    il= res$sumInternalLengths
    test_mcmc <- birthDeathMCMC(
      test,
      maxGrowthRate = 4,
      alpha = 0.05,
      verbose = TRUE,
      nChains = 1,
      nCores = 1,
      chainLength = 2000
    )
    estimator[i] <- ew/il
    estimator_old[i] <- test_mcmc$estimate # Phlyofit
    #res$estimate would be the old one based on the lengths
    estimator_durrett[i] <- ew_durrett/il
  }
  return(list(estimator_old, estimator, estimator_durrett))
}

res_all <- numeric(5)
res_all_old <- numeric(5)
res_p <- numeric(5)

n_all <- c(5,6,7,8,9)
j = 1
# Only works for r = 0.5
for(i in n_all){
  res_test <- repeat_simulation(100, i)
  x1 <- res_test[[1]]
  x2 <- res_test[[2]]
  x3 <- res_test[[3]]
  mse_old <- mean((x1 - 0.5)^2)
  mse_new <- mean((x2 - 0.5)^2)
  mse_p <- mean((x3-0.5)^2)
  res_all[j] <- mse_new
  res_all_old[j] <- mse_old
  res_p[j] <- mse_p
  j <- j+1
}

x <- n_all 

# Create Plot 
res_p = c( 0.48243283, 0.29230783, 0.16264230, 0.10191675, 0.06320459) # Phylofit
res_all = c(0.12922632,0.06071594, 0.04398088, 0.10730434, 0.02440512) # New
plot(x, res_all, type = "b", col = "blue", pch = 19, ylim = range(c(res_all, res_p)),
     xlab = "n", ylab = "MSE")
res_durrett =  c(0.25710729, 0.10272845, 0.06690614 ,0.16368736, 0.03056126) #Durrett
lines(x, res_p, type = "b", col = "green", pch = 17)
lines(x, res_durrett, type = "b", col = "red", pch = 17)

legend("topright", legend = c("Lengths small sample sizes", "Durrett", "Phylofit"),
       col = c("blue", "red", "green"), pch = c(19, 17, 17), lty = 1)
