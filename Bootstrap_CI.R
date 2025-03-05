source("main.R")

# Code to evaluate whether the estimator is asym
asymptotic_distribution <- function(numRep, n, t, r){
  estimator <- rep(0, numRep)
  for(i in 1:numRep){
    a = runif(1, min = r, max = r+1)
    test <- simUltra(a,b = a - r,cloneAge = t,n = n,nTrees = 1,precBits = 1000,addStem = FALSE,nCores = 1)
    res <- cloneRate::internalLengths(test[[1]], alpha = 0.05)
    il= res$sumInternalLengths
    h = test[[2]]
    hat_r = calc_estimator_adapted(n, h)
    print(il)
    print(hat_r)
    estimator[i] <- hat_r[[2]]/il
  }
  return(estimator)
}

r = 1
n = 1000
t = 40
numRep = 100
res_normal = asymptotic_distribution(numRep, n, t, r)


# QQ-Plot 
qqnorm(res_normal)  
qqline(res_normal, col = "red", lwd = 2)  


shapiro_test <- shapiro.test(res_normal)
print(shapiro_test)

ks_test <- ks.test(res_normal, "pnorm", mean = mean(res_normal), sd = sd(res_normal))
print(ks_test)
