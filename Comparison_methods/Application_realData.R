library(cloneRate)
library(ape)
library(ggplot2)
data(realCloneData)
data(longitudinalData)

for(i in 1:42){
  max_temp = Ntip(realCloneData[["cloneTrees"]][[i]])
  max_individuals[i] = max_temp # number of individuals per sample
}

# plot n
hist(max_individuals,
     xlab = "n",
     ylab = "Frequency",
     col = "skyblue",
     breaks = 22,
     xlim = c(0, max(max_individuals) + 20),
     main = "")


# Calculate Coalescence Times
get_coal_times <- function(tree) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required for this function.")
  }
  coal_times <- ape::branching.times(tree)
  return(coal_times)
}

# Apply the different estimation methods to the data
calc_estimator_realdata <- function(numInds){
  estimator_JS = rep(0, numInds) # novel method
  estimator_mcmc = rep(0, numInds) # method by Tanja Stadler
  
  for(i in 1:numInds){
    tree = realCloneData[["cloneTrees"]][[i]] # Real Data
    res = cloneRate::internalLengths(tree, alpha = 0.05)
    n = Ntip(tree)
    il= res$sumInternalLengths
    coal_times = node.depth.edgelength(tree)[]
    h = get_coal_times(tree)
    # New estimator
    estimator_JS[i] = calc_estimator_adapted(n, h)
    test_mcmc <-birthDeathMCMC(tree,maxGrowthRate = 4,alpha = 0.05,verbose = FALSE,nChains = 4,nCores = 1,chainLength = 2000 )
    estimator_mcmc[i] <- test_mcmc$estimate
  }
  return(list(estimator_mcmc, estimator_JS ))
}

est_real_try1 = calc_estimator_realdata(42)

# Calculate the estimators

# Table with n and the constants c_MSE(n)
df <- data.frame(
  n = c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
        30, 40, 50, 60, 70, 80, 90, 100, 23, 27, 29, 34, 48, 46, 53, 76, 78, 71, 83, 109
        ,22, 25, 65, 51, 38),
  value3 = c(0.355, 0.49, 0.58, 0.63, 0.67, 0.70, 0.80, 0.81, 0.82,
             0.83, 0.84, 0.84, 0.85, 0.86, 0.86, 0.87,
             0.92, 0.93, 0.93, 0.95, 0.96, 0.96, 0.96, 0.96,
             0.83, 0.86, 0.87, 0.88, 0.91, 0.89, 0.92, 
             0.94, 0.94, 0.94,0.94, 0.96, 0.88, 0.90, 0.94, 0.93, 0.92)
)

n_values <- max_individuals

lookup_values <- df$value3[match(n_values, df$n)]

result <- lookup_values * (n_values-2) * (n_values - 1)

result_all = rep(0, 42)
for(i in 1:42){
  result_all[i] = result[i] * est_real_try1[2][[1]][i]
}

# Plot the results
est_real_try_copy = est_real_try1
est_real_try_copy[[2]] = result_all

# Sort the individuals according to n
ord <- order(max_individuals)   
est_real_try_copy_c1 <- est_real_try_copy[[1]][ord]    
est_real_try_copy_c2 <- est_real_try_copy[[2]][ord]    
final =  est_real_try_copy    #
final[[1]] = est_real_try_copy_c1
final[[2]] = est_real_try_copy_c2


names_final = names(realCloneData[["cloneTrees"]])[ord]
n_ind_sorted <- max_individuals[ord]

res_sorted = cbind(est_real_try_copy_c2, est_real_try_copy_c1)


colors <- c("#e41a1c", "#377eb8")

par(mar = c(8, 4, 4, 2) + 0.1)
bp <- barplot(
  t(res_sorted),
  beside = TRUE,
  col = colors,
  xaxt = "n",           
  xlab = "",
  ylab = "Estimator"
)

legend(
  "topright",
  legend = c(
    expression(hat(r)[MSE]),
    expression(hat(r)[MCMC])),
  fill = colors
)

x_labels <- paste0(names_final, " (", n_ind_sorted, ")")

axis(
  1,                  
  at = colMeans(bp),     
  labels = x_labels,
  las = 2,             
  cex.axis = 0.7         
)

﻿
