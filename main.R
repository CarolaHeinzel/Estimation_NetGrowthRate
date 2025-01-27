library(gtools) # Für Permutationen
library(pracma) # Für die Factorial-Funktion

# Werte initialisieren
h_values <- c(1, 2, 3, 4)

# Alle Permutationen

is_valid_permutation <- function(h, k) {
  if (k > 1) {
    #print(h[k])
    #print(h[1:(k-1)])
    if (h[k] > max(h[1:(k-1)])) {
      return(FALSE)
    }
    if(n-k-1 > 0){
      for (i in 1:(n - k-1)) {
        if (min(h[i], h[i + k]) > max(h[(i + 1):(i + k-1)])) {
          return(FALSE)
        }
      }
    }
    #print(h[n-k])
    #print(h[(n - k+1):(n-1)])
    if (h[n - k] > max(h[(n - k+1):(n-1)])) {
      return(FALSE)
    }
  }
  return(TRUE)
}


k <- 3
h_values = seq(1, 7)
n <- length(h_values) +1
permutations <- permutations(n = length(h_values), r = length(h_values), v = h_values)

# Filtere gültige Permutationen
valid_permutations <- permutations[apply(permutations, 1, function(perm) is_valid_permutation(perm, k)), ]
print(valid_permutations)
ratio <- nrow(valid_permutations) / factorial(n-1)
cat(sprintf("Ratio of Valid Permutations: %.4f\n", ratio))


ew_durrett_adapted <- function(n, tree_data){
  sum_result <- 0
  temp = calculate_leaf_info(tree_data, n)
  h_values <-  seq(1, (n-1))
  permutations <- permutations(n = length(h_values), r = length(h_values), v = h_values)
  for (k in 1:(n - 2)) {
    if(temp$check_leaves$Exists[k] == TRUE){
      #print(k)
      #print(h_values)
      permutations <- permutations(n = length(h_values), r = length(h_values), v = h_values)
      #print(permutations)
      valid_permutations <- permutations[apply(permutations, 1, function(perm) is_valid_permutation(perm, k+1)), ]
      # k = 3 falls k = 2 in R
      den = nrow(valid_permutations)
      en = factorial(n-1)
      value <-  n / (k*(k+1)) * en/(en - den)
      sum_result <- sum_result + value
    }
  }
  return(sum_result)
}


repeat_simulation <- function(num_Rep, n, r, t){
  estimator <- numeric(num_Rep)
  estimator_old <- c(num_Rep)
  estimator_mcmc <- c(num_Rep)
  estimator_new <- c(num_Rep)
  estimator_JS <- c(num_Rep)
  ew_durrett <- ew_durrett(n)
  for(i in 1:num_Rep){
    a = runif(1, min = r, max = r+1)
    test <- simUltra(a,b = a - r,cloneAge = t,n = n,nTrees = 1,precBits = 1000,addStem = FALSE,nCores = 1)
    res <- cloneRate::internalLengths(test[[1]], alpha = 0.05)
    il= res$sumInternalLengths
    h = test[[2]]
    test_mcmc <- birthDeathMCMC(test[[1]],maxGrowthRate = 4,alpha = 0.05,verbose = TRUE,nChains = 4,nCores = 1,chainLength = 2000 )
    estimator_JS[i] = calc_estimator(n, h)
    res_new =  internalLengths(test[[1]], alpha = 0.05)
    ew_new <- ew_durrett_adapted(n, res_new[[2]])
    estimator[i] <-  ew_durrett/il # durrett
    estimator_old[i] <- res$estimate
    estimator_mcmc[i] <- test_mcmc$estimate
    estimator_new[i] <- ew_new/il
  }
  return(list(estimator_old, estimator, estimator_mcmc, estimator_new, estimator_JS))
}


res_all = repeat_simulation(10, 8, 1, 40)


print(mean((res_all[[5]]-1)^2))
print(mean((res_all[[4]]-1)^2))
print(mean((res_all[[3]]-1)^2))


repeat_simulation_T <- function(n_rep, T_list, r, n){
  K = length(T_list)
  results_repeat <- rep(0, K)
  estimator <- numeric(K)
  estimator_old <- c(K)
  estimator_mcmc <- c(K)
  estimator_new <- c(K)
  estimator_JS <- c(K)
  i = 1
  for(t in T_list){
    #print(t)
    temp <- repeat_simulation(n_rep, n, r, t)
    temp_df <- data.frame(
      list1 = temp[[1]],                    
      list2 = temp[[2]],                   
      list3 = temp[[3]],                    
      list4 = temp[[4]],                    
      list5 = temp[[5]]                   
    )
    # Ergebnisse in einer .txt-Datei speichern
    output_file = paste0("table_t_",t,"_r_",r,"_n_",n , ".txt")
    write.table(
      temp_df,
      file = output_file,
      sep = "\t",            # Tabulator als Trennzeichen
      row.names = FALSE,     # Keine Zeilennummern
      col.names = TRUE,      # Spaltennamen beibe
    )
    mse1 = mean((temp[[1]]-r)^2)
    mse2 = mean((temp[[2]]-r)^2)
    mse3 = mean((temp[[3]]-r)^2)
    mse4 = mean((temp[[4]]-r)^2)
    mse5 = mean((temp[[5]]-r)^2)
    estimator_old[i] = mse1
    estimator[i] = mse2
    estimator_mcmc[i] = mse3
    estimator_new[i] = mse4
    estimator_JS[i] = mse5
    i = i+1
  }
  return(list(estimator_old, estimator, estimator_mcmc, estimator_new, estimator_JS))
}
temp = repeat_simulation_T(100, c(40, 100), 1, 8)
