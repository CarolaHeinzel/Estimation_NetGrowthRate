# Calculates the Estimator with the new idea, i.e. \hat r_{new}
source("Simulation_self.R")
is_valid_permutation <- function(h, k, n) {
  if (k > 1) {
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
  
    if (h[n - k] > max(h[(n - k+1):(n-1)])) {
      return(FALSE)
    }
  }
  return(TRUE)
}

ew_durrett_adapted <- function(n, tree_data){
  sum_result <- 0
  temp = calculate_leaf_info(tree_data, n)
  h_values <-  seq(1, (n-1))
  permutations <- permutations(n = length(h_values), r = length(h_values), v = h_values)
  for (k in 1:(n - 2)) {
    if(temp$check_leaves$Exists[k] == TRUE){
      print(k)
      valid_permutations <- permutations[apply(permutations, 1, function(perm) is_valid_permutation(perm, k+1, n)), ]
      ratio <- nrow(valid_permutations) / factorial(n-1)
      
      value <-  n / (k*(k+1)) *1/(1-ratio)
      sum_result <- sum_result + value
    }
  }
  return(sum_result)
}
