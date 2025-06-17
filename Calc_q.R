simUltra_infty <- function(a, b, n, nTrees = 1,
                           precBits = 1000, addStem = FALSE, nCores = 1) {
  # Store runtime for each tree
  ptm <- proc.time()
  # Convert params to high precision mpfr
  a_mpfr <- mpfr(a, precBits)
  b_mpfr <- mpfr(b, precBits)
  net_mpfr <- a_mpfr - b_mpfr
  n_mpfr <- mpfr(n, precBits)
  one <- mpfr(1, precBits)
  
  # Draw Y = y from the inverse CDF
  uniform_rv <- mpfr(stats::runif(1, min = 0, max = 1), precBits)
  y_mpfr <- uniform_rv**(1/n)/(1-uniform_rv**(1/n))
  # Generate the coalescence times
  coal_times_mpfr <- sapply(rep(y_mpfr, n - 1), inv_cdf_coal_times_inf,
                            net = net_mpfr,
                            a = a_mpfr, precBits = precBits
  )
  # Convert back to normal numeric (no longer need high precision)
  coal_times <- suppressWarnings(sapply(coal_times_mpfr, Rmpfr::asNumeric))
  # print(coal_times)
  # Convert coal times into tree by randomly merging lineages
  tree <- coal_to_tree(coal_times)
  
  # Add stem starting the tree from zero, rooting the tree appropriately
  if (addStem) {
    tree$edge[tree$edge > n] <- tree$edge[tree$edge > n] + 1
    tree$edge <- rbind(c(n + 1, n + 2), tree$edge)
    tree$edge.length <- c(max(coal_times), tree$edge.length)
    tree$Nnode <- tree$Nnode + 1
  }
  # Add metadata for making the tree
  runtime <- proc.time()[["elapsed"]] - ptm[["elapsed"]]
  tree$metadata <- data.frame(
    "r" = a - b, "a" = a, "b" = b,
    "n" = n, "runtime_seconds" = runtime, "addStem" = addStem
  )
  # Return the tree created from the coalescence times drawn from Lambert distribution
  return(list(tree, coal_times))
}
compute_list_properties <- function(liste, n) {
  normal_list = liste
  max_value <- max(normal_list)
  mean_value <- mean(normal_list)
  
  # ∑_i ∑_j (H_i - H_j)^+
  double_sum <- sum(outer(normal_list, normal_list, function(x, y) pmax(0, x - y)))
  
  return(list(Maximum = max_value, Mittelwert = mean_value, d_sum = double_sum))
}


calc_estimator_adapted <- function(n, h){
  
  values <- compute_list_properties(h, n)
  est <-  values$d_sum
  return(1/est)
  
}

calc_prefactor <- function(n, numRep, r){
  ew_CH <- rep(0, numRep)
  ew_CH_f <- rep(0, numRep)
  
  for(i in 1:numRep){
    a = runif(1, min = r, max = r+1)
    test <- simUltra_infty(a,b = a - r,n = n,nTrees = 1,precBits = 1000,addStem = FALSE,nCores = 1)
    
    h = test[[2]]
    # JS 2
    ew_CH[i] = calc_estimator_adapted(n, h)
    ew_CH_f[i] = calc_estimator_adapted(n, test_f[[2]])
    
  }
  v1 = ew_CH
  y1 <- rep(r, numRep) 
  alpha_opt_CH <- sum(v1 * y1) / sum(v1^2) # pre factor
  quantile_1 <- quantile(ew_CH, probs = 0.025) # 0.05 Quantile
  quantile_3 <- quantile(ew_CH, probs = 0.975) # 0.95 Quantile
  q1 = quantile_1/mean(ew_CH)
  q3 = quantile_3/mean(ew_CH)
  # q1 and q3 are the constants for the 0.025 and the 0.975 quantile, respectively.
  return(list(alpha_opt_CH, q1, q3, ew_CH))
}

# Example application
test_p = calc_prefactor(5, 10, 0.5)
print(test_p)
