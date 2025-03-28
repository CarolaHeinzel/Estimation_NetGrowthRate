# Calculates the Estimator with the new idea

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
    #print(h[n-k])
    #print(h[(n - k+1):(n-1)])
    if (h[n - k] > max(h[(n - k+1):(n-1)])) {
      return(FALSE)
    }
  }
  return(TRUE)
}

k <- 4
h_values = seq(1, 5)
n <- length(h_values) +1
permutations <- permutations(n = length(h_values), r = length(h_values), v = h_values)
valid_permutations <- permutations[apply(permutations, 1, function(perm) is_valid_permutation(perm, k, n)), ]
ratio <- nrow(valid_permutations) / factorial(n-1)
cat(sprintf("Ratio of Valid Permutations: %.4f\n", ratio))


calculate_leaf_info <- function(tree_data, n) {
  # Calculate the leaves
  leaves <- setdiff(tree_data$Node, tree_data$Parent)
  
  get_leaves <- function(node, tree) {
    children <- tree$Node[tree$Parent == node]
    if (length(children) == 0) {
      return(node)
    } else {
      return(unlist(lapply(children, get_leaves, tree = tree)))
    }
  }
  
  results <- data.frame(Edge = integer(), Leaves = integer(), Total_Length = numeric())
  # Calculate the lenght of branches that support k = 2,.., n-1 leaves
  for (i in 1:nrow(tree_data)) {
    node <- tree_data$Node[i]
    edge_length <- tree_data$Edge_length[i]
    
    descendant_leaves <- get_leaves(node, tree_data)
    
    num_leaves <- length(descendant_leaves)
    total_length <- edge_length
    
    results <- rbind(
      results,
      data.frame(Edge = node, Leaves = num_leaves, Total_Length = total_length)
    )
  }
  
  summed_results <- aggregate(Total_Length ~ Leaves, data = results, FUN = sum)
  
  check_leaves <- data.frame(Leaves = 1:(n - 1), Exists = FALSE, Total_Length_Positive = FALSE)
  
  for (k in check_leaves$Leaves) {
    # Check wheter k occurs
    exists <- any(results$Leaves == k)
    positive_length <- any(results$Leaves == k & results$Total_Length > 0)
    
    
    check_leaves$Exists[check_leaves$Leaves == k] <- exists
    check_leaves$Total_Length_Positive[check_leaves$Leaves == k] <- positive_length
  }
  return(list(results = results, summed_results = summed_results, check_leaves = check_leaves))
}

ew_durrett_adapted <- function(n, tree_data){
  sum_result <- 0
  temp = calculate_leaf_info(tree_data, n)
  h_values <-  seq(1, (n-1))
  permutations <- permutations(n = length(h_values), r = length(h_values), v = h_values)
  value1 = 0
  number_all = 0
  number_ex = 0
  res_th = 0
  for (k in 1:(n - 2)) {
    if(k == 1){
      number_all = number_all + 1
      valid_permutations <- permutations[apply(permutations, 1, function(perm) is_valid_permutation(perm, k+1, n)), ]
      ratio <- nrow(valid_permutations) / factorial(n-1)
      res_th = 1/(1-ratio)
    }
    if(temp$check_leaves$Exists[k] == TRUE){
    #if(2>1){
      if(k ==1){
        number_ex = number_ex + 1
      }
      
      valid_permutations <- permutations[apply(permutations, 1, function(perm) is_valid_permutation(perm, k+1, n)), ]
      ratio <- nrow(valid_permutations) / factorial(n-1)
      if(ratio != 0){
        value <-  n / (k*(k+1)) #*1/(ratio)
      }else{
        value <-  n / (k*(k+1)) 
        
      }
      #print(1/(1-ratio))
      sum_result <- sum_result + value
    }
    value1 <- value1 + n / (k*(k+1))
  }
  return(list(sum_result, value1, number_all, number_ex, res_th))
}

rep_ew_durrett_adapted <- function(numRep, n){
  number_ex = rep(0, numRep)
  number_all = rep(0, numRep)
  for(i in 1:numRep){
    test =simUltra(3,1,100,n) # lambda = 2
    res_new = internalLengths(test[[1]])

    res_temp = ew_durrett_adapted(n, res_new[[2]])
    number_all[i] = res_temp[[3]]
    number_ex[i] = res_temp[[4]]
    res_th_temp = res_temp[[5]]
    if(res_th_temp!= 0){
      res_th = res_th_temp
    }
  }
  return(list(res_th, sum(number_ex), sum(number_all)))
}


res_bootstrap = rep_ew_durrett_adapted(100, 5)
print(res_bootstrap)
n = 6
test =simUltra(3,1,100,n)
res_new = internalLengths(test[[1]])
input = calculate_leaf_info(res_new[[2]],n)
print(res_new[[2]])
print(res_new[[1]]$sumInternalLengths)
res = ew_durrett_adapted(n, res_new[[2]])
print(res)
