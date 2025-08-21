# Calculation of the new estimator with constant c(n), 
# i.e. the result is \hat r_{MSE}, \hat r_{Bias} or \hat r_{Inv}, depending
# on the choice of the constant c

# Calculates ∑_i ∑_j (H_i - H_j)^+
# h is a list that contains the coalescence times
compute_list_properties <- function(h) {
  double_sum <- sum(outer(h, h, function(x, y) pmax(0, x - y)))
  return(double_sum)
}

# Calculates the estimator with the constant c (which has to be determined)
# n is the samples size
calc_estimator <- function(h, c){
  n = length(h) + 1
  values <- compute_list_properties(h)
  est <-  1/(n-1)*values*c
  return(est)
}

# Exaple application
h = c(1,4,5,2,3)
c = 1
r = calc_estimator(h, c)
