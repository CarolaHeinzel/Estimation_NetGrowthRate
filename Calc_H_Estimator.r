library(Rmpfr)


# Funktion zur Berechnung von t
inv_cdf_coal_times <- function(y, net, a, alpha, precBits) {
  one <- mpfr(1, precBits)
  rv <- mpfr(stats::runif(1, min = 0, max = 1), precBits)
  print(rv)
  phi <- (alpha * rv) / (a * (one - alpha * (one - y)))
  return((-one / net) * log((one - y * a * phi) / (one + (net - y * a) * phi)))
}


calc_H <- function(a, b, T_coal, n, precBits){
  
  cloneAge_mpfr <- mpfr(T_coal, precBits)
  a_mpfr <- mpfr(a, precBits)
  b_mpfr <- mpfr(b, precBits)
  net_mpfr <- a_mpfr - b_mpfr
  n_mpfr <- mpfr(n, precBits)
  one <- mpfr(1, precBits)
  
  # Calculate a in Lambert (2018)
  alpha_mpfr <- (a_mpfr * (exp(net_mpfr * cloneAge_mpfr) - one)) / (a_mpfr * exp(net_mpfr * cloneAge_mpfr) - b_mpfr)
  if (alpha_mpfr == one) {
    stop("alpha value is equal to 1 due to insufficient machine precision. This
          will lead to NaN coalescence times. Increase param 'precBits'.")
  }
  # Calculate Y = y 
  uniform_rv <- mpfr(stats::runif(1, min = 0, max = 1), precBits)
  y_mpfr <- ((one - alpha_mpfr) * (uniform_rv**(one / n_mpfr))) / (one - alpha_mpfr * uniform_rv**(one / n_mpfr))

  # Calculate the y n-1 times
  h_list <- rep(0, n-1)
  lambda = a
  r = a-b
  for(i in 1:n-1){
    h_list[i] = inv_cdf_coal_times(y_mpfr, net_mpfr , a_mpfr, alpha_mpfr, precBits) 
  }
  coal_times <- suppressWarnings(sapply(h_list, Rmpfr::asNumeric))
  
  return(coal_times)
}


h = calc_H(2,1,40, 5, 100)
print(h)

compute_list_properties <- function(liste, n) {
  # Maximum
  normal_list = liste
  max_value <- max(normal_list)
  
  # Mittelwert
  mean_value <- mean(normal_list)
  
  # Doppelsumme: ∑_i ∑_j (l_i - l_j)^+
  double_sum <- sum(outer(normal_list, normal_list, function(x, y) pmax(0, x - y)))
  
  return(list(Maximum = max_value, Mittelwert = mean_value, Doppelsumme = double_sum))
}




calc_estimator <- function(n){
  h = calc_H(2,1,40, 5, 100)
  print(h)
  values <- compute_list_properties(h, n)
  print(values)
  est <- values$Maximum - values$Mittelwert + 1/(n-1)*values$Doppelsumme
  ev_th <- ew_durrett(n)
  
  return(est/ev_th)
  
}
est = calc_estimator(5)
print(est)
