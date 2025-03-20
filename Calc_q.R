source("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\CPP_T_infinity.R")
# We need calc_estimator_adapted, simUltra from the main file

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
  return(list(alpha_opt_CH, q1, q3, ew_CH))
}

test_p = calc_prefactor(15, 1000, 0.5)
print(test_p)

eva_quantile <- function(n, r, t, numRep){
  ew_CH= rep(0, numRep)
  
  c_all = calc_prefactor(n, numRep, r)
  c_1 = c_all[[2]]
  c_3 = c_all[[3]]
  c_alpha = c_all[[1]]
  for(i in 1:numRep){
    a = runif(1, min = r, max = r+1)
    test <- simUltra_infty(a,b = a - r,n = n,nTrees = 1,precBits = 1000,addStem = FALSE,nCores = 1)
    
    h = test[[2]]
    # JS 2
    ew_CH[i] = calc_estimator_adapted(n, h)
  }
  v1 = ew_CH
  y1 <- rep(r, numRep) 
  alpha_opt_CH <- sum(v1 * y1) / sum(v1^2) # pre factor
  quantile_1 <- quantile(ew_CH*c_alpha, probs = 0.025) # 0.05 Quantile
  quantile_3 <- quantile(ew_CH*c_alpha, probs = 0.975) # 0.95 Quantile
  test_all_1 = c_1 *  ew_CH * c_alpha
  test_all_3 = c_3 *  ew_CH * c_alpha
  # Calculate the MSE between the empirical Quantiles and the estimated Quantiles
  mse_1 = mean((test_all_1-quantile_1)^2)
  mse_3 = mean((test_all_3-quantile_3)^2)
  
  return(list(quantile_1, quantile_3, test_all_1, test_all_3, mse_1, mse_3))
}
test_eva = eva_quantile(5, 2, 100, 1000)
test_eva_6 = eva_quantile(6, 2, 100, 1000)
test_eva_7 = eva_quantile(7, 2, 100, 1000)
test_eva_8 = eva_quantile(8, 2, 100, 1000)
test_eva_9 = eva_quantile(9, 2, 100, 1000)
test_eva_10 = eva_quantile(10, 2, 100, 1000)

print(test_eva[[5]]/test_eva[[1]])
print(test_eva[[6]]/test_eva[[2]])


rep_calc_prefactor <- function(n_range, num_Rep, r){
  res_alpha_CH = rep(0, length(n_range))
  res_q1 = rep(0, length(n_range))
  res_q3 = rep(0, length(n_range))

  i = 1
  for(n in n_range){
    res_alpha_temp = calc_prefactor(n, num_Rep, r)
    res_alpha_CH[i] = res_alpha_temp[[1]]
    res_q1[i] = res_alpha_temp[[2]]
    res_q3[i] = res_alpha_temp[[3]]
    i = i+1
  }
  return(list(res_alpha_CH, res_q3, res_q1))
}

res_alpha_all_40 = rep_calc_prefactor(c(5), 10, 2)

res_alpha_all = rep_calc_prefactor(c(5, 6, 7, 8, 9, 10), 1000, 2)
print(res_alpha_all)
# Originale Liste
my_list <- list(
  "JS2" = c(4.351666, 10.129650, 17.656991, 26.240726, 37.836560, 52.301474),
  "Durrett" = c(0.4218443, 0.5910854, 0.6585634, 0.6868385, 0.7190877, 0.7857804),
  "JS" = c(0.4577988, 0.6282367, 0.7242631, 0.7463504, 0.7992449, 0.8503860),
  "NL" = c(0.5913339, 0.7781026, 0.8741731, 0.8617310, 0.8586734, 0.9166453)
)

# Umwandlung in DataFrame
df <- as.data.frame(my_list)
# Ausgabe anzeigen

write.table(res_alpha_all[[1]], "prefactors_infty.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)



