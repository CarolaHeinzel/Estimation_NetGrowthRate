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
  if(n==5){
    c_n <- 12*115/144
    c_2 = 0.9505112
  }else if(n == 6){
    c_n <- 5*4*163/200
    c_2 = 1.0457389
  }else if(n == 7){
    c_n <- 6*5*497/600
    c_2 =1.0557716 
  }else if(n == 8){
    c_n <- 7*6*617/735
    c_2 = 1.0489458 
  }else if(n == 9){
    c_n <- 8*7*13311/15680
    c_2 = 1.0555155 
  }else if(n == 10){
    c_n <- 9*8*15551/18144
    c_2 = 1.064129
  }else{
    c_n = 1
  }
  values <- compute_list_properties(h, n)
  est <-  values$d_sum
  return(list(c_n/est, 1/est, c_n*c_2/est))
  
}


calc_prefactor <- function(n, numRep, r){
  ew_CH <- rep(0, numRep)
  ew_CH_unbiased_r = rep(0, numRep)
  temp = rep(0, numRep)
  for(i in 1:numRep){
    a = runif(1, min = r, max = r+1)
    test <- simUltra_infty(a,b = a - r,n = n,nTrees = 1,precBits = 1000,addStem = FALSE,nCores = 1)
    
    h = test[[2]]
    # JS 2
    temp[i] = calc_estimator_adapted(n, h)[[1]] # with c_n as prefactor
    ew_CH[i] = calc_estimator_adapted(n, h)[[2]]
    ew_CH_unbiased_r[i] = calc_estimator_adapted(n, h)[[3]]
  }
  c_unbiased = r/mean(temp)
  #print(temp)
  #print(mean(temp))
  v1 = ew_CH
  y1 <- rep(r, numRep) 
  alpha_opt_CH <- sum(v1 * y1) / sum(v1^2) # pre factor
  quantile_1 <- quantile(ew_CH, probs = 0.025) # 0.05 Quantile
  quantile_3 <- quantile(ew_CH, probs = 0.975) # 0.95 Quantile
  q1 = quantile_1/mean(ew_CH)
  q3 = quantile_3/mean(ew_CH)
  return(list(alpha_opt_CH, q1, q3, c_unbiased, ew_CH_unbiased_r, temp))
}


rep_calc_prefactor <- function(n_range, num_Rep, r){
  res_alpha_CH = rep(0, length(n_range))
  res_q1 = rep(0, length(n_range))
  res_q3 = rep(0, length(n_range))
  c_unbiased = rep(0, length(n_range))
  mse = rep(0, length(n_range))
  mse_c_n = rep(0, length(n_range))
  i = 1
  for(n in n_range){
    res_alpha_temp = calc_prefactor(n, num_Rep, r)
    res_alpha_CH[i] = res_alpha_temp[[1]]
    res_q1[i] = res_alpha_temp[[2]]
    res_q3[i] = res_alpha_temp[[3]]
    c_unbiased[i] = res_alpha_temp[[4]]
    #print(res_alpha_temp[[5]])
    mse[i] = mean((res_alpha_temp[[5]] - r)^2)
    mse_c_n[i] = mean((res_alpha_temp[[6]] - r)^2)
    
    
    i = i+1
  }
  return(list(res_alpha_CH, res_q3, res_q1, c_unbiased, mse, mse_c_n))
}
res_alpha_all_unbiased_mse = rep_calc_prefactor(c(5, 6, 7, 8, 9, 10), 1000, 0.5)
res_alpha_all_unbiased_mse_1 = rep_calc_prefactor(c(5, 6, 7, 8, 9, 10), 100, 1)

res_alpha_all_unbiased = rep_calc_prefactor(c(6), 1000, 1)


# Daten definieren
values_5 <- c(0.72864421, 0.31987382, 0.38154079, 0.10812473, 0.11494817, 0.09430849)
values_6 <- c(0.81951030, 0.28661739, 0.33727112, 0.10345704, 0.10922361, 0.08728605)

# Kombiniere die Daten zu einer Matrix, sodass jede Spalte eine Gruppe ist
data_matrix <- rbind(values_5, values_6)

# Balkendiagramm erstellen
barplot(data_matrix, beside = TRUE, 
        col = c("orange", "yellow"), 
        legend.text = c("unbiased for r", "unbiased for 1/r"),
        main = "r = 1, T = 40",
        xlab = "n",
        ylab = "MSE",
        names.arg = 5:10)


# Daten definieren
values_5 <- c(0.15026509, 0.06778794, 0.06063515, 0.03657883, 0.03237710, 0.02941234)
values_6 <- c(0.16813078, 0.06181772, 0.05425067, 0.03369801, 0.02957533, 0.02673182)

# Kombiniere die Daten zu einer Matrix
data_matrix <- rbind(values_5, values_6)

barplot(data_matrix, beside = TRUE, 
        col = c("orange", "yellow"), 
        legend.text = c("unbiased for r", "unbiased for 1/r"),
        main = "r = 0.5, T = 40",
        xlab = "n",
        ylab = "MSE",
        names.arg = 5:10)


res_alpha_all = rep_calc_prefactor(c(5, 6, 7, 8, 9, 10), 1000, 2)
res_alpha_all_1 = rep_calc_prefactor(c(5, 6, 7, 8, 9, 10), 1000, 1)
res_alpha_all_2 = rep_calc_prefactor(c(5, 6, 7, 8, 9, 10), 1000, 0.5)
res_alpha_all_3 = rep_calc_prefactor(c(5, 6, 7, 8, 9, 10), 1000, 2)
res_alpha_all_4 = rep_calc_prefactor(c(5, 6, 7, 8, 9, 10), 1000, 2)
res_alpha_all_5 = rep_calc_prefactor(c(5, 6, 7, 8, 9, 10), 1000, 2)
res_alpha_all_6 = rep_calc_prefactor(c(5, 6, 7, 8, 9, 10), 1000, 2)
res_alpha_all_7 = rep_calc_prefactor(c(5, 6, 7, 8, 9, 10), 1000, 2)
res_alpha_all_8 = rep_calc_prefactor(c(5, 6, 7, 8, 9, 10), 1000, 2)
res_alpha_all_9 = rep_calc_prefactor(c(5, 6, 7, 8, 9, 10), 1000, 2)


# Liste aller Vektoren (ersetze die Anzahl der Listen entsprechend deiner Daten)
vector_list <- list(res_alpha_all_1[[1]], res_alpha_all[[1]], res_alpha_all_2[[1]], res_alpha_all_3[[1]], res_alpha_all_4[[1]], res_alpha_all_5[[1]], res_alpha_all_6[[1]], res_alpha_all_7[[1]], res_alpha_all_8[[1]], res_alpha_all_9[[1]])  # Falls es mehr gibt, erweitern

# Berechnung des komponentenweisen Mittelwerts
mean_vector <- Reduce("+", vector_list) / length(vector_list)

# Ausgabe des Ergebnisvektors
print(mean_vector)


vector_matrix <- do.call(rbind, vector_list)
par(mfrow=c(1, 1))  # Passe an, um alle Boxplots in einer Zeile anzuzeigen
boxplot(vector_matrix[,1], 
                 xlab="n = 5", 
                 ylab="c(5)", 
                 col="lightblue", 
                 border="darkblue")
points(rep(1, nrow(vector_matrix)), vector_matrix[, 1], 
       col="red", pch=16) 
