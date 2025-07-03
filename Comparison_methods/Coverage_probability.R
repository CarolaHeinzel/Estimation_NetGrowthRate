q1 = rep(0, 8)
q3 = rep(0,8)
i =1
# Calculate the constants
for(n in c(5,6,7,8,9,10,15)){
  filename = paste0("CPPs_infinity_",n,"_r_1.txt")
  daten <- read.table(filename, header = FALSE, sep = "", col.names = c("X1"))
  q1[i] = quantile(daten$X1, probs = c(0.025))*(n-1)*(n-2)
  q3[i] = quantile(daten$X1, probs = c(0.975))*(n-1)*(n-2)
  i = i+1
}

i = 1
# Empirical coverage probability
rel_f_all = rep(0,8)
# Evaluate the quality of the constants
for(n in c(5,6,7,8,9,10,15,20)){
  t = 40
  r = 1
  filename <- paste0("CPPs_40_n_",n,"_r_",r,".txt")
  daten <- read.table(filename, header = FALSE, sep = "", col.names = c("X1", "X2", "X3"))
  eva_MSE = daten$X1 # Estimator with constant 1
  rel_f = 0

  
  for (hat_r in eva_MSE){
    print(hat_r)
    if(r <= (hat_r)/q1[i]){
      if(r >= (hat_r)/q3[i]){
        rel_f = rel_f + 1 
      }
    }
  }
  rel_f_all[i] = rel_f
  i = i+1  
}

