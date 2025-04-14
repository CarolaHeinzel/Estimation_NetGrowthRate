# Calculate the optimal constants

n = 20
filename_1 <- paste0("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\Simulation_infinity\\CPPs_infinity_1_",n,"_r_1.txt")

# Datei einlesen (ersetze 'datei.csv' durch den tatsächlichen Dateinamen)
daten <- read.csv(filename_1, header = FALSE)
mean(daten$V1)
y1 <- rep(1, length(daten$V1)) 
v1 = daten$V1
alpha_opt_CH <- sum(v1 * y1) / sum(v1^2) 
print(alpha_opt_CH)
print(1/mean(daten$V1))


# Plot the densities of
# the three new estimators, Phylofit, Old one
t = 40
n = 20
r = 0.5

filename <- paste0("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\Simulation_infinity\\CPPs_eva_",t,"_n_", n, "_r_", r, ".txt")


daten <- read.table(filename, header = FALSE, sep = "", col.names = c("X1", "X2", "X3"))
colnames(daten) <- c("X1", "X2", "X3")
if(n ==10){
  c_n  = 50.30939
  c_1 = 56.0427
  c_2 = 9*8*15551/18144
}else if(n==5){
  c_n = 4.185482
  c_1 =  7.074815
  c_2 = 4 *3 *115/144
}else if(n==6){
  c_n = 9.917149
  c_1 =  13.23726
  c_2 = 5*4*163/200
}else if(n==7){
  c_n = 17.22395
  c_1 =  21.16463
  c_2 =6*5*497/600
}else if(n==8){
  c_n = 26.41381
  c_1 =  30.91002
  c_2 = 7*6*617/735	
}else if(n==9){
  c_n = 37.47819
  c_1 =  42.57689
  c_2 = 8*7*13311/15680 
}else if(n == 15){
  c_n =  143.479
  c_1 = 152.3723
  c_2 = 14*13*3873307/4372368
}else if(n == 20){
  c_n =  284.5928
  c_1 = 296.7941
  c_2 = 19*18*1199057081/1326917592
}
# N = 10
eva_MCMC = daten$X3
eva_old = n*daten$X2
eva_MSE = c_n*daten$X1
eva_unbiased = c_1*daten$X1
eva_biased =  c_2*daten$X1

d1 <- density(eva_MCMC)
d2 <- density(eva_old)
d3 <- density(eva_MSE)
d4 <- density(eva_unbiased)
d5 <- density(eva_biased)
#my_lists <- list(list1, list2, list3, list4, list5)
densities <- list(d1, d2, d3, d4, d5)
y_max <- max(sapply(densities, function(d) max(d$y)))
x_max <- max(sapply(densities, function(d) max(d$x)))

# Plot der ersten Dichte (legt die Achsen etc. fest)
plot(d1, col = "red", xlim = c(0, x_max), ylim = c(0, y_max), lwd = 2, main = "Densities of the Estimators for n = 20", xlab = "Value", ylab = "Density")

# Weitere Dichten hinzufügen
lines(d2, col = "blue", lwd = 2)
lines(d3, col = "green", lwd = 2)
lines(d4, col = "purple", lwd = 2)
lines(d5, col = "orange", lwd = 2)
# n= 5
eva_MCMC_5 = daten$X3
eva_old_5 = n*daten$X2
eva_MSE_5 = c_n*daten$X1
eva_unbiased_5 = c_1*daten$X1
eva_biased_5 =  c_2*daten$X1

d1_5 <- density(eva_MCMC_5)
d2_5 <- density(eva_old_5)
d3_5 <- density(eva_MSE_5)
d4_5 <- density(eva_unbiased_5)
d5_5 <- density(eva_biased_5)

lines(d1_5, col = "blue", lwd = 2, lty = 2)
lines(d2_5, col = "blue", lwd = 2, lty = 2)
lines(d3_5, col = "green", lwd = 2, lty = 2)
lines(d4_5, col = "purple", lwd = 2, lty = 2)
lines(d5_5, col = "orange", lwd = 2, lty = 2)

# Legende hinzufügen
legend("topright",
       legend = c("Phylofit", "Johnson", "Minimal MSE", "Unbiased for r", "Unbiased for 1/r"),
       col = c("red", "blue", "green", "purple", "orange"),
       lwd = 2)
abline(v = r, col = "black", lty = 2, lwd = 2)
abline(v = mean(eva_MCMC), col = "red", lty = 2, lwd = 2)
abline(v = mean(eva_old), col = "blue", lty = 2, lwd = 2)
abline(v = mean(eva_new), col = "green", lty = 2, lwd = 2)
abline(v = mean(eva_unbiased), col = "purple", lty = 2, lwd = 2)
abline(v = mean(eva_biased), col = "orange", lty = 2, lwd = 2)

# Calculate the MSE for all n
res_MSE = rep(0, 5)
for(n in c(5,6,7,8,9,10,15,20)){
  t = 40
  r = 1
  
  filename <- paste0("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\Simulation_infinity\\CPPs_eva_",t,"_n_", n, "_r_", r, ".txt")
  daten <- read.table(filename, header = FALSE, sep = "", col.names = c("X1", "X2", "X3"))
  colnames(daten) <- c("X1", "X2", "X3")
  if(n ==10){
    c_n  = 50.30939
    c_1 = 56.0427
    c_2 = 9*8*15551/18144
  }else if(n==5){
    c_n = 4.185482
    c_1 =  7.074815
    c_2 = 4 *3 *115/144
  }else if(n==6){
    c_n = 9.917149
    c_1 =  13.23726
    c_2 = 5*4*163/200
  }else if(n==7){
    c_n = 17.22395
    c_1 =  21.16463
    c_2 =6*5*497/600
  }else if(n==8){
    c_n = 26.41381
    c_1 =  30.91002
    c_2 = 7*6*617/735	
  }else if(n==9){
    c_n = 37.47819
    c_1 =  42.57689
    c_2 = 8*7*13311/15680 
  }else if(n == 15){
    c_n =  143.479
    c_1 = 152.3723
    c_2 = 14*13*3873307/4372368
  }else if(n == 20){
    c_n =  284.5928
    c_1 = 296.7941
    c_2 = 19*18*1199057081/1326917592
  }
  eva_MCMC = daten$X3
  eva_old = n*daten$X2
  eva_MSE = c_n*daten$X1
  eva_unbiased = c_1*daten$X1
  eva_biased =  c_2*daten$X1
  
  mse_MCMC = mean((eva_MCMC - r)^2)
  mse_old = mean((eva_old - r)^2)
  mse_MSE = mean((eva_MSE - r)^2)
  print(mse_MSE)
  mse_unbiased = mean((eva_unbiased - r)^2)
  mse_biased = mean((eva_biased - r)^2)
  res_temp = c(mse_MCMC, mse_old, mse_MSE, mse_unbiased, mse_biased)
  res_MSE = cbind(res_MSE, res_temp)
}
#final_matrix <- t(cbind(first_elements, first_elements_6, first_elements_7, first_elements_8, first_elements_9, first_elements_10))
matrix_restuls1 <- t(res_MSE) #[,-c(4)]

left_column1 = c(5,6,7,8,9,10, 15, 20)
final_matrix1 <- cbind(left_column1, matrix_restuls1[-c(1),])
colnames(final_matrix1) <- c("n","Phylofit", "Johnson", "New", "Unbiased for r", "Unbiased for 1/r" )

rownames(final_matrix1) <- NULL

data = as.data.frame(final_matrix1)

data_long <- melt(data, id.vars = "n", variable.name = "Method", value.name = "Value")

# Create Plot
ggplot(data_long, aes(x = factor(n), y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Comparison of the Estimators for r = 1",
       x = "n",
       y = "MSE") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 14))

r = 1
n = 100
t = 40
filename <- paste0("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\Simulation_infinity\\CPPs_eva_",t,"_n_", n, "_r_", r, ".txt")

read_cn <- read.table("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\Simulation_infinity\\CPPs_infinity_time_1_100_r_1.txt", header = FALSE)
mean(read_cn$V1)
y1 <- rep(1, length(read_cn$V1)) 
v1 = read_cn$V1
alpha_opt_CH <- sum(v1 * y1) / sum(v1^2) 
print(alpha_opt_CH)
daten <- read.table(filename, header = FALSE, sep = "", col.names = c("X1", "X2", "X3"))
colnames(daten) <- c("X1", "X2", "X3")

eva_MCMC = daten$X3
eva_old = n*daten$X2
eva_MSE = alpha_opt_CH*daten$X1

mse_MCMC = mean((eva_MCMC - r)^2)
mse_old = mean((eva_old - r)^2)
mse_MSE = mean((eva_MSE - r)^2)
print(mse_MCMC)
print(mse_old)
print(mse_MSE)
