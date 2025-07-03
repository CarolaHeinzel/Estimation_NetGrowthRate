# Calculate the influence of the constant on the MSE
res_MSE = rep(0, 5)
for(n in c(10)){
  t = 40
  r = 1
  # read the data
  filename <- paste0("CPPs_40_n_",n,"_r_",r,".txt")
 
  daten <- read.table(filename, header = FALSE, sep = "", col.names = c("X1", "X2", "X3"))
  eva_MSE = daten$X1 # this step depends on the data format
}

c_values <- seq(0.5, 0.9, by = 0.01)

# calculate mse(c * x, r)
mse <- sapply(c_values, function(c) mean((c * 9* 8 * eva_MSE - 1)^2))

# Plot the results
plot(c_values, mse, type = "l", col = "black", lwd = 2,
     xlab = "c", ylab = "MSE")
abline(v = 0.78, col = "purple", lty = 2, lwd = 2)
abline(v = 0.72, col = "black", lty = 2, lwd = 2)
abline(v = 0.70, col = "green", lty = 2, lwd = 2)
