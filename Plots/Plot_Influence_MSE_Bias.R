t = 40
n = 10
r = 0.5

filename <- paste0("CPPs_2_",t,"_n_", n, "_r_", r, ".txt")
daten <- read.csv(filename, header = FALSE)

# Define values for constant
c_values <- seq(40, 60.0, by=0.1)  
mse_values <- numeric(length(c_values))  
bias_values <- numeric(length(c_values))  

# MSE 
for (i in seq_along(c_values)) {
  c <- c_values[i]
  mse_values[i] <- mean((c * daten$V1 - r)^2)
}
# Bias
for (i in seq_along(c_values)) {
  c <- c_values[i]
  bias_values[i] <- mean((c * daten$V1 - r))
}

# Plot the Bias
par(mar = c(5, 6, 4, 2))  
plot(c_values, bias_values, type="l", col="blue", lwd=2,
     xlab="c(10,40, 0.5)", ylab="Bias")
abline(v = 50.30939, col="black", lty=2)
grid()

# Plot the MSE
par(mar = c(5, 6, 4, 2))  
plot(c_values, mse_values, type="l", col="blue", lwd=2,
     xlab="c(10,40, 0.5)", ylab=expression(MSE[c(10, 40, 0.5)]))
abline(v = 50.30939, col="black", lty=2)
grid()
