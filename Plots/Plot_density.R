filename <- paste0("CPPs_infinity_2_5_r_1.txt") # Load data

daten <- read.csv(filename, header = FALSE)

density_values <- density(daten)

# Plot of density function
plot(density_values, 
     xlab = expression(tilde(c)[MSE](5)), 
     main = "",
     ylab = "Density", 
     col = "blue", 
     lwd = 2,
     xlim = c(0, 1))
