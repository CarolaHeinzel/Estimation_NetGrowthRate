t = 40
n = 10
r = 0.5

filename <- paste0("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\Simulation_infinity\\CPPs_2_",t,"_n_", n, "_r_", r, ".txt")


# Datei einlesen (ersetze 'datei.csv' durch den tatsächlichen Dateinamen)
daten <- read.csv(filename, header = FALSE)



# Wertebereich für c definieren
c_values <- seq(40, 60.0, by=0.1)  # c von 0.1 bis 2 in Schritten von 0.1
mse_values <- numeric(length(c_values))  # Leeres Vektor für MSE-Werte

# Berechnung des MSE für verschiedene c-Werte
for (i in seq_along(c_values)) {
  c <- c_values[i]
  mse_values[i] <- mean((c * daten$V1 - r)^2)
}

for (i in seq_along(c_values)) {
  c <- c_values[i]
  mse_values[i] <- mean((c * daten$V1 - r))
}

par(mar = c(5, 6, 4, 2))  # Erhöht linken Rand für y-Achsenbeschriftung
# MSE in Abhängigkeit von c plotten
plot(c_values, mse_values, type="l", col="blue", lwd=2,
     xlab="c(10,40, 0.5)", ylab="Bias")
abline(v = 50.30939, col="black", lty=2)

grid()

par(mar = c(5, 6, 4, 2))  # Erhöht linken Rand für y-Achsenbeschriftung
# MSE in Abhängigkeit von c plotten
plot(c_values, mse_values, type="l", col="blue", lwd=2,
     xlab="c(10,40, 0.5)", ylab=expression(MSE[c(10, 40, 0.5)]))
abline(v = 50.30939, col="black", lty=2)

grid()
