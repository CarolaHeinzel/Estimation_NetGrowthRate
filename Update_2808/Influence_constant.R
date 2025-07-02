# Calculate the MSE for all n
res_MSE = rep(0, 5)
for(n in c(10)){
  t = 40
  r = 1
  if(n <= 10){
    filename <- paste0("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\Simulation_infinity\\CPPs_eva_40_n_",n,"_r_",r,".txt")
  }else if(n == 15){
    filename <- paste0("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\Simulation_infinity\\CPPs_final_n",n,"_r_05_T40.txt")
    # CPPs_eva_5_
    #filename <- paste0("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\Simulation_infinity\\CPPs_eva_6_40_n_15_r_0.5.txt")
    filename <- paste0("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\Simulation_infinity\\CPPs_final_n",n,"_r_",r,"_T40.txt")
  } else if(n == 20){
    #filename <- paste0("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\Simulation_infinity\\CPPs_eva_40_n_20_r_0.5.txt")
    #filename =  "C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\Simulation_infinity\\CPPs_eva_40_n_20_r_0.5.txt"
    filename <- paste0("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\Simulation_infinity\\CPPs_final_n",n,"_r_",r,"_T40.txt")
  }else{
    filename <- "C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\Simulation_infinity\\CPPs_eva_2_40_n_50_r_1.txt"
    
  } 
  daten <- read.table(filename, header = FALSE, sep = "", col.names = c("X1", "X2", "X3"))
  print(n)
  if(n ==10){
    c_n  = 9*8*0.70
    c_1 = 9*8*0.78
    c_2 = 9*8*15551/18144
  }
  eva_MSE = daten$X1

}

c_values <- seq(0.65, 0.9, by = 0.01) # MAD
c_values <- seq(0.5, 0.9, by = 0.01)

# MSE berechnen: mse(c * x, target)
mse <- sapply(c_values, function(c) mean((c * 9* 8 * eva_MSE - 1)^2))


mse <- sapply(c_values, function(c) mean(abs(c * 9* 8 * eva_MSE - 1)))

plot(c_values, mse, type = "l", col = "black", lwd = 2,
     xlab = "c", ylab = "MSE")


plot(c_values, mse, type = "l", col = "black", lwd = 2,
     xlab = "c", ylab = "Mean Absolute Deviation")

# Zwei senkrechte Linien bei c = 0.9 und c = 1.1 z. B.
abline(v = 0.78, col = "purple", lty = 2, lwd = 2)
abline(v = 0.72, col = "black", lty = 2, lwd = 2)
abline(v = 0.70, col = "green", lty = 2, lwd = 2)

best_c <- c_values[which.min(mse)]
best_mse <- min(mse)


#abline(v = 0.78, col = "red", lty = 2, lwd = 2)

#final_matrix <- t(cbind(first_elements, first_elements_6, first_elements_7, first_elements_8, first_elements_9, first_elements_10))
matrix_restuls1 <- t(res_MSE) #[,-c(4)]
left_column1 = c(5,6,7,8,9,10, 15, 20, 50)
final_matrix1 <- cbind(left_column1, matrix_restuls1[-c(1),])
colnames(final_matrix1) <- c("n","Phylofit", "Johnson", "New", "Unbiased for r", "Unbiased for 1/r" )

rownames(final_matrix1) <- NULL
additional_values <- mean_all # <- Anpassen!

# Neue Spalte in die Matrix einfügen
final_matrix1 <- cbind(final_matrix1, Additional = additional_values)
data = as.data.frame(final_matrix1)
data_long <- melt(data, id.vars = "n", variable.name = "Method", value.name = "Value")

library("ggplot2")
farben = c("#E41A1C", "#377EB8", "#4DAF4A" ,"#984EA3" ,"#FF7F00", "#44AA99")

ggplot(data_long, aes(x = factor(n), y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(#title = "",
    x = "n",
    y = "Mean Absolute Deviation") +
  theme_minimal() +
  scale_y_continuous(
    breaks = seq(0, 6, by = 0.1),
    sec.axis = sec_axis(~., breaks = NULL)) +
  scale_fill_manual(
    values = farben,
    labels = c(
      expression(hat(r)[Phylofit]),
      expression(hat(r)[Johnson]),
      expression(hat(r)[MSE]),
      expression(hat(r)[Bias]),
      expression(hat(r)[Inv]), 
      expression(hat(r)[MLE])
    )
  ) +
  theme(text = element_text(size = 14),
        legend.position = "top")
