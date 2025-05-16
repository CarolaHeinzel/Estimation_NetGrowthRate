# Calculate the MSE for all n
res_MSE = rep(0, 5)
for(n in c(5,6,7,8,9,10,15,20,50)){
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
  }else if(n==5){
    c_n = 4 * 3 * 0.355
    c_1 =  4*3*0.59
    c_2 = 4 *3 *115/144
  }else if(n==6){
    c_n = 5*4*0.49
    c_1 =  5*4*0.66
    c_2 = 5*4*163/200
  }else if(n==7){
    c_n = 6*5*0.58
    c_1 =  6*5*0.71
    c_2 =6*5*497/600
  }else if(n==8){
    c_n = 7*6*0.63
    c_1 =  7*6*0.74
    c_2 = 7*6*617/735	
  }else if(n==9){
    c_n = 8*7*0.67
    c_1 =  8*7*0.76
    c_2 = 8*7*13311/15680 
  }else if(n == 15){
    c_n =  14*13*0.79
    c_1 = 14*13*0.84
    c_2 = 14*13*3873307/4372368
  }else if(n == 20){
    c_n =  19*18*0.83
    c_1 = 19*18*0.83
    c_2 = 19*18*1199057081/1326917592
  }else if(n == 50){
    c_n =  49*48*0.92
    c_1 = 49*48*0.93
    c_2 = 49*48*137971924020914703586969/145779053479731685069056
  }
  eva_MCMC = daten$X3
  eva_old = n*daten$X2
  eva_MSE = c_n*daten$X1
  eva_unbiased = c_1*daten$X1
  eva_biased =  c_2*daten$X1
  mse_MCMC = mean(abs((eva_MCMC - r)))
  mse_old = mean(abs((eva_old - r)))
  mse_MSE = mean(abs((eva_MSE - r)))
  print(mse_MSE)
  mse_unbiased = mean(abs((eva_unbiased - r)))
  mse_biased = mean(abs((eva_biased - r)))
  res_temp = c(mse_MCMC, mse_old, mse_MSE, mse_unbiased, mse_biased)
  print(res_temp)
  res_MSE = cbind(res_MSE, res_temp)
}
#final_matrix <- t(cbind(first_elements, first_elements_6, first_elements_7, first_elements_8, first_elements_9, first_elements_10))
matrix_restuls1 <- t(res_MSE) #[,-c(4)]
left_column1 = c(5,6,7,8,9,10, 15, 20, 50)
final_matrix1 <- cbind(left_column1, matrix_restuls1[-c(1),])
colnames(final_matrix1) <- c("n","Phylofit", "Johnson", "New", "Unbiased for r", "Unbiased for 1/r" )

rownames(final_matrix1) <- NULL

data = as.data.frame(final_matrix1)
data_long <- melt(data, id.vars = "n", variable.name = "Method", value.name = "Value")


ggplot(data_long, aes(x = factor(n), y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Comparison of the Estimators for r = 1",
       x = "n",
       y = "Mean Absolute Deviation ") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 14)) +
  # scale_y_break(c(0.51, 0.8), scales = 0.01) +
  scale_y_continuous(
    breaks = seq(0, 6, by = 0.1),
    sec.axis = sec_axis(~., breaks = NULL)) +  # Verhindert die sekundäre Achse (rechte Seite)  )
  theme(legend.position = "top")  # Legende oben rechts


