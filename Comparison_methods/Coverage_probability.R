filename_1 <- paste0("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\Simulation_infinity\\CPPs_infinity_1_5_r_1.txt")

res_coverage <- rep(0, 8)
r = 1
i = 1
for(n in c(5,6,7,8,9,10, 15, 20)){
  
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
  }else if(n == 50){
    c_n =  284.5928
    c_1 = 296.7941
    c_2 = 19*18*1199057081/1326917592
  }
  filename_2 <- paste0("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\Simulation_infinity\\CPPs_infinity_1_",n,"_r_1.txt")
  
  daten1 <- read.csv(filename_2, header = FALSE)
  print(nrow(daten1))
  # Calculate the quantiles for the data
  q = quantile(daten1$V1, probs = c(0.025, 0.975))
  if(n <= 10){
    filename <- paste0("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\Simulation_infinity\\CPPs_eva_40_n_",n,"_r_",r,".txt")
    
  }else{
    filename <- paste0("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\Simulation_infinity\\CPPs_final_n",n,"_r_",r,"_T40.txt")
    
  }
  
  # Datei einlesen (ersetze 'datei.csv' durch den tatsächlichen Dateinamen)
  daten <- read.table(filename, header = FALSE, sep = "", col.names = c("X1", "X2", "X3"))
  print(nrow(daten))
  q_pf = q["2.5%"]/mean(daten1$V1)
  q3_pf = q["97.5%"]/mean(daten1$V1)
  q1 = daten$X1*q_pf*c_n
  q3 = daten$X1*q3_pf*c_n
  print(q_pf)
  print(q3_pf)
  t_true = r
  covered <- (t_true >= q1) & (t_true <= q3)
  
  # Empirische Coverage Rate
  coverage_rate <- mean(covered)
  print(paste("Empirical coverage rate:", round(coverage_rate, 3)))
  
  res_coverage[i] = coverage_rate
  i=i+1
}

# Rahmenfarbe der Balken


barplot(res_coverage - 0.7,
        ylim = c(0, 0.3),                      # y-Achse startet bei 0.7
        names.arg = c(5,6,7,8,9,10,15,20),     # x-Achse
        xlab = "n",
        ylab = "Coverage Probability",
        col = "blue",
        border = "black",
        space = 0.2,
        axes = FALSE)                         # Achsen deaktivieren

# Manuelle y-Achse hinzufügen, beginnend bei 0.7
axis(2, at = seq(0, 0.3, by = 0.05),            # Achsenwerte von 0 bis 0.3
     labels = seq(0.7, 1, by = 0.05),            # Zeigt auf der Y-Achse 0.7 bis 1
     las = 1)  







