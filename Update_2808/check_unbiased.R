library(gtools)
library(pracma) 
library(cloneRate)
library(Rmpfr)
# This script implements estimators for 1/r
# Especially, it checks whether the estimators are unbiased for 1/r and whether the calculation 
# of c(n) is correct
source("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\own_idea.R")
source("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Simulation_self.R")

# Important to calculate the expected value according to Durrett
ew_durrett <- function(n){
  sum_result <- 0
  for (k in 2:(n - 1)) {
    value <-  n / (k*(k-1))
    sum_result <- sum_result + value 
  }
  return(sum_result)
}

compute_list_properties <- function(liste, n) {
  normal_list = liste
  max_value <- max(normal_list)
  mean_value <- mean(normal_list)
  
  # ∑_i ∑_j (H_i - H_j)^+
  double_sum <- sum(outer(normal_list, normal_list, function(x, y) pmax(0, x - y)))
  
  return(list(Maximum = max_value, Mittelwert = mean_value, d_sum = double_sum))
}

# Calculates \hat r_{JS}
calc_estimator_unbiased <- function(n, h){
  values <- compute_list_properties(h, n)
  est <-  1/(n-1)*values$d_sum + values$Maximum - values$Mittelwert 
  ev_th <- ew_durrett(n)
  return(est/ev_th)
}

# Calculates \hat r_{JS, 2}
calc_estimator_adapted_unbiased <- function(n, h){
  if(n==5){
    c_n <- 12*115/144
  }else if(n == 6){
    c_n <- 5*4*163/200
  }else if(n == 7){
    c_n <- 6*5*497/600
  }else if(n == 8){
    c_n <- 7*6*617/735
  }else if(n == 9){
    c_n <- 8*7*13311/15680
  }else if(n == 10){
    c_n <- 9*8*15551/18144
  }else{
    c_n = 1
  }
  values <- compute_list_properties(h, n)
  est <-  values$d_sum
  return(est/c_n)
  
}

# repeat the simulation num_Rep times
repeat_simulation <- function(num_Rep, n, r, t){
  estimator <- numeric(num_Rep)
  estimator_mcmc <- c(num_Rep)
  estimator_new <- c(num_Rep)
  estimator_JS <- c(num_Rep)
  estimator_JS_adapted <- c(num_Rep)
  ew_durrett <- ew_durrett(n)
  ew_CH <- c(num_Rep)
  set.seed(22)
  for(i in 1:num_Rep){
    a = runif(1, min = r, max = r+1)
    test <- simUltra(a,b = a - r,cloneAge = t,n = n,nTrees = 1,precBits = 1000,addStem = FALSE,nCores = 1)
    res <- cloneRate::internalLengths(test[[1]], alpha = 0.05)
    il= res$sumInternalLengths
    h = test[[2]]
    estimator_JS[i] = calc_estimator_unbiased(n, h)
    estimator_JS_adapted[i] = calc_estimator_adapted_unbiased(n, h)
    res_new =  internalLengths(test[[1]], alpha = 0.05)
    estimator[i] <-  il/ew_durrett # durrett
  }
  # durrett, \hat r_{JS}, \hat r_{new, 2}
  return(list(estimator, estimator_JS, estimator_JS_adapted))
}

# repeat the simulation for different r
# save the estimators in a .txt file
repeat_simulation_T <- function(n_rep, t, r_list, n){
  i = 1
  for(r in r_list){
    temp <- repeat_simulation(n_rep, n, r, t)
    temp_df <- data.frame(
      list1 = temp[[1]],                    
      list2 = temp[[2]],                   
      list3 = temp[[3]]               
    )
    output_file = paste0("table_test_unbiased_t",t,"_r_",r,"_n_",n , ".txt")
    write.table(
      temp_df,
      file = output_file,
      sep = "\t",          
      row.names = FALSE,    
      col.names = TRUE,      
    )
    i = i+1
  }
  return("Finished")
}
temp_newunbiased_6 = repeat_simulation_T(1000, 100, c(0.5, 2), 6) #22
temp_new_unbiased_5 = repeat_simulation_T(1000, 100, c(0.5, 2), 5) #22
temp_new_unbiased_10 = repeat_simulation_T(1000, 100, c(0.5, 2), 10) #22
temp_new_unbiased_7 = repeat_simulation_T(1000, 100, c(0.5, 2), 7) #22
temp_new_unbiased_8 = repeat_simulation_T(1000, 100, c(0.5, 2), 8) #22
temp_new_unbiased_9 = repeat_simulation_T(1000, 100, c(0.5, 2), 9) #22

data_test_new_unbiased_5 <- read.table(paste0("C:\\Users\\carol\\OneDrive\\Dokumente\\table_test_unbiased_t100_r_0.5_n_5.txt"), header = TRUE, sep = "\t") 
data_test_new_unbiased_6 <- read.table(paste0("C:\\Users\\carol\\OneDrive\\Dokumente\\table_test_unbiased_t100_r_0.5_n_6.txt"), header = TRUE, sep = "\t") 
data_test_new_unbiased_7 <- read.table(paste0("C:\\Users\\carol\\OneDrive\\Dokumente\\table_test_unbiased_t100_r_0.5_n_7.txt"), header = TRUE, sep = "\t") 
data_test_new_unbiased_8 <- read.table(paste0("C:\\Users\\carol\\OneDrive\\Dokumente\\table_test_unbiased_t100_r_0.5_n_8.txt"), header = TRUE, sep = "\t") 
data_test_new_unbiased_9 <- read.table(paste0("C:\\Users\\carol\\OneDrive\\Dokumente\\table_test_unbiased_t100_r_0.5_n_9.txt"), header = TRUE, sep = "\t") 
data_test_new_unbiased_10 <- read.table(paste0("C:\\Users\\carol\\OneDrive\\Dokumente\\table_test_unbiased_t100_r_0.5_n_10.txt"), header = TRUE, sep = "\t") 

mean(data_test_new_unbiased_5$list1)
mean(data_test_new_unbiased_5$list2)
mean(data_test_new_unbiased_5$list3)

mean(data_test_new_unbiased_10$list1)
mean(data_test_new_unbiased_10$list2)
mean(data_test_new_unbiased_10$list3)

n_values <- 5:10  # values for n
means_list1 <- c()
means_list2 <- c()
means_list3 <- c()

for (n in n_values) {
  data_name <- paste0("data_test_new_unbiased_", n)  
  data <- get(data_name)  
  means_list1 <- c(means_list1, mean(data$list1))
  means_list2 <- c(means_list2, mean(data$list2))
  means_list3 <- c(means_list3, mean(data$list3))
}

results_matrix <- matrix(c(means_list1, means_list2, means_list3), 
                         nrow = length(n_values), 
                         ncol = 3, 
                         byrow = FALSE)

# Werte definieren (3 Werte pro x-Wert)
values <- t(results_matrix)
colnames(values) <- NULL
rownames(values) <- NULL

# X-Achsen-Beschriftungen
x_labels <- c("5", "6", "7", "8", "9", "10")

# Farben für die Gruppen
colors <- c("skyblue", "tomato", "seagreen")

# Balkendiagramm zeichnen (gruppiert)
barplot(values, beside = TRUE, names.arg = x_labels, ylab = "Empirical Mean",xlab = "n",col = colors, border = "black", ylim = c(0, max(values) + 0.5))

# Gestrichelte schwarze Linie bei y = 1
abline(h = 0.5, lty = 2, col = "black")

# Legende hinzufügen
legend("topright", legend = c("Durrett", "JS", "JS 2"), fill = colors, border = "black")


temp_new_unbiased_5_1000 = repeat_simulation_T(1000, 100, c(0.5), 5) #22
data_test_new_unbiased_5_1000 <- read.table(paste0("C:\\Users\\carol\\OneDrive\\Dokumente\\table_test_unbiased_t100_r_0.5_n_5.txt"), header = TRUE, sep = "\t") 
mean(data_test_new_unbiased_5_1000$list1)


#%%%
# Evaluate whether the caluation of c(n) is correct:
# Approximate the expected value of E((H_{1,n,\infty} - H_{2,n,\infty})^+)
# by the empirical mean

calc_e <- function(n, numRep, r, t){
  res_e = 0
  for(i in 1:numRep){
    # Simulate the CPP
    a = runif(1, min = r, max = r+1)
    test <- simUltra(a,b = a - r,cloneAge = t,n = n,nTrees = 1,precBits = 1000,addStem = FALSE,nCores = 1)
    h = test[[2]]
    # print(h)
    # Calc (h_1-h_2)^+
    temp = h[1] - h[2]
    if(temp > 0){
      res_e = res_e + temp
    }
  }
  return(res_e/numRep)
}
r = 1
example_ev = calc_e(5, 1000, r, 100)
example_ev_6 = calc_e(6, 1000, r, 100)
example_ev_7 = calc_e(7, 1000, r, 100)
example_ev_8 = calc_e(8, 1000, r, 100)
example_ev_9 = calc_e(9, 1000, r, 100)
example_ev_10 = calc_e(10, 1000, r, 100)

print(example_ev)


# Werte definieren (3 Werte pro x-Wert)
values <- matrix(c(example_ev, 1.4097222222222223, 115/144, 
                   example_ev_6, 1.3358333333333332, 163/200, 
                   example_ev_7, 1.2850000000000001, 497/600, 
                   example_ev_8, 1.2477891156462584, 617/735, 
                   example_ev_9, 1.2193239795918367, 13311/15680,
                   example_ev_10, 1.196819885361552, 15551/18144), 
                 nrow = 3, byrow = FALSE)

# X-Achsen-Beschriftungen
x_labels <- c("5", "6", "7", "8", "9", "10")

# Farben für die Gruppen
colors <- c("skyblue", "tomato", "seagreen")
par(mgp = c(3, 1, 0)) 
# Balkendiagramm zeichnen (gruppiert)
barplot(values, beside = TRUE, names.arg = x_labels, xlab = "n", ylab = expression(E((H[1,n,inf] - H[2,n,inf])^"+")), col = colors, border = "black", ylim = c(0, max(values) + 0.5))


# Legende hinzufügen
legend("topright", legend = c("Empirical Mean", "JS", "CH"), fill = colors, border = "black")
