library(gtools)
library(pracma) 
library(cloneRate)
library(Rmpfr)
library(reshape2) 
source("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Simulation_self.R")


compute_list_properties <- function(liste, n) {
  normal_list = liste

  # ∑_i ∑_j (H_i - H_j)^+
  double_sum <- sum(outer(normal_list, normal_list, function(x, y) pmax(0, x - y)))
  
  return(double_sum)
}

calc_estimator_adapted <- function(n, h){
  values <- compute_list_properties(h, n)
  return(1/values)
}

# repeat the simulation num_Rep times
repeat_simulation <- function(num_Rep, n, r, t){
  estimator_mcmc <- c(num_Rep)
  estimator_old <- c(num_Rep)
  estimator_JS <- c(num_Rep)
  set.seed(r+1)
  for(i in 1:num_Rep){
    print(i)
    a = runif(1, min = r, max = r+1)
    test <- simUltra(a,b = a - r,cloneAge = t,n = n,nTrees = 1,precBits = 1000,addStem = FALSE,nCores = 1)
    
    res <- cloneRate::internalLengths(test[[1]], alpha = 0.05)
    il= res$sumInternalLengths
    h = test[[2]]
    
    # New estimator
    estimator_JS[i] = calc_estimator_adapted(n, h)#[[1]]

    # Estimator based on Johnson et al.
    estimator_old[i] = 1/il
    
    # Phylofit
    test_mcmc <-birthDeathMCMC(test[[1]],maxGrowthRate = 4,alpha = 0.05,verbose = TRUE,nChains = 4,nCores = 1,chainLength = 2000 )
    estimator_mcmc[i] <- test_mcmc$estimate
  }
  filename <- paste0("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\Simulation_infinity\\CPPs_eva_",t,"_n_", n, "_r_", r, ".txt")
  est_all = cbind(estimator_JS, estimator_old, estimator_mcmc)
  write.table(est_all, file = filename, row.names = FALSE, col.names = FALSE)
  
  # Yubo, Phylofit, New estimator
  return(list(estimator_old, estimator_mcmc, estimator_JS))
}
system.time({
t_new_1 = repeat_simulation(100, 50, 1, 40)
})

system.time({
  t_new_1_6 = repeat_simulation(19, 100, 1, 40)
})
system.time({
  t_new_1_7 = repeat_simulation(1000, 7, 1, 40)
})

system.time({
  t_new_8_1 = repeat_simulation(1000, 8, 1, 40)
})

system.time({
  t_new_9_1 = repeat_simulation(100, 9, 1, 40)
})

system.time({
  t_new_10_1= repeat_simulation(100, 10, 1, 40)
})

system.time({
  t_new_10_1= repeat_simulation(100, 15, 1, 40)
})

system.time({
  t_new_10_1= repeat_simulation(100, 20, 1, 40)
})

system.time({
  t_new_8 = repeat_simulation(1000, 9, 0.5, 40)
})


system.time({
  t_new_9 = repeat_simulation(100, 15, 0.5, 40)
})

system.time({
  t_new_10 = repeat_simulation(100, 20, 0.5, 40)
})
# Calculate the MSE
print(mean((t_new[[1]]-0.5)^2))
print(mean((t_new[[2]]-0.5)^2))
print(mean((t_new[[3]]*5.574809 -0.5)^2))

# repeat the simulation for different r
# save the estimators in a .txt file
repeat_simulation_T <- function(n_rep, t, r_list, n){
  # include the prefactors
  res_pre = read.table("C:\\Users\\carol\\OneDrive\\Dokumente\\prefactors_infty.txt")
  K = length(r_list)
  estimator_JS2 <- numeric(K)
  estimator_old <- c(K)
  estimator_mcmc <- c(K)
  i = 1
  c_JS2 = as.numeric(res_pre[n-3,1])
  
  for(r in r_list){
    temp <- repeat_simulation(n_rep, n, r, t)
    temp_df <- data.frame(
      list1 = temp[[1]],                    
      list2 = temp[[2]],                   
      list3 = temp[[3]]*c_JS2                   
    )
    # in list6, we have to include the prefactor afterwards
    output_file = paste0("eva_hatr_100_t",t,"_r_",r,"_n_",n , ".txt")
    write.table(
      temp_df,
      file = output_file,
      sep = "\t",          
      row.names = FALSE,    
      col.names = TRUE,      
    )
    mse1 = mean((temp[[1]]-r)^2)
    mse2 = mean((temp[[2]]-r)^2)
    mse3 = mean((temp[[3]]*c_JS2-r)^2)
    print(c_JS2)
    estimator_old[i] = mse1 # Old
    estimator_mcmc[i] = mse2 # Phylofit
    estimator_JS2[i] = mse3 # JS 2 MSE
    i = i+1
  }
  return(list(estimator_old, estimator_mcmc, estimator_JS2))
}
# Hierfür: funktioniert für JS idee am besten!
temp_new = repeat_simulation_T(10, 40, c(0.5), 6) #22
print(temp_new)


data_test_new <- read.table(paste0("C:\\Users\\carol\\OneDrive\\Dokumente\\eva_hatr_100_t40_r_0.5_n_5.txt"), header = TRUE, sep = "\t") 
print(mean((data_test_new$list3*4.35/5.57-0.5)^2))



y_vals <- rep(0, 100)
x_vals <- seq(2, 8, length.out = 100)
for(i in 1:100){
  y_vals[i] = mean((data_test_new$list3*x_vals[i]/5.574809-0.5)^2)
}
data_test_new_r1 <- read.table(paste0("C:\\Users\\carol\\OneDrive\\Dokumente\\eva_hatr_100_t40_r_1_n_5.txt"), header = TRUE, sep = "\t") 

y_vals_g = rep(0, 100)
for(i in 1:100){
  y_vals_g[i] = mean((data_test_new_r1$list3*x_vals[i]/5.574809-1)^2)
}

par(mar = c(5, 6, 4, 2) + 0.1) 
plot(x_vals, y_vals, type = "l", col = "red", lwd = 2, xlab = "c(n)", 
     ylab = expression(MSE[c(n)](hat(r), r)), 
     main = "n = 5", ylim = c(0, 0.55))
lines(x_vals, y_vals_g, col = "blue", lwd = 2)  # Zweite Kurve hinzufügen
abline(v = 5.574809, col = "black", lwd = 2, lty = 2)  # Gepunktete Linie
legend("topright", legend = c("r = 0.5", "r = 1"), col = c("red", "blue"), lwd = 2)

#start_time <- Sys.time()
temp = repeat_simulation_T(100, 40, c(0.5), 5) #n+1 ist der seed
#end_time <- Sys.time()
#print(end_time - start_time)
temp_6 = repeat_simulation_T(100, 40, c(0.5), 6) 
temp_7 = repeat_simulation_T(100, 40, c(0.5), 7) 
temp_8 = repeat_simulation_T(100, 40, c(0.5), 8)
temp_9 = repeat_simulation_T(100, 40, c(0.5), 9)
temp_10 = repeat_simulation_T(100, 40, c(0.5),10)

i = 1
first_elements <- sapply(temp, function(x) x[i])
first_elements_6 <- sapply(temp_6, function(x) x[i])
first_elements_7 <- sapply(temp_7, function(x) x[i])
first_elements_8 <- sapply(temp_8, function(x) x[i])
first_elements_9 <- sapply(temp_9, function(x) x[i])
first_elements_10 <- sapply(temp_10, function(x) x[i])
final_matrix <- t(cbind(first_elements, first_elements_6, first_elements_7, first_elements_8, first_elements_9, first_elements_10))
matrix_restuls1 <- final_matrix #[,-c(4)]

left_column1 = c(5,6,7,8,9,10)
final_matrix1 <- cbind(left_column1, matrix_restuls1)
colnames(final_matrix1) <- c("n","Johson", "Phylofit", "New" )

rownames(final_matrix1) <- NULL

data = as.data.frame(final_matrix1)

data_long <- melt(data, id.vars = "n", variable.name = "Method", value.name = "Value")

# Create Plot
ggplot(data_long, aes(x = factor(n), y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Comparison of the Estimators for r = 0.5",
       x = "n",
       y = "MSE") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 14))


