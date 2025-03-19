library(gtools)
library(pracma) 
library(cloneRate)
library(Rmpfr)
library(reshape2) 
# This code compares the estimators \hat r, \hat r_{durrett} and \hat r_{Johnson} by using the MSE
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

  set.seed(n)
  for(i in 1:num_Rep){
    a = runif(1, min = r, max = r+1)
    test <- simUltra(a,b = a - r,cloneAge = t,n = n,nTrees = 1,precBits = 1000,addStem = FALSE,nCores = 1)
    res <- cloneRate::internalLengths(test[[1]], alpha = 0.05)
    il= res$sumInternalLengths
    h = test[[2]]
    # New estimator
    estimator_JS[i] = calc_estimator_adapted(n, h)[[1]]
    # Estimator based on Johnson et al.
    estimator_old[i] = n/il
    # Phylofit
    test_mcmc <-birthDeathMCMC(test[[1]],maxGrowthRate = 4,alpha = 0.05,verbose = TRUE,nChains = 4,nCores = 1,chainLength = 2000 )
    estimator_mcmc[i] <- test_mcmc$estimate
  }
  # Yubo, Phylofit, New estimator
  return(list(estimator_old, estimator_mcmc, estimator_JS))
}

t_new = repeat_simulation(100, 5, 0.5, 40)

# repeat the simulation for different r
# save the estimators in a .txt file
repeat_simulation_T <- function(n_rep, t, r_list, n){
  K = length(r_list)
  estimator_johnsons <- numeric(K)
  estimator <- c(K)
  estimator_mcmc <- c(K)
  
  i = 1
  for(r in r_list){
    temp <- repeat_simulation(n_rep, n, r, t)
    temp_df <- data.frame(
      list1 = temp[[1]],                    
      list2 = temp[[2]],                   
      list3 = temp[[3]],                    
    )
    output_file = paste0("eva_hatr_t",t,"_r_",r,"_n_",n , ".txt")
    write.table(
      temp_df,
      file = output_file,
      sep = "\t",          
      row.names = FALSE,    
      col.names = TRUE,      
    )
    
    if(n==5){
      c_CH = 4.769173
    }else if(n==6){
      c_CH = 9.833995
    }else if(n==7){
      c_CH = 17.367272
    }else if(n==8){
      c_CH = 26.940960
    }else if(n==9){
      c_CH = 36.827248
    }else if(n==10){
      c_CH = 51.135133
    }
      
    mse1 = mean((temp[[1]]-r)^2)
    mse2 = mean((temp[[2]]-r)^2)
    mse3 = mean((temp[[3]]*c_CH-r)^2)    
    
    estimator[i] = mse3 # New Estimator
    estimator_mcmc[i] = mse2 # Phylofit
    estimator_johnson[i] = mse1 # Old estimator

    i = i+1
  }
  return(list(estimator_johson, estimator_mcmc, estimator))
}

temp_new = repeat_simulation_T(10, 40, c(0.5), 6) #22
print(temp_new)


data_test_new <- read.table(paste0("C:\\Users\\carol\\OneDrive\\Dokumente\\table_test_new_t40_r_0.5_n_6.txt"), header = TRUE, sep = "\t") 

temp = repeat_simulation_T(100, 40, c(0.5, 1), 5) 
temp_6 = repeat_simulation_T(100, 40, c(0.5, 1), 6) 
temp_7 = repeat_simulation_T(100, 40, c(0.5, 1), 7) 
temp_8 = repeat_simulation_T(100, 40, c(0.5, 1), 8)
temp_9 = repeat_simulation_T(100, 40, c(0.5, 1), 9)
temp_10 = repeat_simulation_T(100, 40, c(0.5, 1),10)

# Plot the MSE                        
i = 2
first_elements <- sapply(temp, function(x) x[i])
first_elements_6 <- sapply(temp_6, function(x) x[i])
first_elements_7 <- sapply(temp_7, function(x) x[i])
first_elements_8 <- sapply(temp_8, function(x) x[i])
first_elements_9 <- sapply(temp_9, function(x) x[i])
first_elements_10 <- sapply(temp_10, function(x) x[i])
final_matrix <- t(cbind(first_elements, first_elements_6, first_elements_7, first_elements_8, first_elements_9, first_elements_10))
matrix_restuls1 <- final_matrix 

left_column1 = c(5,6,7,8,9,10)
final_matrix1 <- cbind(left_column1, matrix_restuls1)
colnames(final_matrix1) <- c("n","Old", "Phylofit", "New", )

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

