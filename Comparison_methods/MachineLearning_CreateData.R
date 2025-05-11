repeat_simulation_train_data <- function(num_Rep, r, n, t){
res_ml <- data.frame(
    col1 = numeric(0), 
    col2 = numeric(0), 
    col3 = numeric(0), 
    col4 = numeric(0), 
    col5 = numeric(0), 
    col6 = numeric(0)
  )
  
  # Neue Zeile zu DataFrame hinzufügen
  estimator_JS <- c(num_Rep)
  #set.seed(3)
  #set.seed(2)
  for(i in 1:num_Rep){
    a = runif(1, min = r, max = r+1)
    test <- simUltra(a,b = a - r,cloneAge = t,n = n,nTrees = 1,precBits = 1000,addStem = FALSE,nCores = 1)
    res <- cloneRate::internalLengths(test[[1]], alpha = 0.05)
    il= res$sumInternalLengths
    il_ex = res$sumExternalLengths 
    h = test[[2]]
    res_ml <- rbind(res_ml, c(il, il_ex, h))
    # New estimator
    estimator_JS[i] = calc_estimator_adapted(n, h)#[[1]]

  }
  #filename <- paste0("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\Simulation_infinity\\CPPs_eva_4_",t,"_n_", n, "_r_", r, ".txt")
  #est_all = cbind(estimator_JS, estimator_old, estimator_mcmc)
  #write.table(est_all, file = filename, row.names = FALSE, col.names = FALSE)
  
  # Yubo, Phylofit, New estimator
  return(list(res_ml, estimator_JS))
}

test = repeat_simulation_train_data(1000, 1, 5, 40)


create_train <- function(numRep, r_vec, t, n){
  i = 0
  for(r in r_vec){
    test = repeat_simulation_train_data(numRep, r,  n, t)
    my_list = test[[2]]
    quantile_90 <- quantile(my_list, 0.9)
    # Positionen, die größer sind als das 90%-Quantil
    positions <- which(my_list > quantile_90)
    df = test[[1]]
    df$new_column <- 0
    # Setze in den Positionen mit 1
    df$new_column[positions] <- 1
    colnames(df) <- c("a", "b", "c", "d", "e", "f", "Class")
    if(i == 0){
      df_all = df
      i = i+1
    }else{
      df_all = rbind(df_all, df)
    }
  }
  return(df_all)
}


test1 = create_train(1000, c(0.8), 40, 5)


write.csv(test1, "C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\ML\\Test_Data_ml.csv", row.names = FALSE)


ape::plot.phylo(test[[1]],
                direction = "downwards", show.tip.label = FALSE
)
