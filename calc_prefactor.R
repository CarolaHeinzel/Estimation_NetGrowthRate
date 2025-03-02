source("C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\main.R")


calc_prefactor <- function(n, num_Rep, t, r){
  ew_CH <- c(num_Rep)
  for(i in 1:num_Rep){
    a = runif(1, min = r, max = r+1)
    test <- simUltra(a,b = a - r,cloneAge = t,n = n,nTrees = 1,precBits = 1000,addStem = FALSE,nCores = 1)
    h = test[[2]]
    ew_CH[i] = calc_estimator_adapted(n, h)[[2]]
  }
  v = ew_CH
  y <- rep(r, num_Rep) 
  alpha_opt <- sum(v * y) / sum(v^2)
  return(alpha_opt)
}

rep_calc_prefactor <- function(n_range, num_Rep, t, r){
  res_alpha = rep(0, length(n_range))
  i = 1
  for(n in n_range){
    res_alpha[i] = calc_prefactor(n, num_Rep, t, r)
    i = i+1
  }
  return(res_alpha)
}


res_alpha_all = rep_calc_prefactor(c(5,6,7,8,9,10), 1000, 100, 1)
print(res_alpha_all)