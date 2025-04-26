# Determines the number of digits of the constants for which we are sure that they are correct.

# Load data
n = 5
for(n in c(6,7,8,9,10,11,12,13,14,15,16, 17, 18, 19, 20)){
  input_file <- paste0("/home/ch1158/CPPs_infinity_", n, "_r_1.txt")
  
  # Load Data
  data <- scan(input_file, what = numeric())
  
  # Divide data 
  
  chunk_size <- 1e6
  
  K  = 100 # Number of repetitions
  r = 1
  res <- c()
  for (k in 1:K) {
    # Randomly select data
    shuffled_data <- sample(data)
    
    
    for (i in 1:10) {
      start_idx <- (i - 1) * chunk_size + 1
      end_idx <- i * chunk_size
      chunk <- shuffled_data[start_idx:end_idx]
      
      # Calculate the optimal values
      v = chunk
      y <- rep(r, chunk_size) 
      alpha_opt <- sum(v * y) / sum(v^2)
      res = c(res, alpha_opt)
    }
  }
  print(min(res/((n-1)*(n-2))))
  print(max(res/((n-1)*(n-2))))
}
