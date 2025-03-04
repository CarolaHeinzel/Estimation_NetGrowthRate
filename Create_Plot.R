library(ggplot2)
library(reshape2)

# Example Data
data <- data.frame(
  n = c(5, 6),
  Durrett = c(0.3235,1),
  Phylofit = c(0.2607, 0.4387),
  NL = c(0.2076, 0.8304),
  JS2 = c(0.2684, 1.0738),
  JSMSE = c(0.0979, 0.3918)
)

data = as.data.frame(final_matrix1)

data_long <- melt(data, id.vars = "n", variable.name = "Method", value.name = "Value")

# Create Plot for the MSE
ggplot(data_long, aes(x = factor(n), y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Comparison of the Estimators for r = 1",
       x = "n",
       y = "MSE") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 14))

# Create the boxplots
n_values <- c(5, 6,7,8,9,10)  
r_values <- c(1, 0.5)      
t_values <- c(40, 40)    
save_path <- "C:\\Users\\Results\\" 
for (n in n_values) {
  for (r in r_values) {
    for (t in t_values) {
      # Data with the estimators
      data <- read.table(paste0("table_test_t_",t,"_r_",r, "_n_",n,".txt"), header = TRUE, sep = "\t") 
      
      # Dateiname of the plot
      filename <- paste0(save_path, "boxplot_new_", t, "_r_", r, "_n_", n, ".png")
      png(filename = filename, width = 800, height = 600)
      
      par(mar = c(5, 6, 4, 2))
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
      boxplot(
        data$list1,data$list2,  data$list4,  data$list5, data$list6* c_CH,
        names = c("Durrett", "Phylofit", "JS", "JS 2", "JS MSE"),  # Names of the boxplots
        col = c("lightblue", "lightgreen", "yellow", "orange", "grey"),             # Colors of the boxplots
        main = paste("Comparison of the  estimators for n =", n, ", r =", r, ", T =", t),  # Title of the plot
        ylab = "Estimators"                                              # Label for the y-axis
      ) 
      
      points(1, mean(data$list1), col = "black", pch = 19, cex = 1.5) 
      points(2, mean(data$list2), col = "black", pch = 19, cex = 1.5) 
      #points(3, mean(1/data$list3), col = "black", pch = 19, cex = 1.5)
      points(3, mean(data$list4), col = "black", pch = 19, cex = 1.5) 
      points(4, mean(data$list5), col = "black", pch = 19, cex = 1.5) 
      points(5, mean(data$list6)*c_CH, col = "black", pch = 19, cex = 1.5) 
      
      abline(h = r, col = "black", lwd = 2, lty = 2) 
      
    
      dev.off()
      
    }
  }
}
