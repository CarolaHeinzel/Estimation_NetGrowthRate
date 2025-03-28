# install.packages("pheatmap")
# This code plots the boxplots and the MSE

library(pheatmap)
# Boxplots
# Read the data
data <- read.table("C:\\Users\\carol\\OneDrive\\Dokumente\\table_test_t_40_r_1_n_5.txt", header = TRUE, sep = "\t") 
# determine the parameters
r = 1
n = 5
t = 40
png(filename = paste0("boxplot_new_all_t_",t,"_r_",r,"_n_",n , ".png"), width = 1200, height = 600)
dev.off()

# Werte für n, r und T definieren
n_values <- c(5, 6,7,8,9,10)  # Beispielwerte für n
r_values <- c(0.5)       # Beispielwerte für r
t_values <- c(40)    # Beispielwerte für T
save_path <- "C:\\Users\\carol\\OneDrive\\Desktop\\Promotion\\San Diego\\Results\\"  # Windows
# Iteration über alle Kombinationen von n, r und T
for (n in n_values) {
  for (r in r_values) {
    for (t in t_values) {
      data <- read.table(paste0("C:\\Users\\carol\\OneDrive\\Dokumente\\eva_hatr_100_t",t,"_r_",r, "_n_",n,".txt"), header = TRUE, sep = "\t") 
      
      # Dateiname generieren
      filename <- paste0(save_path, "boxplot_new_Tinfinity_", t, "_r_", r, "_n_", n, ".png")
      
      # PNG-Device öffnen
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
        data$list1,data$list2,  data$list3,
        #  # durrett, mcmc, meiner, JS, JS mit Fakotr, JS mit erwartungstreu
        names = c("Johnson", "Phylofit", "New"),  # Names of the boxplots
        col = c("red", "blue", "green", "purple", "orange"),             # Colors of the boxplots
        main = paste("Comparison of the  estimators for n =", n, ", r =", r, ", T =", t),  # Title of the plot
        ylab = "Estimators"                                              # Label for the y-axis
      ) 
      
      points(1, mean(data$list1), col = "black", pch = 19, cex = 1.5) # Mittelwert als roten Punkt einfügen
      points(2, mean(data$list2), col = "black", pch = 19, cex = 1.5) # Mittelwert als roten Punkt einfügen
      points(3, mean(data$list3), col = "black", pch = 19, cex = 1.5) # Mittelwert als roten Punkt einfügen
      #points(4, mean(data$list4), col = "black", pch = 19, cex = 1.5) # Mittelwert als roten Punkt einfügen
      #points(4, mean(data$list5), col = "black", pch = 19, cex = 1.5) # Mittelwert als roten Punkt einfügen
      #points(5, mean(data$list6)*c_CH, col = "black", pch = 19, cex = 1.5) # Mittelwert als roten Punkt einfügen
      
      abline(h = r, col = "black", lwd = 2, lty = 2) # Farbe, Breite und Typ der Linie
      
      # Speichern und Grafik-Device schließen
      dev.off()
      
      # Ausgabe zur Kontrolle
      print(paste("Plot gespeichert:", filename))
    }
  }
}

# MSE

# 1) Calculate the MSE for the nine estimators

mse <- function(x, r){
  mean((x - r)^2)
}

calc_mse <- function(r, n, t, res_alpha_all){
  # Read the data
  input = paste0("C:\\Users\\carol\\OneDrive\\Dokumente\\table_test_new_t",t,"_r_",r,"_n_",n , ".txt")
  print(input)
  input = paste0("C:\\Users\\carol\\OneDrive\\Dokumente\\table_test_t_",t,"_r_",r,"_n_",n,".txt")

  data = read.table(input, header = TRUE)
  print(data)
  # Calculate the MSE for every entry, i.e have 9 values for the MSE in the end!
  # res_alpha_all: JS2 Durrett JS  NL
  # table:"Durrett", "Phylofit", "NL", "JS, "JS 2", "JS MSE"
  alpha_JS2 = res_alpha_all[[1]][n-4]
  alpha_durrett = res_alpha_all[[2]][n-4]
  alpha_JS = res_alpha_all[[3]][n-4]
  alpha_NL = res_alpha_all[[4]][n-4]
  # MSE given estimators
  mse_values <- sapply(data, mse, r = r)
  print(mse_values)
  names(mse_values) <- c("Durrett", "Phylofit", "NL", "JS", "JS2", "JS2_MSE")
  #mse_values["JS2_MSE"] = mean((data$list5*alpha_JS2 - r)**2)
  mse_values["durrett_MSE"] = mean((data$list1*alpha_durrett - r)**2)
  mse_values["JS_MSE"] = mean((data$list4*alpha_JS - r)**2)
  mse_values["NL_MSE"] = mean((data$list3*alpha_NL - r)**2)
  
  return(mse_values)
}
test = calc_mse(0.5, 5, 40, res_alpha_all)
print(test)

first_elements <- sapply(temp, function(x) x[i])
first_elements_6 <- sapply(temp_6, function(x) x[i])
first_elements_7 <- sapply(temp_7, function(x) x[i])
first_elements_8 <- sapply(temp_8, function(x) x[i])
first_elements_9 <- sapply(temp_9, function(x) x[i])
first_elements_10 <- sapply(temp_10, function(x) x[i])
final_matrix <- t(cbind(first_elements, first_elements_6, first_elements_7, first_elements_8, first_elements_9, first_elements_10))
matrix_restuls1 <- final_matrix[,-c(4)]

# 2) Put the MSE in the correct format 
left_column1 = c(5,6,7,8,9,10)
final_matrix1 <- cbind(left_column1, matrix_restuls1)
colnames(final_matrix1) <- c("n","Phylofit", "Durrett", "Durrett_MSE","NL", "NL_MSE", "JS", "JS_MSE", "JS_2", "JS_2__MSE")
rownames(final_matrix1) <- NULL
data = as.data.frame(final_matrix1)

data_long <- melt(data, id.vars = "n", variable.name = "Method", value.name = "Value")

# Create Plot for the MSE
ggplot(data_long, aes(x = factor(n), y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Comparison of the Estimators for r = 0.5",
       x = "n",
       y = "MSE") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 14))

