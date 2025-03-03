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

# Create Plot
ggplot(data_long, aes(x = factor(n), y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Comparison of the Estimators for r = 1",
       x = "n",
       y = "MSE") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +
  theme(text = element_text(size = 14))
