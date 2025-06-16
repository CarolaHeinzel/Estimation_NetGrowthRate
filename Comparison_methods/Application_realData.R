library(cloneRate)

data(realCloneData)
data(longitudinalData)


ape::plot.phylo(cloneRate::realCloneData[["fullTrees"]][["PD34493"]],
                   direction = "downwards", show.tip.label = FALSE
                 )
n_ind = length(realCloneData[["fullTrees"]]) # Number of Individuals
max_individuals = rep(0, n_ind)

for(i in 1:n_ind){
  max_temp = Ntip(realCloneData[["cloneTrees"]][[i]])
  max_individuals[i] = max_temp
}
library(ape)
# Calculate Coalescence Times
ind = 1
tree =  realCloneData[["cloneTrees"]][[ind]]
coal_times <- node.depth.edgelength(tree)

# Print Coalescence Times
coal_times[(Ntip(realCloneData[["cloneTrees"]][[ind]])+1):(2*realCloneData[["cloneTrees"]][[ind]]$Nnode+1)]


# Apply the different esimtation methods to the data
calc_estimator_realdata <- function(numInds){
  estimator_JS = rep(0, numInds)
  estimator_mcmc = rep(0, numInds)
  estimator_old = rep(0, numInds)
  
  for(i in 1:numInds){
      i = 5
      tree = realCloneData[["cloneTrees"]][[i]] # Real Data
      res = cloneRate::internalLengths(tree, alpha = 0.05)
      n = Ntip(tree)
      il= res$sumInternalLengths
      coal_times = node.depth.edgelength(tree)[]
      h = coal_times[(Ntip(realCloneData[["cloneTrees"]][[ind]])+1):(2*realCloneData[["cloneTrees"]][[ind]]$Nnode+1)]
      
      # New estimator
      estimator_JS[i] = calc_estimator_adapted(n, h)
      
      # Estimator based on Johnson et al.
      estimator_old[i] = n/il
      
      # Phylofit
      test_mcmc <-birthDeathMCMC(tree,maxGrowthRate = 4,alpha = 0.05,verbose = TRUE,nChains = 4,nCores = 1,chainLength = 2000 )
      estimator_mcmc[i] <- test_mcmc$estimate
  }
  return(list(estimator_mcmc,estimator_old, estimator_JS ))
}

est_real = calc_estimator_realdata(2)
print(est_real[[3]][5]*161.7437)


# Plot the results
library(ggplot2)

# 0.1371847, 0.145, 0.196
# Erstelle die Daten
data <- data.frame(
  Method = rep(c("Phylofit", "Johnson", "New"), each = 2),  # Die Gruppen Phylofit, Johnson, New
  Wert = c(0.145, 0.2650246 , 0.196, 0.2944672, 0.1371847, 0.2375439),  # Werte
  Kategorie = rep(c(" PD34493", " PD41305"), 3)  # Kategorien
)

# Erstelle das Balkendiagramm
ggplot(data, aes(x = Kategorie, y = Wert, fill = Method, group = Method)) +
  geom_bar(stat = "identity", position = "dodge") +  # 'dodge' für nebeneinander liegende Balken
  labs(x = "Individual", y = "Estimator") +
  theme_minimal() +
  scale_fill_manual(values = c("Phylofit" = "lightblue", "Johnson" = "lightgreen", "New" = "green"))  # Farben für die Gruppen
