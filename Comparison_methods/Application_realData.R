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

# Calculate new estimator
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

# Apply the different esimtation methods to the data
calc_estimator_realdata <- function(numInds, c){
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
      estimator_JS[i] = c * calc_estimator_adapted(n, h)
      
      # Estimator based on Johnson et al.
      estimator_old[i] = n/il
      
      # Phylofit
      test_mcmc <-birthDeathMCMC(tree,maxGrowthRate = 4,alpha = 0.05,verbose = TRUE,nChains = 4,nCores = 1,chainLength = 2000 )
      estimator_mcmc[i] <- test_mcmc$estimate
  }
  return(list(estimator_mcmc,estimator_old, estimator_JS ))
}

est_real = calc_estimator_realdata(2, c)


