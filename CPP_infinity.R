# Simulate ultrametric birth and death branching trees for T = \infty
library(Rmpfr)

simUltra_infty <- function(a, b, n, nTrees = 1,
                     precBits = 1000, addStem = FALSE, nCores = 1) {
  # Store runtime for each tree
  ptm <- proc.time()
  # Convert params to high precision mpfr
  a_mpfr <- mpfr(a, precBits)
  b_mpfr <- mpfr(b, precBits)
  net_mpfr <- a_mpfr - b_mpfr
  n_mpfr <- mpfr(n, precBits)
  one <- mpfr(1, precBits)
  
  # Draw Y = y from the inverse CDF
  uniform_rv <- mpfr(stats::runif(1, min = 0, max = 1), precBits)
  y_mpfr <- uniform_rv**(1/n)/(1-uniform_rv**(1/n))
  # Generate the coalescence times
  coal_times_mpfr <- sapply(rep(y_mpfr, n - 1), inv_cdf_coal_times_inf,
                            net = net_mpfr,
                            a = a_mpfr, precBits = precBits
  )
  
  # Convert back to normal numeric (no longer need high precision)
  coal_times <- suppressWarnings(sapply(coal_times_mpfr, Rmpfr::asNumeric))
  
  # Convert coal times into tree by randomly merging lineages
  tree <- coal_to_tree(coal_times)
  
  # Add stem starting the tree from zero, rooting the tree appropriately
  if (addStem) {
    tree$edge[tree$edge > n] <- tree$edge[tree$edge > n] + 1
    tree$edge <- rbind(c(n + 1, n + 2), tree$edge)
    tree$edge.length <- c(max(coal_times), tree$edge.length)
    tree$Nnode <- tree$Nnode + 1
  }
  
  # Add metadata for making the tree
  runtime <- proc.time()[["elapsed"]] - ptm[["elapsed"]]
  tree$metadata <- data.frame(
    "r" = a - b, "a" = a, "b" = b,
    "n" = n, "runtime_seconds" = runtime, "addStem" = addStem
  )
  
  # Return the tree created from the coalescence times drawn from Lambert distribution
  return(tree)
}
  

tree = simUltra_infty(2,1,5)
print(tree)



coal_to_tree <- function(coal_times) {
  # coal_times must be a vector of numbers
  if (!inherits(coal_times, "numeric") | length(coal_times) < 2) {
    stop("coal_times input to coal_to_tree() function must a numeric vector")
  }
  
  # Get number of tips, n
  n <- length(coal_times) + 1
  
  # Set number of total edges and initialize edge and edge.length
  numEdges <- 2 * n - 2
  edge.length <- rep(0, numEdges)
  edge <- matrix(NA, nrow = numEdges, ncol = 2)
  
  # Fill heights vec to keep track of height of nodes
  heights <- rep(0, numEdges - 1)
  
  # Sort coalescence times (smallest to largest)
  coal_times_sorted <- sort(coal_times, decreasing = FALSE)
  possibleChildren <- as.integer(c(1:n))
  currentNode <- as.integer(2 * n - 1)
  
  # Loop through n-1 internal nodes
  for (i in 1:(n - 1)) {
    # Sample the children
    children <- sample(possibleChildren, size = 2, replace = FALSE)
    
    # Go to next open row
    row <- which(is.na(edge[, 1]))[c(1, 2)]
    
    # Fill second column with children and first with node
    edge[row, 2] <- children
    edge[row, 1] <- currentNode
    
    # Set edge.length as diff. between coal time and children height
    edge.length[row] <- coal_times_sorted[i] - heights[children]
    
    # Set height of current node
    heights[currentNode] <- coal_times_sorted[i]
    
    # Add current node to list of possible children
    possibleChildren <- c(possibleChildren[!possibleChildren %in% children], currentNode)
    
    # Move on to the next current node
    currentNode <- currentNode - 1L
  }
  
  # Make tree as list
  tree <- list(edge = edge, edge.length = edge.length, Nnode = as.integer(n - 1))
  tree$tip.label <- sample(paste0("t", c(1:n)), replace = FALSE)
  
  # Set class
  class(tree) <- "phylo"
  
  return(tree)
}

inv_cdf_coal_times_inf <- function(y, net, a, precBits) {
  one <- mpfr(1, precBits)
  rv <- mpfr(stats::runif(1, min = 0, max = 1), precBits)
  phi <-(y*a*net - rv*y*a*(net - y*a))
  return((-one / net) * log(rv/(1-rv)*1/y))
}

