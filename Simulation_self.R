
calculate_leaf_info <- function(tree_data, n) {
  # Calculate the leaves
  leaves <- setdiff(tree_data$Node, tree_data$Parent)
  
  get_leaves <- function(node, tree) {
    children <- tree$Node[tree$Parent == node]
    if (length(children) == 0) {
      return(node)
    } else {
      return(unlist(lapply(children, get_leaves, tree = tree)))
    }
  }
  
  results <- data.frame(Edge = integer(), Leaves = integer(), Total_Length = numeric())
  # Calculate the lenght of branches that support k = 2,.., n-1 leaves
  for (i in 1:nrow(tree_data)) {
    node <- tree_data$Node[i]
    edge_length <- tree_data$Edge_length[i]
    
    descendant_leaves <- get_leaves(node, tree_data)
    
    num_leaves <- length(descendant_leaves)
    total_length <- edge_length
    
    results <- rbind(
      results,
      data.frame(Edge = node, Leaves = num_leaves, Total_Length = total_length)
    )
  }
  
  summed_results <- aggregate(Total_Length ~ Leaves, data = results, FUN = sum)
  
  check_leaves <- data.frame(Leaves = 1:(n - 1), Exists = FALSE, Total_Length_Positive = FALSE)
  
  for (k in check_leaves$Leaves) {
    # Check wheter k occurs
    exists <- any(results$Leaves == k)
    positive_length <- any(results$Leaves == k & results$Total_Length > 0)
    
    
    check_leaves$Exists[check_leaves$Leaves == k] <- exists
    check_leaves$Total_Length_Positive[check_leaves$Leaves == k] <- positive_length
  }
  return(list(results = results, summed_results = summed_results, check_leaves = check_leaves))
}


inputCheck <- function(tree, alpha) {
  # Must be of class phylo
  if (!inherits(tree, "phylo")) {
    stop("Tree must be of class phylo. Use as.phylo function to convert if the
    formatting is correct. Otherwise, see ape package documentation
    https://cran.r-project.org/web/packages/ape/ape.pdf")
  }
  
  # Only works for ultrametric trees
  if (!ape::is.ultrametric(tree)) {
    stop("Tree is not ultrametric. internalLengths, and maxLike fns.
        should only be used with ultrametric trees. For usage with mutation trees,
        use sharedMuts fn.")
  }
  
  # Make sure alpha is reasonable
  if (alpha < 0 | alpha > 1) {
    stop("alpha must be between 0 and 1")
  }
  if (alpha > 0.25) {
    warning(paste0("We calulate 1-alpha confidence intervals. The given confidence
            intervals, with alpha = ", alpha, " correspond to ", 1 - alpha, "%
            confidence intervals, which will be very narrow."))
  }
  return(NULL)
}


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

inv_cdf_coal_times <- function(y, net, a, alpha, precBits) {
  one <- mpfr(1, precBits)
  rv <- mpfr(stats::runif(1, min = 0, max = 1), precBits)
  phi <- (alpha * rv) / (a * (one - alpha * (one - y)))
  return((-one / net) * log((one - y * a * phi) / (one + (net - y * a) * phi)))
}


simUltra <- function(a, b, cloneAge, n, nTrees = 1,
                     precBits = 1000, addStem = FALSE, nCores = 1) {
  # Store runtime for each tree
  ptm <- proc.time()
  
  # Set nTrees equal to "SKIP_TESTS" on recursive runs to skip repeated checks
  if (nTrees == "SKIP_TESTS") {
    nTrees <- 1
  } else {
    # Make sure length of params is either one or equal to 'nTrees'
    if (!all(unlist(lapply(list(a, b, cloneAge, n), length)) == 1 |
             unlist(lapply(list(a, b, cloneAge, n), length)) == nTrees)) {
      stop(paste0("Input parameters must be length 1 or length equal to the value
                of param 'nTrees', which is ", nTrees))
    }
    
    
    # If a is length 1, make it a vector of length nTrees for first argument to mapply
    if (length(a) == 1) {
      a <- rep(a, nTrees)
    } else {
      a <- a
    }
  }
  
  # Call recursively to generate nTrees if nTrees > 1
  if (nTrees > 1) {
    # Parallelize if "parallel" pkg avail. (user must set nCores > 1 explicitly)
    if (requireNamespace("parallel", quietly = TRUE)) {
      return.list <- parallel::mcmapply(simUltra,
                                        a = a,
                                        b = b,
                                        cloneAge = cloneAge, n = n, precBits = precBits,
                                        addStem = addStem, nTrees = "SKIP_TESTS", nCores = 1,
                                        mc.cores = nCores, SIMPLIFY = FALSE
      )
    } else {
      return.list <- mapply(
        FUN = simUltra, a = a,
        b = b, cloneAge = cloneAge,
        n = n, precBits = precBits, addStem = addStem,
        nTrees = "SKIP_TESTS", SIMPLIFY = FALSE
      )
    }
    return(return.list)
  }
  
  # Convert params to high precision mpfr
  cloneAge_mpfr <- mpfr(cloneAge, precBits)
  a_mpfr <- mpfr(a, precBits)
  b_mpfr <- mpfr(b, precBits)
  net_mpfr <- a_mpfr - b_mpfr
  n_mpfr <- mpfr(n, precBits)
  one <- mpfr(1, precBits)
  
  # Define alpha and ensure it doesn't map to 1 at given precision
  alpha_mpfr <- (a_mpfr * (exp(net_mpfr * cloneAge_mpfr) - one)) / (a_mpfr * exp(net_mpfr * cloneAge_mpfr) - b_mpfr)
  if (alpha_mpfr == one) {
    stop("alpha value is equal to 1 due to insufficient machine precision. This
          will lead to NaN coalescence times. Increase param 'precBits'.")
  }
  
  # Draw Y = y from the inverse CDF
  uniform_rv <- mpfr(stats::runif(1, min = 0, max = 1), precBits)
  y_mpfr <- ((one - alpha_mpfr) * (uniform_rv**(one / n_mpfr))) / (one - alpha_mpfr * uniform_rv**(one / n_mpfr))
  
  # Generate the coalescence times
  coal_times_mpfr <- sapply(rep(y_mpfr, n - 1), inv_cdf_coal_times,
                            net = net_mpfr,
                            a = a_mpfr, alpha = alpha_mpfr, precBits = precBits
  )
  
  # Convert back to normal numeric (no longer need high precision)
  coal_times <- suppressWarnings(sapply(coal_times_mpfr, Rmpfr::asNumeric))
  
  # Convert coal times into tree by randomly merging lineages
  tree <- coal_to_tree(coal_times)
  
  # Add stem starting the tree from zero, rooting the tree appropriately
  if (addStem) {
    tree$edge[tree$edge > n] <- tree$edge[tree$edge > n] + 1
    tree$edge <- rbind(c(n + 1, n + 2), tree$edge)
    tree$edge.length <- c(cloneAge - max(coal_times), tree$edge.length)
    tree$Nnode <- tree$Nnode + 1
  }
  
  # Sanity checks (coalescence times must match and be less than cloneAge)
  stopifnot(all(coal_times <= cloneAge))
  if (!all(round(coal_times, 4) %in% round(ape::branching.times(tree), 4))) {
    stop("Unexpected error: coalescence times not matching between drawn times
         and output tree using ape::branching.times()")
  }
  
  # Add metadata for making the tree
  runtime <- proc.time()[["elapsed"]] - ptm[["elapsed"]]
  tree$metadata <- data.frame(
    "r" = a - b, "a" = a, "b" = b, "cloneAge" = cloneAge,
    "n" = n, "runtime_seconds" = runtime, "addStem" = addStem
  )
  
  # Return the tree created from the coalescence times drawn from Lambert distribution
  return(list(tree, coal_times))
}


internalLengths <- function(tree, alpha = 0.05) {
  # time function
  ptm <- proc.time()
  
  # If we have a list of phylo objects instead of a single phylo objects, call recursively
  if (inherits(tree, "list") & !inherits(tree, "phylo")) {
    # Call function recursively on all trees in list, then combine results into one data.frame
    return.df <- do.call(rbind, lapply(tree, internalLengths, alpha = alpha))
    return.df$cloneName_result <- names(tree)
    return(return.df)
  }
  # Get number of tips
  n <- ape::Ntip(tree)
  # Check if tree has stem
  nodes <- tree$edge[tree$edge > n]
  if (1 %in% table(nodes)) {
    hasStem <- TRUE
    stemNode <- as.numeric(names(which(table(nodes) == 1)))
  } else {
    hasStem <- FALSE
  }
  # Get the number of direct descendants from a node, identifying the nodes with > 2
  countChildren <- table(tree$edge[, 1])
  # Get list of descendants from each internal node
  descendant_df <- data.frame(
    "Node" = (n + 2):max(tree$edge), "Parent" = NA,
    "Edge_length" = NA
  )
  # Find parent and edge length preceding each internal node
  for (k in descendant_df$Node) {
    descendant_df$Edge_length[descendant_df$Node == k] <- tree$edge.length[which(tree$edge[, 2] == k)]
    descendant_df$Parent[descendant_df$Node == k] <- tree$edge[which(tree$edge[, 2] == k), 1]
  }
  # If tree has a stem, remove stem from calculation
  if (hasStem) {
    descendant_df <- descendant_df[!descendant_df$Parent == stemNode, ]
  }
  # The sum of edge lengths in descendant_df is equal to the total internal lengths
  intLen <- sum(descendant_df$Edge_length)
  # Calculate growth rate and confidence intervals
  growthRate <- n / intLen
  growthRate_lb <- growthRate * (1 + stats::qnorm(alpha / 2) / sqrt(n))
  growthRate_ub <- growthRate * (1 - stats::qnorm(alpha / 2) / sqrt(n))
  # Calculate total external lengths
  extLen <- sum(tree$edge.length[tree$edge[, 2] %in% c(1:n)])
  
  # Estimate clone age. If tree has stem, take tree age, otherwise estimate by adding 1/r
  if (hasStem) {
    cloneAgeEstimate <- max(ape::branching.times(tree))
  } else {
    cloneAgeEstimate <- max(ape::branching.times(tree)) + 1 / growthRate
  }
  # Get runtime (including all tests)
  runtime <- proc.time() - ptm
  
  result.df <- data.frame(
    "lowerBound" = growthRate_lb, "estimate" = growthRate,
    "upperBound" = growthRate_ub, "cloneAgeEstimate" = cloneAgeEstimate,
    "sumInternalLengths" = intLen,
    "sumExternalLengths" = extLen, extIntRatio = extLen / intLen,
    "n" = n, "alpha" = alpha, "runtime_s" = runtime[["elapsed"]],
    "method" = "lengths"
  )
  # added descendant_df
  return(list(result.df, descendant_df))
}

n = 6
test =simUltra(3,1,100,n)
res_new = internalLengths(test[[1]])
input = calculate_leaf_info(res_new[[2]],n)
print(res_new[[2]])
res = ew_durrett_adapted(n, res_new[[2]])
print(res)
rep_new <- function(numRep, n, t, r){
  est = rep(0, numRep)
  for(i in 1:numRep){
    a = runif(1, min = r, max = r+1)
    
    test =simUltra(a,a-r,t,n)
    res_new = internalLengths(test[[1]])
    input = calculate_leaf_info(res_new[[2]],n)
    #print(res_new[[2]])
    res = ew_durrett_adapted(n, res_new[[2]])
    est[i] = res
  }
  return(mean(est))
  
}
print(rep_new(100, 9, 40, 2))

