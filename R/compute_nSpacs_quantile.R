compute_nSpacs_quantile <- function(SpaCoObject, nSim, PC_criterion = "percent",
                           PC_value = .9, SpacQuantile = 0.05)
{
  require(parallel)
  data = SpaCoObject@data
  neighbours = SpaCoObject@neighbours
  GraphLaplacian = SpaCoObject@GraphLaplacian
  Lambdas = SpaCoObject@Lambdas
  GraphLaplacian <- SpaCoObject@GraphLaplacian
  n <- nrow(data)
  p <- ncol(data)
  W <- sum(neighbours)
  preFactor <- (n - 1)/(2 * n * W)
  #Input check: Check if number of desired number is larger than number of genes
  if(PC_criterion == "number")
  {
    if(PC_value > p)
    {
      PC_value <- p
      warning("Desired number of principal components is larger than number of genes. Using number of genes instead.")
    }
  }
  #Center data
  data_centered <- scale(data, scale = FALSE)
  #Scale data using spatial scalar product
  GeneANorms <- sqrt(preFactor * colSums(data_centered * (GraphLaplacian %*% data_centered)))
  data_centered_GL_scaled <- sweep(data_centered, 2, GeneANorms, "/")
  #Perform initial PCA for dimension reduction
  VarMatrix <- (1 / (n - 1)) * t(data_centered_GL_scaled) %*% data_centered_GL_scaled
  InitialPCA <- svd(VarMatrix)
  if(PC_criterion == "percent")
  {
    if(PC_value == 1)
    {
      nEigenVals <- p
    }else
    {
      nEigenVals <- min(which(cumsum(InitialPCA$d)/sum(InitialPCA$d) > PC_value))
    }
  }else
  {
    nEigenVals <- PC_value
  }
  data_reduced <- t(data_centered_GL_scaled %*% InitialPCA$v[,1:nEigenVals])
  data_reduced <- t(scale(t(data_reduced)))
  simSpacCFunction <- function(i)
  {
    shuffleOrder <- sample(ncol(GraphLaplacian), ncol(GraphLaplacian))
    permutationMatrix <- diag(length(shuffleOrder))[shuffleOrder,]
    R_x_Shuffled <- preFactor * data_reduced %*% t(permutationMatrix) %*%
      GraphLaplacian %*% permutationMatrix %*% t(data_reduced)
    Svd_Rx_Shuffled <- svd(R_x_Shuffled)
    K <- which.max(outer(Svd_Rx_Shuffled$d[length(Svd_Rx_Shuffled$d):1], Lambdas, "<"))-1
    if (K<=0) {
      K <- Inf
    }
    return(K)
  }
  numcores = detectCores()
  cl = makeCluster(numcores - 3)
  clusterExport(cl,
                list("data_reduced", "GraphLaplacian",
                     "simSpacCFunction", "preFactor","Lambdas"),
                envir = environment())
  #Apply Gene Score Function to all genes
  results_all <- t(parSapply(cl, 1:nSim, simSpacCFunction))
  stopCluster(cl)

  nSpacs <- quantile(results_all, SpacQuantile)
  return(nSpacs)
}
