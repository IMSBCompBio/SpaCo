#' Title
#'
#' @param x
#' @param neighbourindexmatrix
#'
#' @return
#' @export
#'
#' @examples
compute_C <- function(x, neighbourindexmatrix)
{
  n <- nrow(neighbourindexmatrix)
  W <- sum(neighbourindexmatrix)
  activeIndices <- which(neighbourindexmatrix != 0)
  i <- activeIndices %% n
  i[which(i == 0)] <- n
  j <- ceiling(activeIndices / n)
  C <- (((n - 1) / (2 * n * W)) *
          sum(neighbourindexmatrix[activeIndices] * (x[i] - x[j])^2)) / stats::var(x)
  return(C)
}
#' Title
#'
#' @param X
#' @param type
#' @param W_Matrix
#'
#' @return
#' @export
#'
#' @examples
compute_A <- function(X, type = "C", W_Matrix) {
  stopifnot(type %in% c("C", "I", "1/V"))
  X_length <- ncol(X)
  A_x <- matrix(0, nrow = nrow(X), ncol = nrow(X))
  W <- sum(W_Matrix)
  if(type == "C") {
    res <- Map(function(i) {
      d_i <- sum(W_Matrix[i,])
      neighborLoci <- which(W_Matrix[i,] != 0)
      x_i_bar <- 1 / d_i * rowSums(X[,neighborLoci] *
                                     matrix(W_Matrix[i,neighborLoci], ncol = length(neighborLoci),
                                            nrow = nrow(X), byrow = TRUE))
      2 * d_i * X[,i] %*% t(X[,i] - x_i_bar)
    }, 1:X_length)
    A_x <- Reduce("+", res)
    A_x <- A_x * ((X_length - 1) / (2 * X_length * W))
  }

  return(A_x)
}
#' Title
#'
#' @param x
#' @param A
#'
#' @return
#' @export
#'
#' @examples
normalizeA <- function(x, A)
{
  return((x / c(normA(x, A))))
}
#' Title
#'
#' @param x
#' @param A
#'
#' @return
#' @export
#'
#' @examples
normA <- function(x, A)
{
  return(sqrt(scalarProductA(x, x, A)))
}
scalarProductA <- function(x, y, A)
{
  n <- nrow(A)
  W <- -(1/2)*sum(A - diag(diag(A)))
  preFactor <- (n-1)/(2*n*W)
  return(preFactor * t(x) %*% A %*% y)
}
#' Title
#'
#' @param v
#' @param u
#' @param A
#'
#' @return
#' @export
#'
#' @examples
projAFunction <- function(v, u, A)
{
  return(c(scalarProductA(v, u, A)/scalarProductA(u, u, A)) * u)
}
#' Title
#'
#' @param v
#' @param ONB
#' @param A
#'
#' @return
#' @export
#'
#' @examples
projASubspaceFunction <- function(v, projMatrix, preFactor)
{
  projection <- preFactor * projMatrix %*% v
  return(projection)
}
#' Title
#'
#' @param X
#' @param A
#' @param tol
#'
#' @return
#' @export
#'
#' @examples
orthogonalizeA <- function(X, A, tol = .Machine$double.eps^0.5)
{
  m <- nrow(X)
  n <- ncol(X)
  if (m < n)
    stop("No. of rows of 'A' must be greater or equal no. of colums.")
  Q <- matrix(0, m, n)
  Norms <- rep(0, n)
  for (k in 1:n) {
    Q[, k] <- X[, k]
    if (k > 1) {
      for (i in 1:(k - 1)) {
        Q[, k] <- Q[, k] - projAFunction(Q[, k], Q[, i], A)
      }
    }
    Norms[k]  <- normA(Q[, k], A)
    if (abs(Norms[k]) <= tol)
      stop("Matrix 'A' does not have full rank.")
    Q[, k] <- Q[, k] / c(Norms[k])
  }
  return(list(Q = Q))
}
#' compute_nSpacs
#'
#' @param SpacoObject
#' @param nSim number of simulations i.e. random permutation matrices
#' @param PC_criterion criterion on which to select number of principal components for initial covariance matrix reconstruction; either "number" to select a number of PCs or "percent" to select number of PCs to explain specified amount of data variance
#' @param PC_value Value to specify number of PCs or desired level of explained variance, see "PC_criterion"
#' @param SpacQuantile quantile of simulated Spac distribution used for determining cutoff
#'
#' @return
#' @export
#'
#' @examples
compute_nSpacs <- function(SpacoObject, nSim, PC_criterion = "percent",
                           PC_value = .9, SpacQuantile = 0.5)
{
  require(parallel)
  data = SpacoObject@data
  neighbours = SpacoObject@neighbours
  GraphLaplacian = SpacoObject@GraphLaplacian
  Lambdas = SpacoObject@Lambdas
  GraphLaplacian <- SpacoObject@GraphLaplacian
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
    return(Svd_Rx_Shuffled$d[length(Svd_Rx_Shuffled$d)])
  }
  numcores = detectCores()
  cl = makeCluster(numcores - 3)
  clusterExport(cl,
                list("data_reduced", "GraphLaplacian",
                     "simSpacCFunction", "preFactor"),
                envir = environment())
  #Apply Gene Score Function to all genes
  results_all <- t(parSapply(cl, 1:nSim, simSpacCFunction))
  stopCluster(cl)

  nSpacs <- min(which(Lambdas > quantile(results_all, SpacQuantile)))
  return(nSpacs)
}


#' Title
#'
#' @param SpaCoObject Spaco object to compute profiles of
#' @param nSpacs number of spac's to use to compute the smoothed profiles. Should be determined by the compute_nSpacs function.
#'
#' @return
#' @export
#'
#' @examples
compute_smoothed_profiles <- function(SpaCoObject, nSpacs) {
smoothed <- SpaCoObject@spacs %*% t(SpaCoObject@projection)
slot(SpaCoObject,"smoothed") <- as.data.frame(t(smoothed))
return(SpaCoObject)
}


