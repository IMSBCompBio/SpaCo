#' SCA_function
#'
#' @param data gene expression data matrix; p genes as columns, n loci as rows
#' @param neighbourindexmatrix n x n matrix showing neighborhood weight of loci
#' @param PC_criterion criterion on which to select number of principal components for initial covariance matrix reconstruction; either "number" to select a number of PCs or "percent" to select number of PCs to explain specified amount of data variance
#' @param PC_value Value to specify number of PCs or desired level of explained variance, see "PC_criterion"
#' @param orthogonalResult Logical value to specify if Spacos should be orthogonalized to form a ONB; since transformation of eigenvalues results in non-orthogonal Spacos
#' @compute_projections Boolean if meta genen projections should be computed. May increas runtime significantly.
#' @return
#' @export
#'
#' @examples
#' @import pracma
#' @import MASS
#' @import dplyr
#' @import ggplot2
#'
SCA_function <- function(SpaCoObject, PC_criterion = "percent",
                         PC_value = .8, orthogonalResult = FALSE, compute_projections = F)
{
  require(pracma)
  require(MASS)
  require(dplyr)
  require(ggplot2)
   if(!PC_criterion %in% c("percent", "number"))
  {
    stop("PC_criterion must be either \"percent\" or \"number\".")
  }
  if(PC_criterion == "percent" & (PC_value <= 0 | PC_value > 1))
  {
    stop("Desired level of PC explained variance must be between 0 and 1.")
  }
  if(PC_criterion == "number" & (PC_value <= 0 | PC_value%%1 != 0))
  {
    stop("Desired number of PCs must be a positive integer.")
  }

  #n = Number of Loci; p = number of genes
  data <- SpaCoObject@data
  neighbourindexmatrix <-SpaCoObject@neighbours
  n <- nrow(data)
  p <- ncol(data)
  W <- sum(neighbourindexmatrix)
  #Input check: Check if number of desired number is larger than number of genes
  if(PC_criterion == "number")
  {
    if(PC_value > p)
    {
      PC_value <- p
      warning("Desired number of principal components is larger than number of genes. Using number of genes instead.")
    }
  }
  #Compute Graph Laplacian
  GraphLaplacian <- -neighbourindexmatrix
  diag(GraphLaplacian) <- rowSums(neighbourindexmatrix)
  GraphLaplacian <- 2 * GraphLaplacian

  #Center data
  data_centered <- scale(data, scale = FALSE)
  #Scale data using spatial scalar product
  data_centered_GL_scaled <- t(apply(data_centered, 2, normalizeA, A = GraphLaplacian))
  #Perform initial PCA for dimension reduction
  VarMatrix <- (1 / (n - 1)) * data_centered_GL_scaled %*% t(data_centered_GL_scaled)
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
  data_reduced <- t(t(data_centered_GL_scaled) %*% InitialPCA$v[,1:nEigenVals])
  #Compute Covariance matrix
  Bx <- var(t(data_reduced))
  # #Computing Geary's C for  principal components
  # PCA_C <- rep(0,nEigenVals)
  # for(k in 1:nEigenVals)
  # {
  #   PCA_C[k] <- compute_C(data_reduced[k,], neighbourindexmatrix)
  # }
  # PCA_C <- PCA_C[order(PCA_C)]
  EigenVals <- diag(Bx)
  L <- diag(sqrt(EigenVals))
  L_Minus <- diag(1/sqrt(EigenVals))

  #Compute numerator of Rayleigh Quotient
  A_x <- (n-1)/(2*n*W) * data_reduced %*% GraphLaplacian %*% t(data_reduced)

  #Compute test statistic matrix
  R_x <- L_Minus %*% A_x %*% t(L_Minus)

  #Compute SVD of R_x
  Svd_Rx <- svd(R_x)
  #Reverse order
  PCs_Rx <- Svd_Rx$u[,nEigenVals:1]
  Lambdas <- Svd_Rx$d[nEigenVals:1]
  #Reconstruct selected principal components into original basis before
  #SVD of test statistic matrix
  PCs_Rx_OriginalBasis <- L_Minus %*% PCs_Rx

  #Create orthogonal basis from transformed PCs
  if(orthogonalResult)
  {
    ONB <- gramSchmidt(PCs_Rx_OriginalBasis)$Q
  }else
  {
    ONB <- PCs_Rx_OriginalBasis
  }

  ONB_OriginalBasis <- t(t(ONB) %*% t(InitialPCA$v[,1:nEigenVals]))
  rownames(ONB_OriginalBasis) <- colnames(data_centered)
  colnames(ONB_OriginalBasis) <- paste0("spac_",1:ncol(ONB_OriginalBasis))
  # Spaco_C <- Svd_Rx$d[nEigenVals:1]
  #
  # nSpacos <- min(which(Spaco_C > PCA_C))

  outputList <- list(ONB_OriginalBasis, Lambdas,
                     #nSpacos,
                     R_x)
  #return(outputList)
  slot(SpaCoObject, "spacs") <- ONB_OriginalBasis
  print("computing projections this may take a while")
  if (compute_projections) {
    message("computing projections this may take a while")
    slot(SpaCoObject, "projection") <- t(t(ONB_OriginalBasis) %*% t(data))
  }
  slot(SpaCoObject, "Lambdas") <- Lambdas
  #slot(SpaCoObject, "R_x") <- R_x
  slot(SpaCoObject,"GraphLaplacian") <- GraphLaplacian
  return(SpaCoObject)
}
