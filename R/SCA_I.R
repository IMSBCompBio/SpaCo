#' SCA_function
#'
#' @param data gene expression data matrix; p genes as columns, n loci as rows
#' @param neighbourindexmatrix n x n matrix showing neighborhood weight of loci
#' @param PC_criterion criterion on which to select number of principal components for initial covariance matrix reconstruction; either "number" to select a number of PCs or "percent" to select number of PCs to explain specified amount of data variance
#' @param PC_value Value to specify number of PCs or desired level of explained variance, see "PC_criterion"
#' @param compute_nSpacs Boolean if number of relevant spacs is to be computed. Increases run time significantly
#' @param compute_projections Boolean if meta genen projections should be computed. May increase run time significantly. Default is TRUE
#' @param nSim Number of simulations for computation of spac number
#' @param nSpacQuantile Quantile to use as cutoff for spac number
#'
#' @return
#' Returns a SpaCoObject filled with the result of the spatial component analysis.
#' @export
#'
#' @import methods
#' @import rARPACK
#' @import Rcpp
#' @import RcppEigen
#'
RunSCAI <- function(SpaCoObject, PC_criterion = "percent",
                    PC_value = .8, compute_nSpacs = FALSE,
                    compute_projections = TRUE, nSim = 1000, nSpacQuantile = 0.05)
{
  require(Rcpp)
  require(RcppEigen)
  require(rARPACK)
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
  GraphLaplacian <- neighbourindexmatrix * 1 / W
  #Center data
  data_centered <- scale(data, scale = TRUE)
  #Scale data using spatial scalar product
  GeneANorms <- colSums(data_centered *
                          eigenMapMatMult(GraphLaplacian, data_centered))
  data_centered_GL_scaled <- sweep(data_centered, 2, GeneANorms + 1, "*")
  #Perform initial PCA for dimension reduction
  VarMatrix <- (1 / (n - 1)) * eigenMapMatMult(t(data_centered_GL_scaled), data_centered_GL_scaled)
  InitialPCA <- eigen(VarMatrix, symmetric = TRUE)
  if(PC_criterion == "percent")
  {
    if(PC_value == 1)
    {
      nEigenVals <- p
    }else
    {
      nEigenVals <- min(which(cumsum(InitialPCA$values)/sum(InitialPCA$values) > PC_value))
    }
  }else
  {
    nEigenVals <- PC_value
  }
  data_reduced <- eigenMapMatMult(data_centered_GL_scaled, InitialPCA$vectors[,1:nEigenVals])
  # data_reduced <- eigenMapMatMult(data_reduced, t(InitialPCA$vectors[,1:nEigenVals]))
  data_reduced <- scale(data_reduced)

  #Compute test statistic matrix
  R_x <- eigenMapMatMult(t(data_reduced), eigenMapMatMult(GraphLaplacian, data_reduced))

  #Compute SVD of R_x
  Eigen_Rx <- eigen(R_x)
  #Reverse order
  PCs_Rx <- Eigen_Rx$vectors[,1:nEigenVals]
  Lambdas <- Eigen_Rx$values[1:nEigenVals]
  #Reconstruct selected principal components into original basis before
  #SVD of test statistic matrix

  if(compute_nSpacs)
  {
    message("computing number of releveant spacs")
    # GLEigen <- eigen(GraphLaplacian)
    # # nonZeroIdx <- which(abs(GLEigen$values) > 1e-8)
    # GLEigenVals <- GLEigen$values#[1:k]
    # L <- eigenMapMatMult(t(GLEigen$vectors), diag(sign(GLEigenVals) * sqrt(1/abs(GLEigenVals))))
    simSpacCFunction <- function(i)
    {
      shuffleOrder <- sample(ncol(GraphLaplacian), ncol(GraphLaplacian))
      # tmpMat <- eigenMapMatMult(t(data_reduced[shuffleOrder,]), L)
      R_x_Shuffled <- eigenMapMatMult(t(data_reduced[shuffleOrder,]),
                                      eigenMapMatMult(GraphLaplacian,
                                                      data_reduced[shuffleOrder,]))
      # R_x_Shuffled <- eigenMapMatMult(tmpMat, t(tmpMat)) * preFactor
      largestEigVal <- eigs_sym(R_x_Shuffled, 1, which = "LM")$values
      return(largestEigVal)
    }
    #Sample nSim minimal projection eigenvalues
    results_all <- t(sapply(1:nSim, simSpacCFunction))

    nSpacs <- min(which(Lambdas < quantile(results_all, 1 - nSpacQuantile)))
    slot(SpaCoObject, "nSpacs") <- nSpacs
  }

  ONB_OriginalBasis <- InitialPCA$vectors[,1:nEigenVals] %*%
    PCs_Rx
  rownames(ONB_OriginalBasis) <- colnames(data_centered)
  colnames(ONB_OriginalBasis) <- paste0("spac_",1:ncol(ONB_OriginalBasis))

  slot(SpaCoObject, "spacs") <- ONB_OriginalBasis
  if (compute_projections) {
    message("computing projections this may take a while")
    slot(SpaCoObject, "projection") <- eigenMapMatMult(data_reduced, PCs_Rx)
    rownames(SpaCoObject@projection) <- rownames(data_centered)
    colnames(SpaCoObject@projection) <- paste0("spac_",1:ncol(ONB_OriginalBasis))
  }
  slot(SpaCoObject, "Lambdas") <- Lambdas
  slot(SpaCoObject,"GraphLaplacian") <- GraphLaplacian
  return(SpaCoObject)
}
