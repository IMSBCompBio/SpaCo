#' SCA_function
#'
#' @param SpaCoObject gene expression data matrix; p genes as columns, n loci as rows
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
RunSCA2 <- function(SpaCoObject, PC_criterion = "percent",
                   PC_value = .8, compute_nSpacs = FALSE,
                   compute_projections = TRUE, nSpacQuantile = 0.05)
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
  #Compute Graph Laplacian
  GraphLaplacian <- -neighbourindexmatrix
  diag(GraphLaplacian) <- rowSums(neighbourindexmatrix)
  GraphLaplacian <- 2 * GraphLaplacian
  #Center data
  data_centered <- scale(data, scale = FALSE)
  #Scale data using spatial scalar product
  GeneANorms <- sqrt((n - 1)/(2 * n * W) * colSums(data_centered *
                                                     eigenMapMatMult(GraphLaplacian, data_centered)))
  data_centered_GL_scaled <- sweep(data_centered, 2, GeneANorms, "/")
  #Perform initial PCA for dimension reduction
  VarMatrix <- (1 / (n - 1)) * eigenMapMatMult(t(data_centered_GL_scaled), data_centered_GL_scaled)
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
  data_reduced <- t(eigenMapMatMult(data_centered_GL_scaled, InitialPCA$v[,1:nEigenVals]))
  data_reduced <- t(scale(t(data_reduced)))
  #Compute test statistic matrix
  R_x <- preFactor * eigenMapMatMult(data_reduced, eigenMapMatMult(GraphLaplacian, t(data_reduced)))
  #Compute SVD of R_x
  Svd_Rx <- svd(R_x)
  #Reverse order
  PCs_Rx <- Svd_Rx$u[,nEigenVals:1]
  Lambdas <- Svd_Rx$d[nEigenVals:1]
  #Reconstruct selected principal components into original basis before
  #SVD of test statistic matrix
  if(compute_nSpacs)
  {
    message("computing number of releveant spacs")
    GLeigen <- eigen(GraphLaplacian)
    CVar <- sum(GLeigen$values^2)
    pVals <- pnorm(Lambdas, mean = 1, sd = sqrt(CVar/(n^2)))
    if (suppressWarnings(max(which(pVals < nSpacQuantile))) == -Inf) {
      warning("There are no significant spac's in the level of ", nSpacQuantile )
      slot(SpaCoObject, "nSpacs") <- as.integer(0)
      } else {
    message("Using quantile level of " ,nSpacQuantile)
    nSpacs <- max(which(pVals < nSpacQuantile))
    slot(SpaCoObject, "nSpacs") <- nSpacs
      }
    }
  ONB_OriginalBasis <- t(t(PCs_Rx) %*% t(InitialPCA$v[,1:nEigenVals]))
  rownames(ONB_OriginalBasis) <- colnames(data_centered)
  colnames(ONB_OriginalBasis) <- paste0("spac_",1:ncol(ONB_OriginalBasis))
  slot(SpaCoObject, "spacs") <- ONB_OriginalBasis
  if (compute_projections) {
    message("computing projections this may take a while")
    slot(SpaCoObject, "projection") <- t(eigenMapMatMult(t(ONB_OriginalBasis), t(data)))
    rownames(SpaCoObject@projection) <- rownames(data_centered)
    colnames(SpaCoObject@projection) <- paste0("spac_",1:ncol(ONB_OriginalBasis))
    var_projections <- apply(SpaCoObject@projection,2,sd)
    slot(SpaCoObject, "projection") <- sweep(SpaCoObject@projection, MARGIN = 2, STATS = var_projections, FUN = "/")
    slot(SpaCoObject, "spacs") <- sweep(SpaCoObject@spacs, MARGIN = 2, STATS = var_projections, FUN = "/")
  }
  slot(SpaCoObject, "Lambdas") <- Lambdas
  slot(SpaCoObject,"GraphLaplacian") <- GraphLaplacian
  return(SpaCoObject)
}
