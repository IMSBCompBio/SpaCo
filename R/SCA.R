#' SCA_function
#'
#' @param data gene expression data matrix; p genes as columns, n loci as rows
#' @param neighbourindexmatrix n x n matrix showing neighborhood weight of loci
#' @param PC_criterion criterion on which to select number of principal components for initial covariance matrix reconstruction; either "number" to select a number of PCs or "percent" to select number of PCs to explain specified amount of data variance
#' @param PC_value Value to specify number of PCs or desired level of explained variance, see "PC_criterion"
#' @param set_nspacs Boolean if number of relevant spacs is to be computed. Increases run time significantly
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
# Main Function
RunSCA <- function(SpaCoObject,
                   PC_criterion = "percent",
                   PC_value = 0.95,
                   set_nspacs = NULL,
                   nSim = 1000,
                   nSpacQuantile = 0.05,
                   reducedSpots = FALSE,
                   nReduce = 1000) {
  # Check for required libraries
  requiredLibraries <- c("Rcpp", "RcppEigen", "rARPACK")
  lapply(requiredLibraries, require, character.only = TRUE)

  check_RunSCA_args(as.list(environment()))
  # Extract data and neighbors
  data <- SpaCoObject@data
  neighbourIndexMatrix <- SpaCoObject@neighbours
  n <- nrow(data)
  p <- ncol(data)
  W <- sum(SpaCoObject@neighbours)

  # Adjust PC_value if needed
  adjustPCValue <- function(PC_criterion, PC_value, p) {
    if (PC_criterion == "number" && PC_value > p) {
      warning("PC_value is greater than the number of genes; using p instead.")
      return(p)
    }
    PC_value
  }

  PC_value <- adjustPCValue(PC_criterion, PC_value, p)

  if (reducedSpots)
  {
    availableSpots <- 1:nrow(data)
    ASpots <- c()
    while (length(ASpots) < nReduce)
    {
      testSample <- sample(availableSpots, 1)
      testNeighbors <-
        which(SpaCoObject@neighbours[, testSample] != 0)
      ASpots <- unique(c(ASpots, testNeighbors, testSample))
      availableSpots <- setdiff(availableSpots, ASpots)
    }
    reducedDataA <- data[ASpots,]
    coordinates <- SpaCoObject@coordinates
    reducedNeighborsA <- SpaCoObject@neighbours[ASpots, ASpots]
    reducedGraphLaplacianA <-
      computeGraphLaplacian(reducedNeighborsA)
    tmpTrainData <- reducedDataA
    tmpTrainGL <- reducedGraphLaplacianA
  } else
  {
    tmpTrainData <- data
    tmpTrainGL <- computeGraphLaplacian(SpaCoObject@neighbours)
  }

  # Data preprocessing steps
  dataCentered <- scale(tmpTrainData, scale = TRUE)



  pcaResults <- performPCA(dataCentered, PC_criterion, PC_value, n)
  dataReduced <- scale(pcaResults$dataReduced)

  # Compute test statistic matrix
  Rx <- t(dataReduced) %*% tmpTrainGL %*% dataReduced
  # if (class(neighbourIndexMatrix) == "dgCMatrix")
  # {
  #   Rx <- t(dataReduced) %*% tmpTrainGL %*% dataReduced
  # } else
  # {
  #   Rx <-
  #     eigenMapMatMult(t(dataReduced), eigenMapMatMult(tmpTrainGL, dataReduced))
  # }

  # SVD of Rx
  eigenRx <- eigen(Rx)
  PCsRx <- eigenRx$vectors[, 1:pcaResults$nEigenVals]
  lambdas <- eigenRx$values[1:pcaResults$nEigenVals]

  if (is.null(set_nspacs)) {
    nSpacs <-
      computeRelevantSpacs(nSim, 10, dataReduced, tmpTrainGL, lambdas)
    SpaCoObject@nSpacs <- nSpacs
  } else {
    nSpacs <- set_nspacs
  }

  # Compute ONB and projections
  ONBOriginalBasis <-
    pcaResults$initialPCA$vectors[, 1:pcaResults$nEigenVals] %*% PCsRx
  rownames(ONBOriginalBasis) <- colnames(dataCentered)
  colnames(ONBOriginalBasis) <-
    paste0("spac_", 1:ncol(ONBOriginalBasis))
  SpaCoObject@spacs <- ONBOriginalBasis
  tmp <-
    scale(eigenMapMatMult(data,
                          pcaResults$initialPCA$vectors[, 1:pcaResults$nEigenVals]))
  SpaCoObject@projection <- eigenMapMatMult(tmp, PCsRx)
  rownames(SpaCoObject@projection) <- rownames(data)
  colnames(SpaCoObject@projection) <-
    paste0("spac_", 1:ncol(ONBOriginalBasis))
  SpaCoObject@Lambdas <- lambdas
  SpaCoObject@GraphLaplacian <-
    computeGraphLaplacian(SpaCoObject@neighbours)
  return(SpaCoObject)
}
