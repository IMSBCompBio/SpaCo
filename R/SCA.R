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
# Main Function
RunSCA <- function(SpaCoObject,
                   PC_criterion = "percent",
                   PC_value = 0.95,
                   compute_nSpacs = FALSE,
                   nSim = 1000,
                   nSpacQuantile = 0.05,
                   reducedSpots = FALSE,
                   nReduce = 1000) {
  # Check for required libraries
  requiredLibraries <- c("Rcpp", "RcppEigen", "rARPACK")
  lapply(requiredLibraries, require, character.only = TRUE)

  # Validate input parameters
  validateInputs <- function(PC_criterion, PC_value) {
    if (!PC_criterion %in% c("percent", "number")) {
      stop("PC_criterion must be either 'percent' or 'number'.")
    }
    if (PC_criterion == "percent" && (PC_value <= 0 || PC_value > 1)) {
      stop("PC_value must be between 0 and 1 for percent criterion.")
    }
    if (PC_criterion == "number" && (PC_value <= 0 || PC_value %% 1 != 0)) {
      stop("PC_value must be a positive integer for number criterion.")
    }
  }
  validateInputs(PC_criterion, PC_value)

  # Extract data and neighbors
  data <- SpaCoObject@data
  neighbourIndexMatrix <- SpaCoObject@neighbours
  n <- nrow(data)
  p <- ncol(data)
  W <- sum(neighbourIndexMatrix)

  # Adjust PC_value if needed
  adjustPCValue <- function(PC_criterion, PC_value, p) {
    if (PC_criterion == "number" && PC_value > p) {
      warning("PC_value is greater than the number of genes; using p instead.")
      return(p)
    }
    PC_value
  }

  PC_value <- adjustPCValue(PC_criterion, PC_value, p)

  # Compute Graph Laplacian
  computeGraphLaplacian <- function(neighbourIndexMatrix) {
    W <- sum(neighbourIndexMatrix)
    n <- nrow(neighbourIndexMatrix)
    graphLaplacian <- as.matrix(neighbourIndexMatrix / W)
    graphLaplacian <- graphLaplacian + diag(1 / n, nrow = n)
  }

  if(reducedSpots)
  {
    availableSpots <- 1:nrow(data)
    ASpots <- c()
    while(length(ASpots) < nReduce)
    {
      testSample <- sample(availableSpots, 1)
      testNeighbors <- which(neighbourIndexMatrix[,testSample] != 0)
      ASpots <- unique(c(ASpots, testNeighbors, testSample))
      availableSpots <- setdiff(availableSpots, ASpots)
    }
    reducedDataA <- data[ASpots,]
    coordinates <- SpaCoObject@coordinates
    reducedNeighborsA <- neighbourIndexMatrix[ASpots, ASpots]
    reducedGraphLaplacianA <- computeGraphLaplacian(reducedNeighborsA)
    tmpTrainData <- reducedDataA
    tmpTrainGL <- reducedGraphLaplacianA
  }else
  {
    tmpTrainData <- data
    tmpTrainGL <- computeGraphLaplacian(neighbourIndexMatrix)
  }

  # Data preprocessing steps
  dataCentered <- scale(tmpTrainData, scale = TRUE)

  # PCA and dimension reduction
  performPCA <- function(data, criterion, value) {
    varMatrix <- (1 / (n - 1)) * eigenMapMatMult(t(data), data)
    initialPCA <- eigen(varMatrix, symmetric = TRUE)
    nEigenVals <- if (criterion == "percent") {
      if (value == 1) p else min(which(cumsum(initialPCA$values) / sum(initialPCA$values) > value))
    } else {
      value
    }
    list(dataReduced = t(eigenMapMatMult(diag(1/sqrt(initialPCA$values[1:nEigenVals])), eigenMapMatMult(t(initialPCA$vectors[, 1:nEigenVals]), t(data)))),
         nEigenVals = nEigenVals,
         initialPCA = initialPCA)
  }

  pcaResults <- performPCA(dataCentered, PC_criterion, PC_value)
  dataReduced <- scale(pcaResults$dataReduced)

  # Compute test statistic matrix
  Rx <- eigenMapMatMult(t(dataReduced), eigenMapMatMult(tmpTrainGL, dataReduced))

  # SVD of Rx
  eigenRx <- eigen(Rx)
  PCsRx <- eigenRx$vectors[, 1:pcaResults$nEigenVals]
  lambdas <- eigenRx$values[1:pcaResults$nEigenVals]

  # Compute number of relevant SPACs if required
  computeRelevantSpacs <- function(nSim, batchSize, dataReduced, graphLaplacian, lambdas) {
    simSpacFunction <- function(i) {
      shuffleOrder <- sample(ncol(graphLaplacian), ncol(graphLaplacian))
      RxShuffled <- eigenMapMatMult(t(dataReduced[shuffleOrder, ]),
                                    eigenMapMatMult(graphLaplacian, dataReduced[shuffleOrder, ]))
      eigs_sym(RxShuffled, 1, which = "LM")$values
    }
    batchSize <- 10
    resultsAll <- replicate(100, simSpacFunction())
    eigValSE <- sd(resultsAll) / sqrt(length(resultsAll))
    eigValCI <- mean(resultsAll) + qt(0.975, df = length(resultsAll) - 1) * eigValSE * c(-1, 1)
    lambdasInCI <- lambdas[lambdas > eigValCI[1] & lambdas < eigValCI[2]]
    if(length(lambdasInCI) > 1)
    {
      for(i in 1:round((nSim - 100) / batchSize))
      {
        batchResult <- replicate(100, simSpacFunction())
        # batchResult <- t(sapply(1:batchSize, simSpacCFunction))
        resultsAll <- c(resultsAll, batchResult)
        # eigValCI <- t.test(resultsAll)$conf.int
        eigValSE <- sd(resultsAll) / sqrt(length(resultsAll))
        eigValCI <- mean(resultsAll) + c(-1,1) *
          qt(0.975, df = length(resultsAll) - 1) * eigValSE
        lambdasInCI <- lambdas[which(lambdas > eigValCI[1] &
                                       lambdas < eigValCI[2])]
        if(length(lambdasInCI) < 2)
        {
          break
        }
      }
    }
    relSpacsIdx <- which(lambdas < mean(resultsAll))
    nSpacs <- if (any(relSpacsIdx)) min(relSpacsIdx) else pcaResults$nEigenVals
    nSpacs
  }


  if (compute_nSpacs) {
    nSpacs <- computeRelevantSpacs(nSim, 10, dataReduced, tmpTrainGL, lambdas)
    SpaCoObject@nSpacs <- nSpacs
  }

  # Compute ONB and projections
  ONBOriginalBasis <- pcaResults$initialPCA$vectors[, 1:pcaResults$nEigenVals] %*% PCsRx
  rownames(ONBOriginalBasis) <- colnames(dataCentered)
  colnames(ONBOriginalBasis) <- paste0("spac_", 1:ncol(ONBOriginalBasis))
  SpaCoObject@spacs <- ONBOriginalBasis
  tmp <-
    scale(eigenMapMatMult(data,
                          pcaResults$initialPCA$vectors[, 1:pcaResults$nEigenVals]))
  SpaCoObject@projection <- eigenMapMatMult(tmp, PCsRx)
  rownames(SpaCoObject@projection) <- rownames(data)
  colnames(SpaCoObject@projection) <- paste0("spac_", 1:ncol(ONBOriginalBasis))
  SpaCoObject@Lambdas <- lambdas
  SpaCoObject@GraphLaplacian <-
    computeGraphLaplacian(SpaCoObject@neighbours)
  return(SpaCoObject)
}
