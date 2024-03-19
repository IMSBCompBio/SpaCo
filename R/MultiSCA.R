#' Multi slide SCA_function
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
RunMultiSCA <- function(SpaCoObject,
                   PC_criterion = "percent",
                   PC_value = 0.9,
                   compute_nSpacs = FALSE,
                   nSim = 1000,
                   nSpacQuantile = 0.05) {
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
  dataCentered <- SpaCoObject@data$merged
  neighbourIndexMatrix <- SpaCoObject@neighbours$merged
  n <- nrow(dataCentered)
  p <- ncol(dataCentered)
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
  computeGraphLaplacian <- function(neighbourIndexMatrix, n) {
    W <- sum(neighbourIndexMatrix)
    graphLaplacian <- as.matrix(neighbourIndexMatrix / W)
    graphLaplacian <- graphLaplacian + diag(1 / n, nrow = n)
  }

  graphLaplacian <- computeGraphLaplacian(neighbourIndexMatrix, n)

  # Data preprocessing steps
  geneANorms <- colSums(dataCentered * eigenMapMatMult(graphLaplacian, dataCentered))
  dataCenteredGLScaled <- sweep(dataCentered, 2, geneANorms, FUN = "*")

  # PCA and dimension reduction
  performPCA <- function(data, criterion, value) {
    varMatrix <- (1 / (n - 1)) * eigenMapMatMult(t(data), data)
    initialPCA <- eigen(varMatrix, symmetric = TRUE)
    nEigenVals <- if (criterion == "percent") {
      if (value == 1) p else min(which(cumsum(initialPCA$values) / sum(initialPCA$values) > value))
    } else {
      value
    }
    list(dataReduced = eigenMapMatMult(data, initialPCA$vectors[, 1:nEigenVals]),
         nEigenVals = nEigenVals,
         initialPCA = initialPCA)
  }

  pcaResults <- performPCA(dataCenteredGLScaled, PC_criterion, PC_value)
  dataReduced <- scale(pcaResults$dataReduced)

  # Compute test statistic matrix
  Rx <- eigenMapMatMult(t(dataReduced), eigenMapMatMult(graphLaplacian, dataReduced))

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
    if(length(lambdasInCI) > 1) {
      totalIterations = round((nSim - 100) / batchSize)
      for(i in 1:totalIterations) {
        batchResult <- replicate(batchSize, simSpacFunction()) # Adjusted to use batchSize instead of hardcoded 100
        resultsAll <- c(resultsAll, batchResult)
        eigValSE <- sd(resultsAll) / sqrt(length(resultsAll))
        eigValCI <- mean(resultsAll) + c(-1, 1) *
          qt(0.975, df = length(resultsAll) - 1) * eigValSE
        lambdasInCI <- lambdas[which(lambdas > eigValCI[1] & lambdas < eigValCI[2])]
        if(length(lambdasInCI) < 2) {
          break
        }
        if(i %% 10 == 0) { # Check if the loop iteration is divisible by 10
          message(sprintf("%d/%d simulations done.", i * batchSize + 100, nSim))
        }
      }
    }
    relSpacsIdx <- which(lambdas < mean(resultsAll))
    nSpacs <- if (any(relSpacsIdx)) min(relSpacsIdx) else pcaResults$nEigenVals
    nSpacs
  }



  if (compute_nSpacs) {
    nSpacs <- computeRelevantSpacs(nSim, 10, dataReduced, graphLaplacian, lambdas)
    SpaCoObject@nSpacs <- nSpacs
  }

  # Compute ONB and projections
  ONBOriginalBasis <- pcaResults$initialPCA$vectors[, 1:pcaResults$nEigenVals] %*% PCsRx
  rownames(ONBOriginalBasis) <- colnames(dataCentered)
  colnames(ONBOriginalBasis) <- paste0("spac_", 1:ncol(ONBOriginalBasis))
  SpaCoObject@spacs <- ONBOriginalBasis
  SpaCoObject@projection$merged <- eigenMapMatMult(dataReduced, PCsRx)
  rownames(SpaCoObject@projection$merged) <- rownames(dataCentered)
  colnames(SpaCoObject@projection$merged) <- paste0("spac_", 1:ncol(ONBOriginalBasis))
  SpaCoObject@Lambdas <- lambdas
  SpaCoObject@GraphLaplacian <- graphLaplacian

  return(SpaCoObject)
}
