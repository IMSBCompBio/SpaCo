#' computes smoothed gene profiles of genes present in the data.
#'
#' @param SpaCoObject Spaco object to compute profiles of.
#'
#'
#' @return smoothed gene profiles in the SpaCoObject.
#' @export
#'
denoise_profiles <- function(SpaCoObject) {
  data <- SpaCoObject@data
  GraphLaplacian <- SpaCoObject@GraphLaplacian

  SpacoProjection <- SpaCoObject@projection[, 1:SpaCoObject@nSpacs]
  if (is(GraphLaplacian, "dgCMatrix"))
  {
    projMatrix <-
      SpacoProjection %*% t(SpacoProjection) %*% GraphLaplacian
  } else
  {
    projMatrix <- eigenMapMatMult(SpacoProjection,
                                  eigenMapMatMult(t(SpacoProjection),
                                                  GraphLaplacian))
  }
  #Center data regarding A-norm
  data_centered <- scale(data)
  if (class(projMatrix) %in% c("dgCMatrix", "dgeMatrix"))
  {
    projection <- projMatrix %*% data
  } else
  {
    projection <- eigenMapMatMult(projMatrix, data)
  }
  colnames(projection) <- colnames(data_centered)
  rownames(projection) <- rownames(data_centered)
  sds <- apply(projection, 2, stats::sd)
  projection <- sweep(projection,
                      MARGIN = 2,
                      STATS = sds,
                      FUN = "/")
  slot(SpaCoObject, "denoised") <- as.data.frame(projection)
  return(SpaCoObject)
}
#' Orthogonalization
#'
#' @param X Projection to orthogonalize
#' @param A GraphLaplacian
#' @param nSpacs number of Spacs to orthogonalize
#' @param tol tolerance
#'
#' @return An orthogonolized A matrix.
#'
#' @keywords internal
#'
.orthogonalizeA <-
  function(X, A, nSpacs, tol = .Machine$double.eps ^ 0.5)
  {
    preFactor <- 1
    m <- nrow(X)
    n <- ncol(X)
    if (m < nSpacs)
      stop("No. of rows of 'A' must be greater or equal no. of colums.")
    Q <- matrix(0, m, n)
    Norms <- rep(0, n)
    for (k in 1:nSpacs) {
      Q[, k] <- X[, k]
      if (k > 1) {
        for (i in 1:(k - 1)) {
          Q[, k] <- Q[, k] -
            rep((Q[, k] %*% A %*% Q[, i]) / (Q[, i] %*% A %*% Q[, i]), m) * Q[, i]
        }
      }
      Norms[k]  <- sqrt(Q[, k] %*% A %*% Q[, k])
      if (abs(Norms[k]) <= tol)
      {
        stop("Matrix 'A' does not have full rank.")
      }
      Q[, k] <- Q[, k] / c(Norms[k])
    }
    return(Q)
  }
#' @keywords internal
computeGraphLaplacian <- function(neighbourIndexMatrix) {
  if (is(neighbourIndexMatrix, "dgCMatrix"))
  {
    W <- sum(neighbourIndexMatrix@x)
    n <- neighbourIndexMatrix@Dim[1]
    neighbourIndexMatrix@x <- neighbourIndexMatrix@x / W
    graphLaplacian <-
      neighbourIndexMatrix + Matrix::Diagonal(n, 1 / n)
  } else
  {
    W <- sum(neighbourIndexMatrix)
    n <- nrow(neighbourIndexMatrix)
    neighbourIndexMatrix <- neighbourIndexMatrix / W
    graphLaplacian <- neighbourIndexMatrix + diag(1 / n, n)
  }
  graphLaplacian
}
#' @keywords internal
# Compute number of relevant SPACs if required
computeRelevantSpacs <-
  function(nSim,
           batchSize = 10,
           dataReduced,
           graphLaplacian,
           lambdas,
           n_pcs = 0) {
    simSpacFunction <- function(i) {
      shuffleOrder <- sample(ncol(graphLaplacian), ncol(graphLaplacian))
      # RxShuffled <-
      #   t(dataReduced[shuffleOrder, ]) %*% graphLaplacian %*% dataReduced[shuffleOrder, ]
      if (is(graphLaplacian, "dgCMatrix"))
      {
        RxShuffled <-
          t(dataReduced[shuffleOrder, ]) %*% graphLaplacian %*% dataReduced[shuffleOrder, ]
      } else
      {
        RxShuffled <- eigenMapMatMult(t(dataReduced[shuffleOrder, ]),
                                      eigenMapMatMult(graphLaplacian, dataReduced[shuffleOrder, ]))
      }
      rARPACK::eigs_sym(RxShuffled, 1, which = "LM")$values
    }
    resultsAll <- replicate(100, simSpacFunction())
    eigValSE <- stats::sd(resultsAll) / sqrt(length(resultsAll))
    eigValCI <-
      mean(resultsAll) +
      stats::qt(0.975, df = length(resultsAll) - 1) * eigValSE * c(-1, 1)
    lambdasInCI <-
      lambdas[lambdas > eigValCI[1] & lambdas < eigValCI[2]]
    if (length(lambdasInCI) > 1)
    {
      for (i in 1:round((nSim - 100) / batchSize))
      {
        batchResult <- replicate(batchSize, simSpacFunction())
        # batchResult <- t(sapply(1:batchSize, simSpacCFunction))
        resultsAll <- c(resultsAll, batchResult)
        # eigValCI <- t.test(resultsAll)$conf.int
        eigValSE <- stats::sd(resultsAll) / sqrt(length(resultsAll))
        eigValCI <- mean(resultsAll) + c(-1, 1) *
          stats::qt(0.975, df = length(resultsAll) - 1) * eigValSE
        lambdasInCI <- lambdas[which(lambdas > eigValCI[1] &
                                       lambdas < eigValCI[2])]
        if (length(lambdasInCI) < 2)
        {
          break
        }
      }
    }
    relSpacsIdx <- which(lambdas < mean(resultsAll))
    nSpacs <-
      if (any(relSpacsIdx))
        min(relSpacsIdx)
    else
      n_pcs
    nSpacs
  }
#' @keywords internal
# PCA and dimension reduction
performPCA <- function(data, criterion, value, n, p) {
  varMatrix <- (1 / (n - 1)) * eigenMapMatMult(t(data), data)
  initialPCA <- eigen(varMatrix, symmetric = TRUE)
  nEigenVals <- if (criterion == "percent") {
    if (value == 1)
      p
    else
      min(which(cumsum(initialPCA$values) / sum(initialPCA$values) > value))
  } else {
    value
  }
  list(
    dataReduced = t(eigenMapMatMult(diag(
      1 / sqrt(initialPCA$values[1:nEigenVals])
    ), eigenMapMatMult(
      t(initialPCA$vectors[, 1:nEigenVals]), t(data)
    ))),
    nEigenVals = nEigenVals,
    initialPCA = initialPCA
  )
}
