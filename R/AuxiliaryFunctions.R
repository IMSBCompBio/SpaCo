#' computed smoothed gene profiles of genes present in the data.
#'
#' @param SpaCoObject Spaco object to compute profiles of.
#'
#'
#' @return smoothed gene profiles in the SpaCoObject.
#' @export
#'
denoise_profiles <- function(SpaCoObject){
  data <- SpaCoObject@data
  GraphLaplacian <- SpaCoObject@GraphLaplacian

  SpacoProjection <- SpaCoObject@projection[,1:SpaCoObject@nSpacs]
  projMatrix <- eigenMapMatMult(SpacoProjection,
                                eigenMapMatMult(t(SpacoProjection),
                                                GraphLaplacian))
  #Center data regarding A-norm
  data_centered <- scale(data)
  projection <- eigenMapMatMult(projMatrix, data)
  colnames(projection) <- colnames(data_centered)
  rownames(projection) <- rownames(data_centered)
  sds <- apply(projection, 2, sd)
  projection <- sweep(projection, MARGIN = 2, STATS = sds, FUN = "/")
  slot(SpaCoObject, "denoised") <- as.data.frame(projection)
  return(SpaCoObject)
}
#' Title
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
.orthogonalizeA <- function(X, A, nSpacs, tol = .Machine$double.eps^0.5)
{
  n <- nrow(A)
  W <- -(1/2)*sum(A - diag(diag(A)))
  # preFactor <- (1)/(2*W)
  preFactor <- 1
  m <- nrow(X)
  n <- ncol(X)
  if (m < n)
    stop("No. of rows of 'A' must be greater or equal no. of colums.")
  Q <- matrix(0, m, n)
  Norms <- rep(0, n)
  for (k in 1:nSpacs) {
    Q[, k] <- X[, k]
    if (k > 1) {
      for (i in 1:(k - 1)) {
        Q[, k] <- Q[, k] -
          c((Q[, k] %*% A %*% Q[, i]) / c(Q[, i] %*% A %*% Q[, i])) * Q[, i]
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
