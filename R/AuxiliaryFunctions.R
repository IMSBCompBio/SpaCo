#' compute geary's C for a given gene.
#'
#' @param feature gene to compute geary's C for.
#' @param SpaCoObject containing the data.
#'
#' @return returns geary's c value.
#' @export
#compute_C <- function(feature = NULL, SpaCoObject) {
#  x <- SpaCoObject@data[, feature]
 # neighbourindexmatrix <- SpaCoObject@neighbours
 # n <- nrow(neighbourindexmatrix)
#  W <- sum(neighbourindexmatrix)
#  activeIndices <- which(neighbourindexmatrix != 0,arr.ind = T)
 # i <- activeIndices %% n
 # i[which(i == 0)] <- n
#  j <- ceiling(activeIndices / n)
 # C <- (((1) / (2 * W)) *
 #         sum(neighbourindexmatrix[activeIndices] * (x[i] - x[j])^2)) / stats::var(x)
#  return(C)
#}
compute_C <- function(feature = NULL, SpaCoObject) {
  x <- SpaCoObject@data[, feature]
  neighbourindexmatrix <- SpaCoObject@neighbours

  # Get the indices of non-zero entries in the neighbour matrix
  activeIndices <- which(neighbourindexmatrix != 0, arr.ind = TRUE)

  # Calculate W
  W <- sum(neighbourindexmatrix)

  # Get i and j indices
  i <- activeIndices[, 1]
  j <- activeIndices[, 2]

  # Compute Geary's C
  C <- (1 / (2 * W)) * sum((x[i] - x[j])^2) / stats::var(x)

  return(C)
}








#' normalizeA
#'
#' @param x placeholder
#' @param A placeholder
#'
#' @return placeholder
#' @keywords internal
normalizeA <- function(x, A, preFactor)
{
  return((x / c(normA(x, A, preFactor))))
}

projASubspaceFunction <- function(v, projMatrix, preFactor)
{
  projection <- preFactor * eigenMapMatMult(projMatrix, v)
  return(projection)
}

#' Title
#'
#' @param x placeholder
#' @param A placeholder
#'
#' @return placeholder
#' @keywords internal
#'
#'
normA <- function(x, A, preFactor)
{
  return(sqrt(scalarProductA(x, x, A, preFactor)))
}
scalarProductA <- function(x, y, A, preFactor)
{
  return(preFactor * t(x) %*% A %*% y)
}
#' Title
#'
#' @param v placeholder
#' @param u placeholder
#' @param A placeholder
#'
#' @return placeholder
#' @keywords internal
#'
#'
projAFunction <- function(v, u, A, preFactor)
{
  return(c(scalarProductA(v, u, A, preFactor)/
             scalarProductA(u, u, A, preFactor)) * u)
}

#' Title
#'
#' @param X Projection to orthogonalize
#' @param A GraphLaplacian
#' @param tol tolerance
#'
#' @return An orthogonolized A matrix.
#'
#' @keywords internal
#'
.orthogonalizeA <- function(X, A, tol = .Machine$double.eps^0.5)
{
  n <- nrow(A)
  W <- -(1/2)*sum(A - diag(diag(A)))
  preFactor <- (1)/(2*W)
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
        Q[, k] <- Q[, k] - projAFunction(Q[, k], Q[, i], A, preFactor)
      }
    }
    Norms[k]  <- normA(Q[, k], A, preFactor)
    if (abs(Norms[k]) <= tol)
      stop("Matrix 'A' does not have full rank.")
    Q[, k] <- Q[, k] / c(Norms[k])
  }
  return(list(Q = Q))
}
.orthogonalizeA2 <- function(X, A, tol = .Machine$double.eps^0.5)
{
  n <- nrow(A)
  W <- -(1/2)*sum(A - diag(diag(A)))
  preFactor <- (1)/(2*W)
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
        Q[, k] <- Q[, k] - projAFunction(Q[, k], Q[, i], A, preFactor)
      }
    }
    Norms[k]  <- normA(Q[, k], A, preFactor)
    if (abs(Norms[k]) <= tol)
      stop("Matrix 'A' does not have full rank.")
    # Q[, k] <- Q[, k] / c(Norms[k])
  }
  return(list(Q = Q))
}



#' computed smoothed gene profiles of genes present in the data.
#'
#' @param SpaCoObject Spaco object to compute profiles of.
#'
#'
#' @return smoothed gene profiles in the SpaCoObject.
#' @export
#'
smooth_profiles <- function(SpaCoObject){
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
  slot(SpaCoObject, "denoised") <- as.data.frame(projection)
  return(SpaCoObject)
}
