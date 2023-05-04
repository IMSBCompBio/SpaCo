#' compute geary's C for a given gene.
#'
#' @param feature gene to compute geary's C for.
#' @param SpaCoObject containing the data.
#'
#' @return returns geary's c value.
#' @export
compute_C <- function(feature = NULL, SpaCoObject) {
  x <- SpaCoObject@data[, feature]
  neighbourindexmatrix <- SpaCoObject@neighbours
  n <- nrow(neighbourindexmatrix)
  W <- sum(neighbourindexmatrix)
  activeIndices <- which(neighbourindexmatrix != 0)
  i <- activeIndices %% n
  i[which(i == 0)] <- n
  j <- ceiling(activeIndices / n)
  C <- (((1) / (2 * W)) *
          sum(neighbourindexmatrix[activeIndices] * (x[i] - x[j])^2)) / stats::var(x)
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
  Spacos_norm <- sweep(SpaCoObject@spacs, 2, sqrt(SpaCoObject@Lambdas), "/")
  Spacos <- Spacos_norm[,1:SpaCoObject@nSpacs]
  n <- nrow(data)
  W <- sum(SpaCoObject@neighbours)
  preFactor <- (1)/(2*W)
  GraphLaplacian <- SpaCoObject@GraphLaplacian

  #Compute metagene expression profiles
  SpacoProjection <- t(eigenMapMatMult(t(Spacos), t(data)))
  ONB <- .orthogonalizeA(SpacoProjection, SpaCoObject@GraphLaplacian)$Q
  projMatrix <- eigenMapMatMult(ONB, eigenMapMatMult(t(ONB), GraphLaplacian))
  #Create orthonormal basis for metagene space

  #Center data regarding A-norm
  data_centered <- scale(apply(data, 2, normalizeA, A = SpaCoObject@GraphLaplacian, preFactor), scale = FALSE)
  projection <- matrix(nrow = nrow(data_centered),ncol=ncol(data_centered))
  for (i in seq_along(1:ncol(data))){
    gene <- data_centered[,i]
    projection[,i] <- projASubspaceFunction(gene, projMatrix, preFactor)###smothed profile
    #projection[,i] <- projection[,i]/norm(projection[,i], type = "2")
  }
  colnames(projection) <- colnames(data_centered)
  rownames(projection) <- rownames(data_centered)
  slot(SpaCoObject, "smoothed") <- as.data.frame(projection)
  return(SpaCoObject)
}
