#' Title
#'
#' @param x
#' @param neighbourindexmatrix
#'
#' @return
#' @export
#'
#' @examples
compute_C <- function(x, neighbourindexmatrix)
{
  n <- nrow(neighbourindexmatrix)
  W <- sum(neighbourindexmatrix)
  activeIndices <- which(neighbourindexmatrix != 0)
  i <- activeIndices %% n
  i[which(i == 0)] <- n
  j <- ceiling(activeIndices / n)
  C <- (((n - 1) / (2 * n * W)) *
          sum(neighbourindexmatrix[activeIndices] * (x[i] - x[j])^2)) / stats::var(x)
  return(C)
}

#' Title
#'
#' @param x
#' @param A
#'
#' @return
#' @export
#'
#' @examples
normalizeA <- function(x, A, preFactor)
{
  return((x / c(normA(x, A, preFactor))))
}
#' Title
#'
#' @param x
#' @param A
#'
#' @return
#' @export
#'
#' @examples
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
#' @param v
#' @param u
#' @param A
#'
#' @return
#' @export
#'
#' @examples
projAFunction <- function(v, u, A, preFactor)
{
  return(c(scalarProductA(v, u, A, preFactor)/
             scalarProductA(u, u, A, preFactor)) * u)
}
#' Title
#'
#' @param v
#' @param ONB
#' @param A
#'
#' @return
#' @export
#'
#' @examples
projASubspaceFunction <- function(v, projMatrix, preFactor)
{
  projection <- preFactor * eigenMapMatMult(projMatrix, v)
  return(projection)
}
#' Title
#'
#' @param X
#' @param A
#' @param tol
#'
#' @return
#' @export
#'
#' @examples
orthogonalizeA <- function(X, A, tol = .Machine$double.eps^0.5)
{
  n <- nrow(A)
  W <- -(1/2)*sum(A - diag(diag(A)))
  preFactor <- (n-1)/(2*n*W)
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



#' Title
#'
#' @param SpaCoObject Spaco object to compute profiles of
#' @param nSpacs number of spac's to use to compute the smoothed profiles. Should be determined by the compute_nSpacs function.
#'
#' @return
#' @export
#'
#' @examples
compute_smoothed_profiles <- function(SpaCoObject, nSpacs) {
smoothed <- SpaCoObject@spacs %*% t(SpaCoObject@projection)
slot(SpaCoObject,"smoothed") <- as.data.frame(t(smoothed))
return(SpaCoObject)
}


