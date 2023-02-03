library(lsa)
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
          sum(neighbourindexmatrix[activeIndices] * (x[i] - x[j])^2)) / var(x)
  return(C)
}
#' Title
#'
#' @param X
#' @param type
#' @param W_Matrix
#'
#' @return
#' @export
#'
#' @examples
compute_A <- function(X, type = "C", W_Matrix) {
  stopifnot(type %in% c("C", "I", "1/V"))
  X_length <- ncol(X)
  A_x <- matrix(0, nrow = nrow(X), ncol = nrow(X))
  W <- sum(W_Matrix)
  if(type == "C") {
    res <- Map(function(i) {
      d_i <- sum(W_Matrix[i,])
      neighborLoci <- which(W_Matrix[i,] != 0)
      x_i_bar <- 1 / d_i * rowSums(X[,neighborLoci] *
                                     matrix(W_Matrix[i,neighborLoci], ncol = length(neighborLoci),
                                            nrow = nrow(X), byrow = TRUE))
      2 * d_i * X[,i] %*% t(X[,i] - x_i_bar)
    }, 1:X_length)
    A_x <- Reduce("+", res)
    A_x <- A_x * ((X_length - 1) / (2 * X_length * W))
  }

  return(A_x)
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
normalizeA <- function(x, A)
{
  return((x / rep(normA(x, A), length(x))))
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
normA <- function(x, A)
{
  return(sqrt(scalarProductA(x, x, A)))
}
scalarProductA <- function(x, y, A)
{
  n <- nrow(A)
  W <- -(1/2)*sum(A - diag(diag(A)))
  preFactor <- (n-1)/(2*n*W)
  return(preFactor * (t(x) %*% A) %*% y)
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
projAFunction <- function(v, u, A)
{
  return(c(scalarProductA(v, u, A)/scalarProductA(u, u, A)) * u)
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
projASubspaceFunction <- function(v, ONB, A)
{
  n <- nrow(A)
  W <- -(1/2)*sum(A - diag(diag(A)))
  preFactor <- (n-1)/(2*n*W)
  projectionCoeffs <- preFactor * t(v) %*% A %*% ONB
  projection <- ONB %*% t(projectionCoeffs)
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
        Q[, k] <- Q[, k] - projAFunction(Q[, k], Q[, i], A)
      }
    }
    Norms[k]  <- normA(Q[, k], A)
    if (abs(Norms[k]) <= tol)
      stop("Matrix 'A' does not have full rank.")
    Q[, k] <- Q[, k] / c(Norms[k])
  }
  return(list(Q = Q))
}
