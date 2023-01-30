orthogonalizeA <- function(X, A, tol = .Machine$double.eps^0.5)
{
  #Orthogonalize columns of matrix X regarding scalar product <v,v> = v^TAv
  m <- nrow(X)
  n <- ncol(X)
  if (m < n)
    stop("No. of rows of 'A' must be greater or equal no. of colums.")
  Q <- matrix(0, m, n)
  R <- matrix(0, n, n)
  for (k in 1:n) {
    Q[, k] <- X[, k]
    if (k > 1) {
      for (i in 1:(k - 1)) {
        R[i, k] <- t(Q[, i]) %*% A %*% Q[, k]
        Q[, k] <- Q[, k] - R[i, k] * Q[, i]
      }
    }
    R[k, k] <- normA(Q[, k], A)
    if (abs(R[k, k]) <= tol)
      stop("Matrix 'A' does not have full rank.")
    Q[, k] <- Q[, k]/R[k, k]
  }
  return(list(Q = Q, R = R))
}
