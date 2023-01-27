library(lsa)
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
compute_A <- function(X, type = "C", W_Matrix)
{
  if(!type %in% c("I", "C", "1/V"))
  {
    stop("type must be either \"I\", \"C\" or \"1/V\"")
  }
  if(type == "1/V")
  {
    A_x <- diag(1, nrow = nrow(X))
  }else
  {
    X_length <- ncol(X)
    A_x <- matrix(0, nrow = nrow(X), ncol = nrow(X))
    W <- sum(W_Matrix)
    if(type == "I")
    {
      for(i in 1:X_length)
      {
        neighborLoci <- which(W_Matrix[i,] != 0)
        x_i_bar <- rowSums(X[,neighborLoci] *
                             matrix(W_Matrix[i, neighborLoci], ncol = length(neighborLoci),
                                    nrow = nrow(X), byrow = TRUE))
        A_x <- A_x + X[,i] %*% t(x_i_bar)
      }
      A_x <- A_x * (-1 / W)
    }else if(type == "C")
    {
      for(i in 1:X_length)
      {
        d_i <- sum(W_Matrix[i,])
        neighborLoci <- which(W_Matrix[i,] != 0)
        x_i_bar <- 1 / d_i * rowSums(X[,neighborLoci] *
                                       matrix(W_Matrix[i,neighborLoci], ncol = length(neighborLoci),
                                              nrow = nrow(X), byrow = TRUE))
        A_x <- A_x + 2 * d_i * X[,i] %*% t(X[,i] - x_i_bar)
      }
      A_x <- A_x * ((X_length - 1) / (2 * X_length * W))
    }
  }
  return(A_x)
}
normalizeA <- function(x, A)
{
  return((x / rep(normA(x, A), length(x))))
}
normA <- function(x, A)
{
  return(sqrt(scalarProductA(x, x, A)))
}
scalarProductA <- function(x, y, A)
{
  n <- nrow(A)
  W <- -(1/2)*sum(A - diag(diag(A)))
  preFactor <- (n-1)/(2*n*W)
  return(preFactor * t(x) %*% A %*% y)
}
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

