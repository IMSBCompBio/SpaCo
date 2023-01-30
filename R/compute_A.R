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
