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
      A_x <- matrix(0, nrow(X), ncol(X))
      X_length <- ncol(X)
      apply(1:X_length, MARGIN = 1,
            function(i) {
              neighborLoci <- which(W_Matrix[i,] != 0)
              x_i_bar <- rowSums(X[,neighborLoci] *
                                   matrix(W_Matrix[i, neighborLoci], ncol = length(neighborLoci),
                                          nrow = nrow(X), byrow = TRUE))
              A_x <<- A_x + X[,i] %*% t(x_i_bar)
            })
      A_x <- A_x * (-1 / W)
    }else if(type == "C")
    {
      A_x <- matrix(0, ncol = ncol(X), nrow = ncol(X))
      X_length <- ncol(X)
      apply(1:X_length, 1, function(i) {
        d_i <- sum(W_Matrix[i,])
        neighborLoci <- which(W_Matrix[i,] != 0)
        x_i_bar <- 1 / d_i * rowSums(X[,neighborLoci] *
                                       matrix(W_Matrix[i,neighborLoci], ncol = length(neighborLoci),
                                              nrow = nrow(X), byrow = TRUE))
        A_x <- A_x + 2 * d_i * X[,i] %*% t(X[,i] - x_i_bar)
        return(A_x)
      })
      A_x <- A_x * ((X_length - 1) / (2 * X_length * W))
    }
  }
  return(A_x)
}
