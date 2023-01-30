scalarProductA <- function(x, y, A)
{
  n <- nrow(A)
  W <- -(1/2)*sum(A - diag(diag(A)))
  preFactor <- (n-1)/(2*n*W)
  return(preFactor * t(x) %*% A %*% y)
}
