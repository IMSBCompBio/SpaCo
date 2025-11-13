#' @keywords internal
computeGraphLaplacian <- function(neighbourIndexMatrix) {
  # W <- sum(neighbourIndexMatrix@x)
  # n <- neighbourIndexMatrix@Dim[1]
  # neighbourIndexMatrix@x <- neighbourIndexMatrix@x / W
  # neighbourIndexMatrix + Matrix::Diagonal(n, 1 / n)
  if (is(neighbourIndexMatrix, "dgCMatrix"))
  {
    W <- sum(neighbourIndexMatrix@x)
    n <- neighbourIndexMatrix@Dim[1]
    neighbourIndexMatrix@x <- neighbourIndexMatrix@x / W
    graphLaplacian <-
      neighbourIndexMatrix + Matrix::Diagonal(n, 1 / n)
  } else
  {
    W <- sum(neighbourIndexMatrix)
    n <- nrow(neighbourIndexMatrix)
    neighbourIndexMatrix <- neighbourIndexMatrix / W
    graphLaplacian <- neighbourIndexMatrix + diag(1 / n, n)
  }
  graphLaplacian
}
