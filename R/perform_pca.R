#' @keywords internal
# PCA and dimension reduction
performPCA <- function(data, criterion, value, n) {
  varMatrix <- (1 / (n - 1)) * eigenMapMatMult(t(data), data)
  initialPCA <- eigen(varMatrix, symmetric = TRUE)
  nEigenVals <- if (criterion == "percent") {
    if (value == 1)
      p
    else
      min(which(
        cumsum(initialPCA$values) / sum(initialPCA$values) > value
      ))
  } else {
    value
  }
  list(
    dataReduced = t(eigenMapMatMult(
      diag(1 / sqrt(initialPCA$values[1:nEigenVals])), eigenMapMatMult(t(initialPCA$vectors[, 1:nEigenVals]), t(data))
    )),
    nEigenVals = nEigenVals,
    initialPCA = initialPCA
  )
}
