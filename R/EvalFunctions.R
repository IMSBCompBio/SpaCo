#' Compute the spatial variable genes of a SpaCoObject after runinng runSCA
#'
#' @param SpaCoObject SpaCoObject to compute spatially variable genes of.
#' @param adjustMethod method for p-value adjustment. See p.adjust function.
#'
#' @return returns a data frame of spatial variable genes and their p-Values.
#' @export
#'
#'
#' @import mgcv
#'
SVGTest <- function(SpaCoObject, adjustMethod = "holm") {
  require(mgcv)
  if(!is.null(SpaCoObject@GraphLaplacian_B) && (ncol(SpaCoObject@GraphLaplacian_B) > 0))
  {
    GraphLaplacian <- SpaCoObject@GraphLaplacian_B
    projection <- SpaCoObject@projection_B
    projection <-
      .orthogonalizeA(projection, GraphLaplacian, SpaCoObject@nSpacs)
    data <- SpaCoObject@data_B
    S <- projection[,1:SpaCoObject@nSpacs]
  }else
  {
    GraphLaplacian <- SpaCoObject@GraphLaplacian
    projection <- SpaCoObject@projection
    data <- SpaCoObject@data
    S <- sweep(SpaCoObject@projection[,1:SpaCoObject@nSpacs], 2, sqrt(SpaCoObject@Lambdas[1:SpaCoObject@nSpacs]), "/")
  }
  sigma <- eigenMapMatMult(GraphLaplacian, eigenMapMatMult(S, eigenMapMatMult(t(S), GraphLaplacian)))
  sigmaSVD <- eigen(sigma, symmetric = TRUE)

  # Check if @meta.data is not NULL and has at least one column
  if (!is.null(SpaCoObject@meta.data) && ncol(SpaCoObject@meta.data) > 0) {
    COVERAGE <- SpaCoObject@meta.data[rownames(SpaCoObject@data), "nCount_RNA"]
    if (!is.null(COVERAGE)) {
      data <- cbind(data, COVERAGE)
      colnames(data)[ncol(data)] <- "COVERAGE"
    }
  }

  C <- sigmaSVD$values
  getpVal <- function(gene) {
    gene <- gene / c(sqrt((t(gene) %*% GraphLaplacian %*% gene)))
    testStat <- t(gene) %*% sigma %*% gene
    pVal <- psum.chisq(testStat, lb = C[1:SpaCoObject@nSpacs],
                       df = rep(1, SpaCoObject@nSpacs),
                       lower.tail = FALSE)
    return(data.frame(score = testStat, pVal = pVal))
  }

  # Apply the function to each column of the data
  resDf <- t(sapply(1:ncol(data), function(x) getpVal(data[,x])))
  resDf <- as.data.frame(resDf)
  resDf[,1] <- unlist(resDf[,1])
  resDf[,2] <- unlist(resDf[,2])
  rownames(resDf) <- colnames(data)
  resDf[resDf$pVal == 0, "pVal"] <- 2e-25
  resDf$p.adjust = p.adjust(resDf$pVal, method = adjustMethod)

  if (!is.null(SpaCoObject@meta.data) && ncol(SpaCoObject@meta.data) > 0 && resDf["COVERAGE", "p.adjust"] < 0.05) {
    warning("The coverage has been tested as significant")
  }

  return(resDf)
}
