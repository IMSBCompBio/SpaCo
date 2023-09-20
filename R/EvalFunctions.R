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
#' @import CompQuadForm
#'
SVGTest <- function(SpaCoObject, adjustMethod = "holm")
{
  require(mgcv)
  require(CompQuadForm)
  GraphLaplacian <- SpaCoObject@GraphLaplacian
  S <- sweep(SpaCoObject@projection[,1:SpaCoObject@nSpacs],
             2, sqrt(SpaCoObject@Lambdas[1:SpaCoObject@nSpacs]), "/")
  sigma <- eigenMapMatMult(GraphLaplacian, eigenMapMatMult(S, eigenMapMatMult(t(S), GraphLaplacian)))
  sigmaSVD <- eigen(sigma, symmetric = TRUE)
  C <- sigmaSVD$values
  getpVal <- function(gene)
  {
    gene <- scale(gene)
    testStat <- t(gene) %*% sigma %*% gene

    pVal <- psum.chisq(testStat, lb = C[1:SpaCoObject@nSpacs],
                       df = rep(1, SpaCoObject@nSpacs),
                       lower.tail = FALSE)
    # pVal <- farebrother(testStat, C[1:SpaCoObject@nSpacs])
    return(data.frame(score = testStat, pVal = pVal))
  }
  #pVals <- apply(SpaCoObject@data, 2, getpVal)

  # Initialize an empty data frame
  resDf <- data.frame()
  # Apply the function to each column of the data
  resDf <- t(sapply(1:ncol(SpaCoObject@data),
                    function(x) getpVal(SpaCoObject@data[,x])))
  rownames(resDf) <- colnames(SpaCoObject@data)

  resDf$p.adjust = p.adjust(resDf$pVal, method = adjustMethod)
  return(resDf)

}

