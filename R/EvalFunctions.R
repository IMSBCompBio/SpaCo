#' wrapper function to compute spatial variable genes and their empirical p-values
#'
#' @param SpaCoObject
#' @param nSpacs
#' @param nSim
#'
#' @return
#' @export
#'
#'

FindSVG <- function(SpaCoObject, nSpacs, nSim = 1e3)
{
  data <- SpaCoObject@data
  Spacos <- SpaCoObject@spacs
  Spacos <- Spacos[,1:nSpacs]
  n <- nrow(data)
  W <- sum(SpaCoObject@neighbours)
  preFactor <- (n-1)/(2*n*W)
  #Compute metagene expression profiles
  SpacoProjection <- t(eigenMapMatMult(t(Spacos), t(data)))
  #Create orthonormal basis for metagene space
  ONB <- orthogonalizeA(SpacoProjection, SpaCoObject@GraphLaplacian)$Q
  #Center data regarding A-norm
  data_centered <- scale(apply(data, 2, normalizeA, A = SpaCoObject@GraphLaplacian, preFactor), scale = FALSE)
  message("computing emprirical p-values this may take a while.")
  GeneScoresDataP <- .getPScoreSerialWrapper(data_centered, SpaCoObject@GraphLaplacian, ONB,
                                            nSim = nSim, TRUE)
  return(GeneScoresDataP)
}

.getPScoreSerialWrapper <- function(data_centered, GraphLaplacian, ONB, nSim,
                                   bootstrap = FALSE)
{
  #compute preFactor
  n <- nrow(GraphLaplacian)
  W <- -(1/2)*sum(GraphLaplacian - diag(diag(GraphLaplacian)))
  preFactor <- (n-1)/(2*n*W)
  #compute projection matrix
  projMatrix <- eigenMapMatMult(ONB, eigenMapMatMult(t(ONB), GraphLaplacian))
  #Apply Gene Score Function to all genes
  results_all <- t(sapply(1:ncol(data_centered), getSingleGeneScoreAndPVal,
                          data_centered = data_centered, A = GraphLaplacian,
                          nSim = nSim, preFactor = preFactor, projMatrix = projMatrix,
                          bootstrap = bootstrap))
  #Create output dataframe
  genePScoreData <- data.frame(Gene = colnames(data_centered),
                               score = results_all[,1],
                               pVal = results_all[,2],
                               pAdjust = p.adjust(results_all[,2],method = "hochberg"))
  return(genePScoreData)
}
