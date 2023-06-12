#' wrapper function to compute spatial variable genes and their empirical p-values
#'
#' @param SpaCoObject SpaCoObject to compute spatially variable genes for.
#' @param nSpacs number of spatial components to consider. Should be determined during RunSCA.
#' @param nSim Number of simulations to compute empirical p-Values.
#'
#' @return Returns a list of genes with distance scores and p-values.
#' @export
#'
#' @import mgcv
#' @import CompQuadForm
#'

FindSVG <- function(SpaCoObject, nSpacs, nSim = 1e3)
{
  data <- SpaCoObject@data
  Spacos <- SpaCoObject@spacs
  Spacos <- Spacos[,1:nSpacs]
  n <- nrow(data)
  W <- sum(SpaCoObject@neighbours)
  preFactor <- (1)/(2*W)
  #Compute metagene expression profiles
  SpacoProjection <- t(eigenMapMatMult(t(Spacos), t(data)))
  #Create orthonormal basis for metagene space
  ONB <- .orthogonalizeA(SpacoProjection, SpaCoObject@GraphLaplacian)$Q
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
  preFactor <- (1)/(2*W)
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

getSingleGeneScoreAndPVal <- function(geneIdx, data_centered, A, nSim, preFactor, projMatrix,
                                      bootstrap = FALSE)
{
  gc()
  gene <- data_centered[,geneIdx]
  #Project gene onto orthonormal basis ONB
  projection <- projASubspaceFunction(gene, projMatrix, preFactor)###smothed profile
  gene <- gene/norm(gene, type = "2")
  projection <- projection/norm(projection, type = "2")
  score <- as.numeric(normA(gene - projection, A, preFactor))
  resampledCoeffs <- replicate(nSim,
                               sample(length(gene), length(gene),
                                      replace = bootstrap))
  simGenes <- matrix(gene[resampledCoeffs], ncol = nSim)
  if(bootstrap)
  {
    simGenes <- apply(simGenes, 2, function(x) x / norm(x, type = "2"))
  }
  simProjections <- projASubspaceFunction(simGenes, projMatrix, preFactor)
  simProjections <- apply(simProjections, 2,
                          function(x) x / norm(x, type = "2"))
  #overwrite to save storage space
  simProjections <- simGenes - simProjections
  simScores <- sqrt(preFactor * colSums(simProjections *
                                          eigenMapMatMult(A, simProjections)))
  pVal <- max(1, sum(simScores < score))/nSim
  return(c(score, pVal))
}




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
  Q <- sigmaSVD$vectors
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
  for(i in 1:ncol(SpaCoObject@data)){
    temp <- getpVal(SpaCoObject@data[,i])
    tempRowName <- colnames(SpaCoObject@data)[i]
    rownames(temp) <- tempRowName
    resDf <- rbind(resDf, temp)
  }

  resDf$p.adjust = p.adjust(resDf$pVal, method = adjustMethod)
  return(resDf)

}
