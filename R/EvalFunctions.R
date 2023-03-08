#source("AuxiliaryFunctions.R")
library(parallel)
getSingleGeneScoreA <- function(gene, ONB, A)
{
  #Project gene onto orthonormal basis ONB
  projection <- projASubspaceFunction(gene, ONB, A)
  # projection <- linearproj(ONB, gene)$Q
  gene <- normalizeA(gene, A)
  projection <- normalizeA(projection, A)
  score <- acos(abs(scalarProductA(gene, projection, A)))
  return(score)
}


getGeneScoreParallelWrapper <- function(data_centeredLocal, GraphLaplacianLocal, ONBLocal)
{
  #Cluster Setup
  numcores = parallel::detectCores()
  cl = parallel::makeCluster(numcores - 2)
  parallel::clusterExport(cl, list("data_centeredLocal", "GraphLaplacianLocal", "ONBLocal"), envir = environment())
  parallel::clusterEvalQ(cl, {
    source("~/SPACO/R/SCA.R")
    # source("LoadPreprocessData.R")
    source("~/SPACO/R/AuxiliaryFunctions.R")
    # source("PlotFunctions.R")
    source("~/SPACO/R/EvalFunctions.R")
  })
  #Apply Gene Score Function to all genes
  results_all <- t(parallel::parSapply(cl, 1:ncol(data_centeredLocal), getGeneScoreParallelFunction))
  parallel::stopCluster(cl)
  #Create output dataframe
  GeneScoresData <- data.frame(Gene = colnames(data_centeredLocal), score = t(results_all))
  return(GeneScoresData)
}
getGeneScoreParallelFunction <- function(i)
{
  return(getSingleGeneScoreA(data_centeredLocal[,i], ONBLocal, GraphLaplacianLocal))
}
getPValuesParallelWrapper <- function(GraphLaplacianLocal, ONBLocal,
                                      nSimLocal, GeneScoresDataLocal)
{
  #Cluster Setup
  numcores = parallel::detectCores()
  cl = parallel::makeCluster(numcores - 2)
  parallel::clusterExport(cl, list("GraphLaplacianLocal", "ONBLocal"),envir = environment())
  parallel::clusterEvalQ(cl, {
    source("~/SPACO/R/SCA.R")
    # source("LoadPreprocessData.R")
    source("~/SPACO/R/AuxiliaryFunctions.R")
    # source("PlotFunctions.R")
    source("~/SPACO/R/EvalFunctions.R")
  })
  #Simulate genes and get score
  simScores <- parallel::parSapply(cl, rep(nrow(GraphLaplacianLocal), nSimLocal), simGeneScoreFunction)
  parallel::stopCluster(cl)
  #Order simulated scores and compute p-values as fraction of simulated genes
  #with lower score
  getPVal <- function(geneScore, simScores)
  {
    return(sum(simScores < geneScore)/length(simScores))
  }
  GeneScoresData$pVal <- sapply(GeneScoresData$score, getPVal, simScores = simScores)
  return(GeneScoresData)
}
simGeneScoreFunction <- function(nSim)
{
  simGene <- stats::rnorm(nSim)
  return(getSingleGeneScoreA(simGene, ONBLocal, GraphLaplacianLocal))
}

getSingleGeneScoreAndPVal <- function(gene, A, nSim, preFactor, projMatrix,
                                      bootstrap = FALSE)
{
  gc()
  #Project gene onto orthonormal basis ONB
  projection <- projASubspaceFunction(gene, projMatrix, preFactor)
  gene <- gene/norm(gene, type = "2")
  projection <- projection/norm(projection, type = "2")
  score <- as.numeric(normA(gene - projection, A))
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
                                          A %*% simProjections))
  pVal <- max(1, sum(simScores < score))/nSim
  return(c(score, pVal))
}
getPScoreParallelFunction <- function(iGene)
{
  return(getSingleGeneScoreAndPVal(data_centered[,iGene], GraphLaplacian,
                                   nSim, preFactor, projMatrix, bootstrap))
}
getPScoreParallelWrapper <- function(data_centered, GraphLaplacian, ONB, nSim,
                                     bootstrap = FALSE)
{
  #compute preFactor
  n <- nrow(GraphLaplacian)
  W <- -(1/2)*sum(GraphLaplacian - diag(diag(GraphLaplacian)))
  preFactor <- (n-1)/(2*n*W)
  #compute projection matrix
  projMatrix <- ONB %*% t(ONB) %*% GraphLaplacian
  #Cluster Setup
  numcores = detectCores()
  cl = makeCluster(numcores - 2)
  clusterExport(cl,
                list("data_centered", "GraphLaplacian", "ONB", "nSim",
                     "preFactor", "projMatrix", "bootstrap"),
                envir = environment())
  parallelDebug <- clusterEvalQ(cl, {
    source("~/SPACO/R/SCA.R")
    # source("LoadPreprocessData.R")
    source("~/SPACO/R/AuxiliaryFunctions.R")
    # source("PlotFunctions.R")
    source("~/SPACO/R/EvalFunctions.R")
  })
  rm(parallelDebug)
  #Apply Gene Score Function to all genes
  results_all <- t(parSapply(cl, 1:ncol(data_centered), getPScoreParallelFunction))
  stopCluster(cl)
  #Create output dataframe
  genePScoreData <- data.frame(Gene = colnames(data_centered),
                               score = results_all[,1],
                               pVal = results_all[,2],
                               pAdjust = p.adjust(results_all[,2],method = "bonferroni"))
  return(genePScoreData)
}
getProjDistance <- function(v, projMatrix, GraphLaplacian, preFactor)
{
  v <- apply(v, 2, function(x) x / norm(x, type = "2"))
  vProj <- projMatrix %*% v
  vProj <- scale(vProj, scale = FALSE)
  vProj <- apply(vProj, 2, function(x) x / norm(x, type = "2"))
  vProj <- v - vProj
  distance <- sqrt(preFactor * colSums(v * GraphLaplacian %*% vProj))
  return(distance)
}
