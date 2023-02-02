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
  numcores = detectCores()
  cl = makeCluster(numcores - 2)
  clusterExport(cl, list("data_centeredLocal", "GraphLaplacianLocal", "ONBLocal"), envir = environment())
  clusterEvalQ(cl, {
    source("~/SPACO/R/SCA.R")
    # source("LoadPreprocessData.R")
    source("~/SPACO/R/AuxiliaryFunctions.R")
    # source("PlotFunctions.R")
    source("~/SPACO/R/EvalFunctions.R")
  })
  #Apply Gene Score Function to all genes
  results_all <- t(parSapply(cl, 1:ncol(data_centeredLocal), getGeneScoreParallelFunction))
  stopCluster(cl)
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
  numcores = detectCores()
  cl = makeCluster(numcores - 2)
  clusterExport(cl, list("GraphLaplacianLocal", "ONBLocal"),envir = environment())
  clusterEvalQ(cl, {
    source("~/SPACO/R/SCA.R")
    # source("LoadPreprocessData.R")
    source("~/SPACO/R/AuxiliaryFunctions.R")
    # source("PlotFunctions.R")
    source("~/SPACO/R/EvalFunctions.R")
  })
  #Simulate genes and get score
  simScores <- parSapply(cl, rep(nrow(GraphLaplacianLocal), nSimLocal), simGeneScoreFunction)
  stopCluster(cl)
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
  simGene <- rnorm(nSim)
  return(getSingleGeneScoreA(simGene, ONBLocal, GraphLaplacianLocal))
}

