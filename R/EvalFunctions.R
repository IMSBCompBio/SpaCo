#source("AuxiliaryFunctions.R")
library(parallel)
# getSingleGeneScoreAndPVal <- function(gene, A, nSim, preFactor, projMatrix,
#                                       bootstrap = FALSE)
# {
#   gc()
#   #Project gene onto orthonormal basis ONB
#   projection <- projASubspaceFunction(gene, projMatrix, preFactor)
#   gene <- gene/norm(gene, type = "2")
#   projection <- projection/norm(projection, type = "2")
#   score <- as.numeric(normA(gene - projection, A, preFactor))
#   resampledCoeffs <- replicate(nSim,
#                                sample(length(gene), length(gene),
#                                       replace = bootstrap))
#   simGenes <- matrix(gene[resampledCoeffs], ncol = nSim)
#   if(bootstrap)
#   {
#     simGenes <- apply(simGenes, 2, function(x) x / norm(x, type = "2"))
#   }
#   simProjections <- projASubspaceFunction(simGenes, projMatrix, preFactor)
#   simProjections <- apply(simProjections, 2,
#                           function(x) x / norm(x, type = "2"))
#   #overwrite to save storage space
#   simProjections <- simGenes - simProjections
#   simScores <- sqrt(preFactor * colSums(simProjections *
#                                           eigenMapMatMult(A, simProjections)))
#   pVal <- max(1, sum(simScores < score))/nSim
#   return(c(score, pVal))
# }
getSingleGeneScoreAndPVal <- function(geneIdx, data_centered, A, nSim, preFactor, projMatrix,
                                      bootstrap = FALSE)
{
  gc()
  gene <- data_centered[,geneIdx]
  #Project gene onto orthonormal basis ONB
  projection <- projASubspaceFunction(gene, projMatrix, preFactor)
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
# getPScoreParallelFunction <- function(iGene)
# {
#   return(getSingleGeneScoreAndPVal(data_centered[,iGene], GraphLaplacian,
#                                    nSim, preFactor, projMatrix, bootstrap))
# }
# getPScoreParallelWrapper <- function(data_centered, GraphLaplacian, ONB, nSim,
#                                      bootstrap = FALSE)
# {
#   #compute preFactor
#   n <- nrow(GraphLaplacian)
#   W <- -(1/2)*sum(GraphLaplacian - diag(diag(GraphLaplacian)))
#   preFactor <- (n-1)/(2*n*W)
#   #compute projection matrix
#   projMatrix <- eigenMapMatMult(ONB, eigenMapMatMult(t(ONB), GraphLaplacian))
#   #Cluster Setup
#   numcores = detectCores()
#   cl = makeCluster(numcores - 2)
#   clusterExport(cl,
#                 list("data_centered", "GraphLaplacian", "ONB", "nSim",
#                      "preFactor", "projMatrix", "bootstrap"),
#                 envir = environment())
#   parallelDebug <- clusterEvalQ(cl, {
#     source("R/SCA.R")
#     # source("LoadPreprocessData.R")
#     source("R/AuxiliaryFunctions.R")
#     # source("PlotFunctions.R")
#     source("R/EvalFunctions.R")
#     library(rARPACK)
#     library(Rcpp)
#     library(RcppEigen)
#     sourceCpp("Rcpp_Functions.cpp")
#   })
#   rm(parallelDebug)
#   #Apply Gene Score Function to all genes
#   results_all <- t(parSapply(cl, 1:ncol(data_centered), getPScoreParallelFunction))
#   stopCluster(cl)
#   #Create output dataframe
#   genePScoreData <- data.frame(Gene = colnames(data_centered),
#                                score = results_all[,1],
#                                pVal = results_all[,2],
#                                pAdjust = p.adjust(results_all[,2],method = "bonferroni"))
#   return(genePScoreData)
# }
getPScoreSerialWrapper <- function(data_centered, GraphLaplacian, ONB, nSim,
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
                               pAdjust = p.adjust(results_all[,2],method = "bonferroni"))
  return(genePScoreData)
}
