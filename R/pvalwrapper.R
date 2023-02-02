PValWrapperFunction <- function(data, SpacoResult, GraphLaplacian, nSim = 1e3)
{
  Spacos <- SpacoResult[[1]]
  Spacos <- Spacos[,1:max(which(SpacoResult[[2]] < 0.5))]
  #Compute metagene expression profiles
  SpacoProjection <- t(t(Spacos) %*% t(data))
  #Create orthonormal basis for metagene space
  ONB <- orthogonalizeA(SpacoProjection, GraphLaplacian)$Q
  #Center data regarding A-norm
  data_centered <- apply(data, 2, normalizeA, A = GraphLaplacian)

  GeneScoresData <- getGeneScoreParallelWrapper(data_centered, GraphLaplacian, ONB)
  GeneScoresDataP <- getPValuesParallelWrapper(GraphLaplacian, ONB,
                                               nSim = nSim, GeneScoresData)
  return(GeneScoresDataP)
}

PValWrapperFunction_object <- function(SpaCoObject, nSim = 1e3)
{
  data <- SpaCoObject@data
  Spacos <- SpaCoObject@spacs
  Spacos <- Spacos[,1:max(which(SpaCoObject@Lambdas < 0.5))]
  #Compute metagene expression profiles
  SpacoProjection <- t(t(Spacos) %*% t(data))
  #Create orthonormal basis for metagene space
  ONB <- orthogonalizeA(SpacoProjection, SpaCoObject@GraphLaplacian)$Q
  #Center data regarding A-norm
  data_centered <- apply(data, 2, normalizeA, A = SpaCoObject@GraphLaplacian)

  GeneScoresData <- getGeneScoreParallelWrapper(data_centered, SpaCoObject@GraphLaplacian, ONB)
  GeneScoresDataP <- getPValuesParallelWrapper(SpaCoObject@GraphLaplacian, ONB,
                                               nSim = nSim, GeneScoresData)
  return(GeneScoresDataP)
}
