PValWrapperFunction <- function(data, SpacoResult, nSpacs, GraphLaplacian, nSim = 1e3)
{
  Spacos <- SpacoResult[[1]]
  Spacos <- Spacos[,1:nSpacs]
  #Compute metagene expression profiles
  SpacoProjection <- t(eigenMapMatMult(t(Spacos), t(data)))
  #Create orthonormal basis for metagene space
  ONB <- orthogonalizeA(SpacoProjection, GraphLaplacian)$Q
  #Center data regarding A-norm
  data_centered <- apply(data, 2, normalizeA, A = GraphLaplacian)

  GeneScoresData <- getGeneScoreParallelWrapper(data_centered, GraphLaplacian, ONB)
  GeneScoresDataP <- getPValuesParallelWrapper(GraphLaplacian, ONB,
                                               nSim = nSim, GeneScoresData)
  return(GeneScoresDataP)
}
PValWrapperFunction_object <- function(SpaCoObject, nSpacs, nSim = 1e3)
{
  sourceCpp("Rcpp_Functions.cpp")
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

  GeneScoresDataP <- getPScoreSerialWrapper(data_centered, SpaCoObject@GraphLaplacian, ONB,
                                              nSim = nSim, TRUE)
  return(GeneScoresDataP)
}
