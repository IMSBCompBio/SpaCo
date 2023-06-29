findEdges <- function(SpaCoObject, SVGTestRes, nDendro = 10)
{
  require(mclust)
  SvgMat <- SpaCoObject@data[,which(SVGTestRes$p.adjust < 0.05)] %>% as.matrix()
  LoGMat <- getLoG(SvgMat, sigma = 5, coordinateData = SpaCoObject@coordinates)
  SobelMat <- getSobel(SvgMat, sigma = 3,
                       coordinateData = SpaCoObject@coordinates)

  LoGClusters <- getGmmClusters(LoGMat, epochs = 100, K = 3, lockedMu = 2)
  SobelClusters <- getGmmClusters(SobelMat, epochs = 100, K = 2, lockedMu = 1)

  negProbsLoG <- lapply(LoGClusters, function(x) x[[1]][,1]) %>% unlist() %>%
    matrix(ncol = length(LoGClusters))
  ZeroProbLoG <- lapply(LoGClusters, function(x) x[[1]][,2]) %>% unlist() %>%
    matrix(ncol = length(LoGClusters))
  posProbsSobel <- lapply(SobelClusters, function(x) x[[1]][,2]) %>% unlist() %>%
    matrix(ncol = length(SobelClusters))
  EdgeProb <- ZeroProbLoG * posProbsSobel

  d1 <- 1 - (1/(nrow(negProbsLoG) - 1)) * eigenMapMatMult(t(scale(EdgeProb)),
                                                          scale(EdgeProb))
  d2 <- 1 - (1/(nrow(negProbsLoG) - 1)) * eigenMapMatMult(t(scale(negProbsLoG)),
                                                          scale(negProbsLoG))
  d <- d1 + (1/4) * d2
  diag(d) <- 0
  colnames(d) <- colnames(SvgMat)
  d <- as.dist(d)
  geneClust <- hclust(d)
  geneClust_sub <- cutree(geneClust, k = nDendro)
  return(list(negProbsLoG, ZeroProbLoG, posProbsSobel, EdgeProb, d, geneClust,
              geneClust_sub))
}
expectation_step <- function(X, paramList)
{
  Mu <- paramList[[1]]
  Sigma <- paramList[[2]]
  Pi <- paramList[[3]]
  K <- length(Mu)
  gamma_nk <- sapply(1:K, function(x) Pi[x] * dnorm(X, Mu[x], Sigma[x]))
  totals <- rowSums(gamma_nk)
  gamma_nk <- gamma_nk/ ifelse(totals == 0, 1, totals)
  retList <- list(gamma_nk, totals)
  return(retList)
}
maximization_step <- function(X, gamma_nk, paramList, lockedMu)
{
  Mu <- paramList[[1]]
  Sigma <- paramList[[2]]
  Pi <- paramList[[3]]
  K <- length(Mu)
  N <- length(X)
  for(k  in 1:K)
  {
    gamma_k <- gamma_nk[,k]
    N_k <- sum(gamma_k)
    Pi[k] <- N_k / N
    if(!k %in% lockedMu)
    {
      Mu[k] <- sum(gamma_k * X) / N_k
    }
    Sigma[k] <- min(sqrt((t(gamma_k * (X - c(Mu[k]))) %*% (X - c(Mu[k]))) / N_k),
                    (1/K) * sd(X))
  }
  retList <- list(Mu, Sigma, Pi)
  return(retList)
}
get_likelihood <- function(X, paramList, totals)
{
  sample_likelihoods <- log(totals)
  llh <- sum(sample_likelihoods)
  retList <- list(llh, sample_likelihoods)
  return(retList)
}
train_gmm <- function(X, K, epochs, lockedMu)
{
  initParams <- list(if(K == 2) c(0, mean(X)) else
    c(mean(X[X < 0]), 0, mean(X[X > 0])),
    rep(sd(X)/K, K),
    rep(1/K, K))
  paramList <- initParams
  for(i in 1:epochs)
  {
    expecStepRet <- expectation_step(X, paramList)
    gamma_nk <- expecStepRet[[1]]
    totals <- expecStepRet[[2]]
    paramList <- maximization_step(X, gamma_nk, paramList, lockedMu)
    # llh <- get_likelihood(X, paramList, totals)[[1]]
    # print(llh)
  }
  fittedVals <- sapply(1:K, function(x) paramList[[3]][x] * dnorm(X, paramList[[1]][x], paramList[[2]][x]))
  fittedVals <- fittedVals / rowSums(fittedVals)
  retList <- list(fitted = fittedVals, params = paramList)
  return(retList)
}
getGmmClusters <- function(data, epochs = 100, K = 3, lockedMu = 2)
{
  gmmRes <- apply(data, 2, function(x)
    train_gmm(x, K, epochs, lockedMu))
  return(gmmRes)
}
findNeighbors <- function(sigma)
{
  grid <- expand.grid((-3 * sigma):(3 * sigma), (-3 * sigma):(3 * sigma))
  distances <- sqrt(grid[,1]^2 + grid[,2]^2)
  neighbors <- grid[distances < 3 * sigma, ]
  neighbors <- neighbors[apply(neighbors, 1, sum) %% 2 == 0,]
  return(neighbors)
}
getKernel <- function(i, neighborMat, kernelEntries, coordinateData)
{
  coordinateStr <- apply(coordinateData, 1, paste, collapse = ";")
  x <- coordinateData[i,1]
  y <- coordinateData[i,2]
  retVec <- rep(0, nrow(coordinateData))
  NIdx <-
    apply(neighborMat, 1,
          function(row) which(coordinateStr == paste(x + row[1], y + row[2],
                                                     sep = ";")))
  NonZeroIdx <- lapply(NIdx, function(entry) length(entry) != 0) %>% unlist()
  retVec[NIdx[NonZeroIdx] %>% unlist()] <- kernelEntries[NonZeroIdx]
  # if(length(NIdx) == nrow(neighborMat))
  # {
  #   retVec[NIdx] <- kernelEntries
  # }
  return(retVec)
}
getKernelMatrixLaplace <- function(sigma, coordinateData)
{
  neighborMat <- findNeighbors(sigma)
  kernelEntries <- apply(neighborMat, 1, function(vec)
  {
    tmpFrac <- -(vec[1]^2 + vec[2]^2)/(2*sigma^2)
    return(-(1 / (sigma^4 * pi)) * exp(tmpFrac) * (1 + tmpFrac))
  }) %>% scale()
  kernelMatrix <- sapply(1:nrow(coordinateData), getKernel,
                         neighborMat = neighborMat,
                         kernelEntries = kernelEntries,
                         coordinateData = coordinateData)
  return(kernelMatrix)
}
getKernelMatrixSobel <- function(sigma, coordinateData, dir = "x")
{
  if(!dir %in% c("x", "y"))
  {
    stop("dir must be either x or y")
  }
  neighborMat <- findNeighbors(sigma)
  kernelEntries <- apply(neighborMat, 1, function(vec)
  {
    tmpFrac <- -(vec[1]^2 + vec[2]^2)/(2*sigma^2)
    return(-(1 / (sigma^4 * 2 * pi)) * vec[1 + (dir == "y")] * exp(tmpFrac))
  }) %>% scale()
  kernelMatrix <- sapply(1:nrow(coordinateData), getKernel,
                         neighborMat = neighborMat,
                         kernelEntries = kernelEntries,
                         coordinateData = coordinateData)
  return(kernelMatrix)
}
getLoG <- function(data, sigma, coordinateData)
{
  kernelMatrix <- getKernelMatrixLaplace(sigma, coordinateData)
  LaplaceMat <- eigenMapMatMult(t(kernelMatrix), scale(data))
  return(LaplaceMat)
}
getSobel <- function(data, sigma, coordinateData)
{
  kernelSobelX <- getKernelMatrixSobel(sigma, coordinateData, "x")
  kernelSobelY <- getKernelMatrixSobel(sigma, coordinateData, "y")
  sobelX <- eigenMapMatMult(t(kernelSobelX), scale(data))
  sobelY <- eigenMapMatMult(t(kernelSobelY), scale(data))
  GradientSquare <- sobelX^2 + sobelY^2
  return(GradientSquare)
}
