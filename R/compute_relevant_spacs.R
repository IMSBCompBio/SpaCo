#' @keywords internal
# Compute number of relevant SPACs if required
computeRelevantSpacs <-
  function(nSim,
           batchSize = 10,
           dataReduced,
           graphLaplacian,
           lambdas) {
    simSpacFunction <- function(i) {
      shuffleOrder <- sample(ncol(graphLaplacian), ncol(graphLaplacian))
      # RxShuffled <-
      #   t(dataReduced[shuffleOrder, ]) %*% graphLaplacian %*% dataReduced[shuffleOrder, ]
      if (is(graphLaplacian, "dgCMatrix"))
      {
        RxShuffled <-
          t(dataReduced[shuffleOrder,]) %*% graphLaplacian %*% dataReduced[shuffleOrder,]
      } else
      {
        RxShuffled <- eigenMapMatMult(t(dataReduced[shuffleOrder,]),
                                      eigenMapMatMult(graphLaplacian, dataReduced[shuffleOrder,]))
      }
      eigs_sym(RxShuffled, 1, which = "LM")$values
    }
    resultsAll <- replicate(100, simSpacFunction())
    eigValSE <- sd(resultsAll) / sqrt(length(resultsAll))
    eigValCI <-
      mean(resultsAll) + qt(0.975, df = length(resultsAll) - 1) * eigValSE * c(-1, 1)
    lambdasInCI <-
      lambdas[lambdas > eigValCI[1] & lambdas < eigValCI[2]]
    if (length(lambdasInCI) > 1)
    {
      for (i in 1:round((nSim - 100) / batchSize))
      {
        batchResult <- replicate(batchSize, simSpacFunction())
        # batchResult <- t(sapply(1:batchSize, simSpacCFunction))
        resultsAll <- c(resultsAll, batchResult)
        # eigValCI <- t.test(resultsAll)$conf.int
        eigValSE <- sd(resultsAll) / sqrt(length(resultsAll))
        eigValCI <- mean(resultsAll) + c(-1, 1) *
          qt(0.975, df = length(resultsAll) - 1) * eigValSE
        lambdasInCI <- lambdas[which(lambdas > eigValCI[1] &
                                       lambdas < eigValCI[2])]
        if (length(lambdasInCI) < 2)
        {
          break
        }
      }
    }
    relSpacsIdx <- which(lambdas < mean(resultsAll))
    nSpacs <-
      if (any(relSpacsIdx))
        min(relSpacsIdx)
    else
      pcaResults$nEigenVals
    nSpacs
  }
