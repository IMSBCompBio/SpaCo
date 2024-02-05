#' computed smoothed gene profiles of genes present in the data.
#'
#' @param SpaCoObject Spaco object to compute profiles of.
#'
#'
#' @return smoothed gene profiles in the SpaCoObject.
#' @export
#'
denoise_profiles <- function(SpaCoObject){
  data <- SpaCoObject@data
  GraphLaplacian <- SpaCoObject@GraphLaplacian

  SpacoProjection <- SpaCoObject@projection[,1:SpaCoObject@nSpacs]
  projMatrix <- eigenMapMatMult(SpacoProjection,
                                eigenMapMatMult(t(SpacoProjection),
                                                GraphLaplacian))
  #Center data regarding A-norm
  data_centered <- scale(data)
  projection <- eigenMapMatMult(projMatrix, data)
  colnames(projection) <- colnames(data_centered)
  rownames(projection) <- rownames(data_centered)
  sds <- apply(projection, 2, sd)
  projection <- sweep(projection, MARGIN = 2, STATS = sds, FUN = "/")
  slot(SpaCoObject, "denoised") <- as.data.frame(projection)
  return(SpaCoObject)
}
