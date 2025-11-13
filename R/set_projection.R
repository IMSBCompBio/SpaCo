#' @param SpaCoObject Fitted object of class {\code{SpaCoObject} as returned from [RunSCA] on which to update the projection slot with new data.
#' @param data gene expression data matrix; p genes as columns, n loci as rows
#' @param neighbourhood neighborhood matrix of the new data. Must be of dimension n times n, where n is the number of columns of \code{data}
#' @return
#' Returns the SpaCoObject with the updated data, neighborhood matrix, and projection matrix.
#' @export
set_projection <- function(SpaCoObject, data, neighborhood) {
  requireNamespace("Matrix")
  SpaCoObject@data <- data
  GL <- computeGraphLaplacian(as(neighborhood, "dgCMatrix"))
  SpaCoObject@GraphLaplacian <- GL
  SpaCoObject@projection <-
    eigenMapMatMult(data, SpaCoObject@spacs)
    # as.matrix(Matrix::t(t(data) %*% GL) %*%
    #             SpaCoObject@spacs[, 1:SpaCoObject@nSpacs])
  SpaCoObject
}
