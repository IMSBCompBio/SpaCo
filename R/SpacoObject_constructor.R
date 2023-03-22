# Define a custom class for the object
setClass("SpaCoObject",
         representation(neighbours = "matrix",
                        data = "matrix",
                        coordinates = "data.frame",
                        pixel_positions_list = "data.frame",
                        data.dir = "character",
                        slice = "character",
                        spacs = "matrix",
                        projection = "matrix",
                        GraphLaplacian="matrix",
                        Lambdas="vector",
                        nSpacs ="integer",
                        smoothed = "data.frame"))


#
#' Create a constructor function that creates an object of class SpaCoObject
#'
#' @param neighbours Binary matrix with weights describing if cells are to be considered neighbours or not depending on the defined distance.
#' @param data Matrix with normalized and scaled gene counts. Rows as cells and genes as columns
#' @param coordinates Matri with the cell coordinates on the slides. Rows and Columns in the 10x Visium case.
#'
#' @return Returns a SpaCoObject with the given slots filled
#' @export
#'
#'
SpaCoObject <- function(neighbours, data, coordinates) {
  new("SpaCoObject", neighbours = neighbours, data = data, coordinates = coordinates)
}
