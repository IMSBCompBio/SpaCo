# Define a custom class for the object
setClass("SpaCoObject",
         representation(neighbours = "matrix",
                        data = "matrix",
                        coordinates= "data.frame",
                        data.dir = "character",
                        slice = "character",
                        spacs = "matrix",
                        projection = "matrix",
                        GraphLaplacian="matrix",
                        Lambdas="vector"))

# Create a constructor function that creates an object of class SpaCoobject
#' Title
#'
#' @param neighbours
#' @param data
#' @param coordinates
#'
#' @return
#' @export
#'
#' @examples
SpaCoObject <- function(neighbours, data, coordinates) {
  new("SpaCoObject", neighbours = neighbours, data = data, coordinates = coordinates)
}
