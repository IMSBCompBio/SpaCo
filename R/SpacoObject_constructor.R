# Define a custom class for the object

setClass("SpaCoObject",
         representation(neighbours = "dgCMatrix",  # Assuming neighbours is sparse
                        data = "matrix",        # Assuming data can be sparse
                        data_B = "matrix",
                        coordinates = "data.frame",
                        pixel_positions_list = "data.frame",
                        data.dir = "character",
                        slice = "character",
                        spacs = "matrix",       # Assuming spacs is sparse
                        projection = "matrix",  # Assuming projection is sparse
                        projection_B = "matrix",
                        GraphLaplacian = "matrix", # Assuming GraphLaplacian is sparse
                        GraphLaplacian_B = "matrix",
                        Lambdas = "numeric",
                        nSpacs = "integer",
                        meta.data = "data.frame",
                        denoised = "data.frame"
                        ))

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
# Create a constructor function that creates an object of class SpaCoObject
# Define the custom class for the object
SpaCoObject <- function(neighbours, data, coordinates, pixel_positions_list) {
  new("SpaCoObject",
      neighbours = as(neighbours, "dgCMatrix"),
      data = data,
      coordinates = coordinates,
      pixel_positions_list = new("data.frame"),
      data.dir = character(0),
      slice = character(0),
      spacs = new("matrix"),
      projection = new("matrix"),
      projection_B = new("matrix"),
      GraphLaplacian = new("matrix"),
      GraphLaplacian_B = new("matrix"),
      Lambdas = numeric(0),
      nSpacs = integer(0),
      meta.data = new("data.frame"),
      denoised = new("data.frame")
  )
}






