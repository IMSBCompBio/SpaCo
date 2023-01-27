# Define a custom class for the object
setClass("SpaCoObject",
         representation(neighbours = "matrix",
                        data = "matrix",
                        coordinates= "data.frame",
                        data.dir = "character",
                        slice = "character",
                        spacs = "matrix",
                        projection = "matrix"))

# Create a constructor function that creates an object of class SpaCoobject
SpaCoObject <- function(neighbours, data, data.dir, slice, coordinates) {
  new("SpaCoObject", neighbours = neighbours, data = data, data.dir = data.dir, slice = slice, coordinates = coordinates)
}

# Create an object of class SpaCoobject
#my_obj <- SpaCoObject(matrix(1:9, nrow = 3, ncol = 3), data.frame(1:9, nrow = 3, ncol = 3), "path/to/data", "slice1")

