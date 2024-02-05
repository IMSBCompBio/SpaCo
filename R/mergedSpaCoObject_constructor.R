setClass("mergedSpaCoObject",
         representation(
           neighbours = "list",  # List to hold one or more matrices
           data = "list",        # List to hold one or more matrices
           coordinates = "list", # List to hold one or more data frames
           pixel_positions_list = "list", # List to hold one or more data frames
           data.dir = "character",
           slice = "character",
           spacs = "matrix",       # Assuming spacs is matrix for now
           projection = "list",    # List to hold one or more matrices
           GraphLaplacian = "matrix", # Assuming GraphLaplacian is sparse
           Lambdas = "numeric",
           nSpacs = "integer",
           meta.data = "list",     # List to hold one or more data frames
           denoised = "list"       # List to hold one or more data frames
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
mergedSpaCoObject <- function(neighbours, data, coordinates, pixel_positions_list, projection, meta.data, denoised) {
  new("SpaCoObject",
      neighbours = ensureList(neighbours, "dgCMatrix"),
      data = ensureList(data, "matrix"),
      coordinates = ensureList(coordinates, "data.frame"),
      pixel_positions_list = ensureList(pixel_positions_list, "data.frame"),
      data.dir = character(0),
      slice = character(0),
      spacs = new("matrix"),
      projection = ensureList(projection, "matrix"),
      GraphLaplacian = new("matrix"),
      Lambdas = numeric(0),
      nSpacs = integer(0),
      meta.data = ensureList(meta.data, "data.frame"),
      denoised = ensureList(denoised, "data.frame")
  )
}

mergeSpaCoObjects <- function(object1, object2, label1 = "dataset1", label2 = "dataset2") {
  # Check if objects are of the correct class
  if (!inherits(object1, "SpaCoObject") || !inherits(object2, "SpaCoObject")) {
    stop("Both objects must be of class 'SpaCoObject'")
  }

  # Append labels to sample names to ensure uniqueness
  object1 <- appendLabelsToNames(object1, label1)
  object2 <- appendLabelsToNames(object2, label2)

  # Proceed with merging as before
  mergedObject <- new("mergedSpaCoObject",
                      neighbours = list(object1@neighbours, object2@neighbours),
                      data = list(object1@data, object2@data),
                      coordinates = list(object1@coordinates, object2@coordinates),
                      pixel_positions_list = list(object1@pixel_positions_list, object2@pixel_positions_list),
                      data.dir = character(0),
                      slice = character(0),
                      spacs = new("matrix"),  # Assuming no merging needed for spacs
                      projection = list(object1@projection, object2@projection),
                      GraphLaplacian = new("matrix"),  # Assuming no merging needed for GraphLaplacian
                      Lambdas = numeric(0),
                      nSpacs = integer(0),
                      meta.data = mergeMetaData(object1@meta.data, object2@meta.data),
                      denoised = mergeLists(object1@denoised, object2@denoised)
  )

  return(mergedObject)
}

appendLabelsToNames <- function(object, label) {
  # Function to add label to the names of samples in various slots based on slot name rules

  slotsToLabel <- c("neighbours", "data", "coordinates", "pixel_positions_list", "projection", "GraphLaplacian", "meta.data", "denoised")

  # Define rules based on slot names
  rowsOnlySlots <- c("data","coordinates","pixel_positions_list","projection", "meta.data", "denoised")  # Slots to only modify row names
  colsOnlySlots <- c(NULL)  # Slots to only modify column names

  for (slot in slotsToLabel) {
    slotContent <- slot(object, slot)
    if (!is.null(slotContent) && ncol(slotContent) > 0) {
      # Modify row names if the slot is not in colsOnlySlots
      if (!slot %in% colsOnlySlots && !is.null(rownames(slotContent))) {
        rownames(slotContent) <- paste(rownames(slotContent), label, sep = "_")
      }
      # Modify column names if the slot is not in rowsOnlySlots
      if (!slot %in% rowsOnlySlots && !is.null(colnames(slotContent))) {
        colnames(slotContent) <- paste(colnames(slotContent), label, sep = "_")
      }
      slot(object, slot) <- slotContent
    }
  }
  return(object)
}

mergeMetaData <- function(df1, df2 ){
  # Step 1: Identify common and unique columns
  commonCols <- intersect(names(df1), names(df2))
  uniqueDf1 <- setdiff(names(df1), commonCols)
  uniqueDf2 <- setdiff(names(df2), commonCols)

  # Step 2: Merge data frames with all rows and appropriate columns
  # First, merge the common columns
  mergedCommon <- merge(df1[commonCols, , drop = FALSE], df2[commonCols, , drop = FALSE], by = "row.names", all = TRUE)

  # Now bind the unique columns back to the merged data frame
  finalMergedDf <- cbind(mergedCommon, df1[uniqueDf1], df2[uniqueDf2])

  # Step 3: Restore the original row names
  rownames(finalMergedDf) <- finalMergedDf$Row.names
  finalMergedDf <- finalMergedDf[, -1]  # Removing the added row names column
  return(mergedDf)

}




