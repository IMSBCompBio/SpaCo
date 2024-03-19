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
                      neighbours = setNames(list(
                        object1@neighbours,
                        object2@neighbours,
                        MergeSparseMatrices(object1@neighbours, object2@neighbours)
                      ), c(label1, label2, "merged")),

                       data = setNames(list(
                        object1@data,
                        object2@data,
                        MergeDataMatrices(object1@data, object2@data)
                      ), c(label1, label2, "merged")),

                      coordinates = setNames(list(
                        object1@coordinates,
                        object2@coordinates
                      ), c(label1, label2)),

                      pixel_positions_list = setNames(list(
                        object1@pixel_positions_list,
                        object2@pixel_positions_list
                      ), c(label1, label2)),
                      data.dir = character(0),
                      slice = character(0),
                      spacs = new("matrix"),  # Assuming no merging needed for spacs
                      projection = setNames(list(
                        object1@projection,
                        object2@projection,
                        new("matrix")
                      ), c(label1, label2, "merged")),
                      GraphLaplacian = new("matrix"),  # Assuming no merging needed for GraphLaplacian
                      Lambdas = numeric(0),
                      nSpacs = integer(0),
                      meta.data = mergeMetaData(object1@meta.data, object2@meta.data),
                      denoised = new("list")
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
  return(finalMergedDf)

}

MergeSparseMatrices <- function(mat1, mat2) {
  require(Matrix)
  # Ensure both matrices are of class dgCMatrix
  if (!inherits(mat1, "dgCMatrix") || !inherits(mat2, "dgCMatrix")) {
    stop("Both matrices must be of class dgCMatrix.")
  }

  # Combine row and column names from both matrices
  row_names1 <- rownames(mat1)
  col_names1 <- colnames(mat1)
  row_names2 <- rownames(mat2)
  col_names2 <- colnames(mat2)

  all_row_names <- union(row_names1, row_names2)
  all_col_names <- union(col_names1, col_names2)

  # Create a new larger sparse matrix
  n_row <- length(all_row_names)
  n_col <- length(all_col_names)
  new_mat <- Matrix(0, n_row, n_col, sparse = TRUE)

  # Set row and column names
  rownames(new_mat) <- all_row_names
  colnames(new_mat) <- all_col_names

  # Map the rows and columns of the original matrices to the new matrix
  # and fill in the values
  if (!is.null(row_names1) && !is.null(col_names1)) {
    new_mat[row_names1, col_names1] <- mat1
  }
  if (!is.null(row_names2) && !is.null(col_names2)) {
    new_mat[row_names2, col_names2] <- mat2
  }

  return(new_mat)
}

#MergeDataMatricesAllGenes <- function(mat1, mat2) {
  # Ensure input are matrices
#  if (!is.matrix(mat1) || !is.matrix(mat2)) {
 #   stop("Both inputs must be matrices.")
#  }

#  mat1 <- scale(mat1, scale = TRUE)
#  mat2 <- scale(mat2, scale = TRUE)

  # Combine unique column names from both matrices
#  all_columns <- unique(c(colnames(mat1), colnames(mat2)))

  # Prepare new matrices with all columns, filled with zeros
#  mat1_expanded <- matrix(0, nrow = nrow(mat1), ncol = length(all_columns))
 # colnames(mat1_expanded) <- all_columns
 # rownames(mat1_expanded) <- rownames(mat1)

#  mat2_expanded <- matrix(0, nrow = nrow(mat2), ncol = length(all_columns))
#  colnames(mat2_expanded) <- all_columns
#  rownames(mat2_expanded) <- rownames(mat2)

  # Fill the expanded matrices with existing data
 # mat1_expanded[, colnames(mat1)] <- mat1
 # mat2_expanded[, colnames(mat2)] <- mat2

  # Combine the rows from both matrices
 # merged_mat <- rbind(mat1_expanded, mat2_expanded)

 # return(merged_mat)
#}

MergeDataMatrices <- function(mat1, mat2) {
  # Ensure input are matrices
  if (!is.matrix(mat1) || !is.matrix(mat2)) {
    stop("Both inputs must be matrices.")
  }

  # Scale the matrices
  mat1 <- scale(mat1, scale = TRUE)
  mat2 <- scale(mat2, scale = TRUE)

  # Find common columns in both matrices
  common_columns <- intersect(colnames(mat1), colnames(mat2))

  # Subset matrices to keep only common columns
  mat1_common <- mat1[, common_columns]
  mat2_common <- mat2[, common_columns]

  # Combine the rows from both matrices
  merged_mat <- rbind(mat1_common, mat2_common)

  return(merged_mat)
}

