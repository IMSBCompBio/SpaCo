#' Read in 10x Visium spatial transcriptomics data
#'
#' @param data_dir Directory containing the H5 file specified by file name and the image data in a sub directory called spatial
#' @param slice Name for the stored image of the tissue slice
#' @param only_var Logical True if you want to keep only the most variable features in the SPACO object. Default True.
#' @param variable_features_n Number of most variable features to keep.
#'
#' @return
#' @export
#'
#' @examples
read_10x_for_spaco <- function(data_dir,
                                slice,
                                filename,
                                only_var = TRUE,
                                variable_features_n = variable_features_n) {
  require(Seurat)
  require(gmodels)
  data <- Load10X_Spatial(data.dir = data_dir,
                          slice  = slice,
                          filename = filename,
                          assay = "RNA",
                          filter.matrix = TRUE,
                          to.upper = FALSE)
  if (only_var) {
    message("performing SCTransform and variable feature selection")

    data <- SCTransform(data, assay = "RNA",
                        verbose = TRUE,
                        variable_features_n = variable_features_n)

  }
  tissue_positions_list <- read.csv(paste(slice, "tissue_positions_list.csv", sep = "/"),
                                    col.names = c("barcode", "tissue", "row", "col", "imagerow", "imagecol"),
                                    header = FALSE)

  rownames(tissue_positions_list) <- tissue_positions_list$barcode
  tissue_positions_list <- tissue_positions_list[tissue_positions_list$tissue == 1, ]
  coords <-  tissue_positions_list[, c("row", "col")]


  distm <- dist(coords, method = "euclidean", upper = TRUE)
  distm <- as.matrix(distm)
  diag(distm) <- Inf
  neighboursindex <- matrix(nrow = length(rownames(distm)),
                            ncol = length(colnames(distm)))
  rownames(neighboursindex) <- rownames(distm)
  colnames(neighboursindex) <- colnames(distm)


  for (i in seq_len(ncol(distm))){
    neighboursindex[, i] <- ifelse(distm[, i] <= 2, yes = 1, no = 0)
  }


  if (only_var) {
    data <- t(as.matrix(GetAssayData(object = data,
                                     assay = "SCT",
                                     slot = "scale.data")))
  } else {
    data <- t(as.matrix(GetAssayData(object = data,
                                     assay = "RNA",
                                     slot = "counts")))
  }


  if (any(colSums(neighboursindex) == 0) == "TRUE") {
    warning("removing spots without direct neighbours")
    data <- data[!rownames(data) %in% names(which(colSums(neighboursindex) == 0)), ]
    coords <-  coords[!rownames(coords) %in% names(which(colSums(neighboursindex) == 0)), ]
    neighboursindex <- neighboursindex[-which(colSums(neighboursindex) == 0),
                                       -which(colSums(neighboursindex) == 0)]

  }



  poslist <- as.data.frame(which(neighboursindex == 1, arr.ind = TRUE))

  edges_DF <- data.frame(pair1 = colnames(neighboursindex)[poslist$col],
                     pair2 = rownames(neighboursindex)[poslist$row])

  return(list(neighboursindex <-  neighboursindex, edges_DF <-  edges_DF, data <-  as.matrix(data), plotLayout <- coords, data_dir <- data_dir, slice <- slice))
}
