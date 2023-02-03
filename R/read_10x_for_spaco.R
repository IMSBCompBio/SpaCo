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
#' @import Seurat
read_10x_for_spaco <- function(data_dir, slice, filename, only_var = TRUE, variable_features_n = variable_features_n) {
  require(Seurat)
  data <- Load10X_Spatial(data.dir = data_dir, slice = slice, filename = filename, assay = "RNA", filter.matrix = TRUE)
  if (only_var) {
    data <- SCTransform(data, assay = "RNA", variable.features.n = variable_features_n)
  }
  tissue_positions_list <- read.csv(paste(slice, "tissue_positions_list.csv", sep = "/"), col.names = c("barcode", "tissue", "row", "col", "imagerow", "imagecol"),row.names = 1,
                                    header = FALSE)
  tissue_positions_list <- tissue_positions_list[tissue_positions_list$tissue == 1, c("row", "col")]
  distm <- as.matrix(dist(tissue_positions_list, method = "euclidean", upper = TRUE))
  diag(distm) <- Inf
  neighboursindex <- distm <= 2
  if (only_var) {
    data <- t(as.matrix(GetAssayData(object = data, assay = "SCT", slot = "scale.data")))
  } else {
    data <- t(as.matrix(GetAssayData(object = data, assay = "RNA", slot = "counts")))
  }
  if (any(colSums(neighboursindex) == 0)) {
    warning("removing cells without any neighbours in defined distance")
    #data <- data[colSums(neighboursindex) != 0 & rowSums(neighboursindex) != 0, ]
    neighboursindex <- neighboursindex[colSums(neighboursindex) != 0, colSums(neighboursindex) != 0]
    data <- data[colnames(neighboursindex), ]
    tissue_positions_list <- tissue_positions_list[rownames(data),]
  }

  LociNames <- colnames(neighboursindex)
  data <- data[match(LociNames, rownames(data)),]

  return(SpaCoObject(neighbours <-  neighboursindex, data <-  as.matrix(data), coordinates<- tissue_positions_list))
}




