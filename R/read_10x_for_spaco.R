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
read_10x_for_spaco <- function(data_dir, slice, filename, variable_features_n = variable_features_n,spatial_file = "tissue_positions_list.csv",vars_to_regress = NULL) {
  require(Seurat)
  data <- Seurat::Load10X_Spatial(data.dir = data_dir, slice = slice, filename = filename, assay = "RNA", filter.matrix = TRUE )

 if(is.null(vars_to_regress)){
    data <- SCTransform(data, assay = "RNA", variable.features.n = variable_features_n)
  }
  else {
    for (i in 1:length(vars_to_regress)) {
      data <- Seurat::PercentageFeatureSet(data, pattern = vars_to_regress[i],col.name = make.names(paste0("percent",vars_to_regress[i])))
          }
    data <- SCTransform(data, assay = "RNA", variable.features.n = variable_features_n, vars.to.regress = make.names(paste0("percent",vars_to_regress)))
  }


  pixel_positions_list <- GetTissueCoordinates(data)

  data <- t(as.matrix(GetAssayData(object = data, assay = "SCT", slot = "scale.data")))

  tissue_positions_list <- read.csv(paste(slice, spatial_file, sep = "/"), col.names = c("barcode", "tissue", "row", "col", "imagerow", "imagecol"),row.names = 1,
                                    header = FALSE)

  tissue_positions_list <- tissue_positions_list[tissue_positions_list$tissue == 1, c("row", "col")]
  distm <- as.matrix(dist(tissue_positions_list, method = "euclidean", upper = TRUE))
  diag(distm) <- Inf
  neighboursindex <- distm <= 2

  if (any(colSums(neighboursindex) == 0)) {
    message("removing cells without any neighbours in defined distance")
    neighboursindex <- neighboursindex[colSums(neighboursindex) != 0, colSums(neighboursindex) != 0]
    data <- data[colnames(neighboursindex), ]
    tissue_positions_list <- tissue_positions_list[rownames(data),]
    pixel_positions_list <- pixel_positions_list[rownames(data),]
  }

  LociNames <- colnames(neighboursindex)
  data <- data[match(LociNames, rownames(data)),]

  SpaCoObject <- SpaCoObject(neighbours <-  neighboursindex, data <-  as.matrix(data), coordinates <- tissue_positions_list)
  slot(SpaCoObject, "pixel_positions_list") <- as.data.frame(pixel_positions_list)
  return(SpaCoObject)
  }





