#' Read in 10x Visium spatial transcriptomics data
#'
#' @param data_dir Directory containing the H5 file specified by file name and the image data in a sub directory called spatial
#' @param slice Name for the stored image of the tissue slice
#' @param variable_features_n Number of most variable features to keep.
#'
#' @return Retuns a ready to run SPaCoObject.
#' @export
#'
#'
#' @import Seurat
read_10x_for_spaco <- function(data_dir,
                               slice,
                               filename,
                               variable_features_n = variable_features_n,
                               spatial_file = spatial_file,
                               vars_to_regress = NULL) {
  require(Seurat)

  # Load data
  data <- Load10X_Spatial(data.dir = data_dir,
                          filename = filename,
                          assay = "RNA",
                          filter.matrix = TRUE)

  # Normalization based on vars_to_regress
  if (is.null(vars_to_regress)) {
    data <- SCTransform(data, assay = "RNA", variable.features.n = variable_features_n)
  } else {
    for (var in vars_to_regress) {
      percent_col_name <- make.names(paste0("percent", var))
      data <-PercentageFeatureSet(data, pattern = var, col.name = percent_col_name)
    }
    regress_vars <- make.names(paste0("percent", vars_to_regress))
    data <- SCTransform(data, assay = "RNA", variable.features.n = variable_features_n, vars.to.regress = regress_vars)
  }

  # Coordinate processing
  pixel_positions_list <- GetTissueCoordinates(data)
  meta <- as.data.frame(data@meta.data)
  data_matrix <- t(as.matrix(GetAssayData(object = data, assay = "SCT", slot = "scale.data")))

  # Read and process tissue positions
  tissue_positions_list <- read.csv(paste(slice, spatial_file, sep = "/"),
                                    col.names = c("barcode", "tissue", "row", "col", "imagerow", "imagecol"),
                                    row.names = 1,
                                    header = TRUE)
  tissue_positions_list <- tissue_positions_list[tissue_positions_list$tissue == 1, c("row", "col")]

  # Distance matrix and neighbor index
  distm <- as.matrix(dist(tissue_positions_list, method = "euclidean", upper = TRUE))
  diag(distm) <- Inf
  neighbours_index <- distm <= 2

  # Handling cells without neighbors
  if (any(colSums(neighbours_index) == 0)) {
    message("Removing cells without any neighbours in defined distance")
    valid_cells <- colSums(neighbours_index) != 0
    neighbours_index <- neighbours_index[valid_cells, valid_cells]
    data_matrix <- data_matrix[colnames(neighbours_index), ]
    tissue_positions_list <- tissue_positions_list[rownames(data_matrix), ]
    pixel_positions_list <- pixel_positions_list[rownames(data_matrix), ]
  }

  # Final data adjustments
  meta <- meta[rownames(data_matrix), ]
  loci_names <- colnames(neighbours_index)
  data_matrix <- data_matrix[match(loci_names, rownames(data_matrix)), ]
  pixel_positions_list <- pixel_positions_list[match(loci_names, rownames(pixel_positions_list)), ]

  # Create SpaCoObject
  SpaCoObject <- SpaCoObject(neighbours = neighbours_index,
                             data = as.matrix(data_matrix),
                             coordinates = tissue_positions_list,
                             pixel_positions_list = pixel_positions_list)
  SpaCoObject@meta.data <- meta
  SpaCoObject@pixel_positions_list <- as.data.frame(pixel_positions_list)

  return(SpaCoObject)
}
