#' Wrapper to transform existing Seurat object into an SpaCoObject.
#'
#' @param Seurat Seurat object to export
#' @param assay Assay to export from the Seurat object. Default is SCT assay.
#' @param n_image Number of the image to export from Seurat object. Only relevant if Seurat object contains multiple images. Default is 1.
#' @param slot Which slot to export data from. Default is scale.data
#'
#' @return
#' @export
#'
#' @examples
seurat_to_spaco <- function(Seurat, assay = "SCT", n_image = 1, slot = "scale.data") {
  data <- t(as.matrix(Seurat::GetAssayData(object = Seurat, assay = Seurat::DefaultAssay(Seurat), slot = "scale.data")))
   if (nrow(data) == 0 && ncol(data) == 0 ) {
     stop("Assay to transform from Seurat is empty.")
   }

  tissue_positions_list <- as.data.frame(Seurat@images[[n_image]]@coordinates)

  tissue_positions_list <- tissue_positions_list[tissue_positions_list$tissue == 1, c("row", "col")]
  distm <- as.matrix(stats::dist(tissue_positions_list, method = "euclidean", upper = TRUE))
  diag(distm) <- Inf
  neighboursindex <- distm <= 2

if (any(colSums(neighboursindex) == 0)) {
  warning("removing cells without any neighbours in defined distance")
  neighboursindex <- neighboursindex[colSums(neighboursindex) != 0, colSums(neighboursindex) != 0]
  data <- data[colnames(neighboursindex), ]
  tissue_positions_list <- tissue_positions_list[rownames(data), ]
}

LociNames <- colnames(neighboursindex)
data <- data[match(LociNames, rownames(data)), ]

return(SpaCoObject(neighbours <-  neighboursindex, data <-  as.matrix(data), coordinates <- tissue_positions_list))
}

#' Export computed spatial components to existing Seurat object.
#'
#' @param SpaCoObject SpaCoObject to export spatial components from.
#' @param Seurat Seurat object to add spatial components to.
#'
#' @return
#' @export
#'
#' @examples
spacs_to_seurat <- function(SpaCoObject, Seurat ) {
  if (all(colnames(Seurat[[Seurat::DefaultAssay(Seurat)]] )%in% rownames(SpaCoObject@data)) == FALSE) {
    stop("Cells without neighbours in defined distance found in Seurat object. Please subset cells first.")
  }
  message("copying projections into reduction slot spaco")
  Seurat[["spaco"]] <- Seurat::CreateDimReducObject(embeddings = SpaCoObject@projection, key = "Spac_", assay = Seurat::DefaultAssay(Seurat))

  return(Seurat)
}

#' Filtering function to remove cells without neighbours in defined distance from existing Seurat object to be conformable with existing SpaCoObject.
#'
#' @param SpaCoObject SpaCoObject to integrate into Seurat object.
#' @param Seurat Seurat object to be filtered.
#'
#' @return
#' @export
#'
#' @examples
subset_non_neighbour_cells <- function(SpaCoObject, Seurat) {
  require(Seurat)

  Seurat <- subset(Seurat, cells = intersect(colnames(Seurat),rownames(SpaCoObject@data)))

  return(Seurat)
}
