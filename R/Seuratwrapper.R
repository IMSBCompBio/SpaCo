#' Title
#'
#' @param Seurat
#' @param assay
#' @param n_image
#' @param slot
#'
#' @return
#' @export
#'
#' @examples
seurat_to_spaco <- function(Seurat, assay = "SCT", n_image = 1, slot = "scale.data") {
  data <- t(as.matrix(GetAssayData(object = Seurat, assay = DefaultAssay(Seurat), slot = "scale.data")))
   if (nrow(data) == 0 && ncol(data) == 0 ) {
     stop("Assay to transform from Seurat is empty.")
   }

  data <- t(as.matrix(GetAssayData(object = Seurat, assay = DefaultAssay(Seurat), slot = "scale.data")))
  tissue_positions_list <- as.data.frame(Seurat@images[[n_image]]@coordinates)

  tissue_positions_list <- tissue_positions_list[tissue_positions_list$tissue == 1, c("row", "col")]
  distm <- as.matrix(dist(tissue_positions_list, method = "euclidean", upper = TRUE))
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

#' Title
#'
#' @param SpaCoObject
#' @param Seurat
#'
#' @return
#' @export
#'
#' @examples
spaco_to_seurat <- function(SpaCoObject, Seurat ) {
  if (all(colnames(Seurat[[DefaultAssay(Seurat)]] )%in% rownames(SpaCoObject@data)) == FALSE) {
    stop("Cells without neighbours in defined distance found in Seurat object. Please subset cells first.")
  }
  message("copying projections into reduction slot spaco")
  Seurat[["spaco"]] <- CreateDimReducObject(embeddings = SpaCoObject@projection, key = "Spac_", assay = DefaultAssay(Seurat))

  return(Seurat)
}

#' Title
#'
#' @param SpaCoObject
#' @param Seurat
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

