#' Wrapper to transform existing Seurat object into an SpaCoObject.
#'
#' @param Seurat Seurat object to export
#' @param assay Assay to export from the Seurat object. Default is SCT assay.
#' @param n_image Number of the image to export from Seurat object. Only relevant if Seurat object contains multiple images. Default is 1.
#' @param slot Which slot to export data from. Default is scale.data.
#'
#' @return Returns a SpaCoObject with all slots filled necessary to perform the spatial component analysis.
#' @export
#'
#'
seurat_to_spaco <- function(Seurat, assay = "SCT", n_image = 1, slot = "scale.data") {
  data <- t(as.matrix(Seurat::GetAssayData(object = Seurat, assay = Seurat::DefaultAssay(Seurat), slot = slot)))
   if (nrow(data) == 0 && ncol(data) == 0 ) {
     stop("Assay to transform from Seurat is empty.")
   }
  pixel_positions_list <- GetTissueCoordinates(Seurat)
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
  pixel_positions_list <- pixel_positions_list[rownames(data),]
}

LociNames <- colnames(neighboursindex)
data <- data[match(LociNames, rownames(data)), ]

SpaCoObject <- SpaCoObject(neighbours <-  neighboursindex, data <-  as.matrix(data), coordinates <- tissue_positions_list)
slot(SpaCoObject, "pixel_positions_list") <- as.data.frame(pixel_positions_list)
return(SpaCoObject)
}

#' transfer computed spatial components to existing Seurat object.
#'
#' @param SpaCoObject SpaCoObject to export spatial components from.
#' @param Seurat Seurat object to add spatial components to.
#'
#' @return Returns a Seurat Object with the spatial components projections in the dimensional reduction slot.
#' @export
#'
#'
spacs_to_seurat <- function(SpaCoObject, Seurat, nSpacs = SpaCoObject@nSpacs ) {
  if (all(colnames(Seurat[[Seurat::DefaultAssay(Seurat)]] )%in% rownames(SpaCoObject@data)) == FALSE) {
    stop("Cells without neighbours in defined distance found in Seurat object. Please subset cells first.")
  }
  message("copying significant projections into reduction slot spaco")
  Seurat[["spaco"]] <- Seurat::CreateDimReducObject(embeddings = SpaCoObject@projection[,1:nSpacs], key = "Spac_", assay = Seurat::DefaultAssay(Seurat))

  return(Seurat)
}

#' Filtering function to remove cells without neighbours in defined distance from existing Seurat object to be conformable with existing SpaCoObject.
#'
#' @param SpaCoObject SpaCoObject to integrate into Seurat object.
#' @param Seurat Seurat object to be filtered.
#'
#' @return Returns a Seurat object with cells filtered to match SpaCoObject.
#' @export
#'
#'
subset_non_neighbour_cells <- function(SpaCoObject, Seurat) {
  require(Seurat)

  Seurat <- subset(Seurat, cells = intersect(colnames(Seurat),rownames(SpaCoObject@data)))

  return(Seurat)
}

#' create_SpaCoObject_from_KNN
#'
#' @param Seurat Seurat object to export kNN-graph from.
#' @param n Number of neighbors to consider.
#' @return Returns a SPaCoObject with the SCT data and the kNN-graph as neighborhood matrix.
#' @export
#'
#'
create_SpaCoObject_from_KNN <- function(Seurat, n=10){
  neighbourindexmatrix <- matrix(data=0,nrow = length(colnames(Seurat)),ncol=length(colnames(Seurat)))
  rownames(neighbourindexmatrix) <- colnames(Seurat)
  colnames(neighbourindexmatrix) <- colnames(Seurat)

  for (i in 1:length(colnames(Seurat))){
    neighbourindexmatrix[colnames(Seurat)[i],TopNeighbors(Seurat[["SCT.nn"]],cell=colnames(Seurat)[i],n=n)] <- 1
    neighbourindexmatrix[TopNeighbors(Seurat[["SCT.nn"]],cell=colnames(Seurat)[i],n=n),colnames(Seurat)[i]] <- 1
  }
  diag(neighbourindexmatrix) <- 0


data_to_spac <- t(as.matrix(Seurat::GetAssayData(object = Seurat, assay = "SCT", slot = "scale.data")))
LociNames <- colnames(neighbourindexmatrix)
data_to_spac <- data_to_spac[match(LociNames, rownames(data_to_spac)),]
tissue_positions_list <- as.data.frame(Embeddings(Seurat, reduction = "umap"))
pixel_positions_list <- as.data.frame(Embeddings(Seurat, reduction = "umap"))
pixel_positions_list <- pixel_positions_list[match(LociNames, rownames(pixel_positions_list)),]

SpaCoObject <-SpaCoObject(neighbours=neighbourindexmatrix, data=data_to_spac, coordinates =pixel_positions_list)
tissue_positions_list <- as.data.frame(Embeddings(Seurat, reduction = "umap"))
tissue_positions_list <- tissue_positions_list[rownames(data_to_spac),]
colnames(tissue_positions_list) <- c("imagecol","imagerow")
slot(SpaCoObject, "pixel_positions_list") <- as.data.frame(tissue_positions_list)
return(SpaCoObject)
}






