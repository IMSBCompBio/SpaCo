Seurat_to_Spaco <- function(Seurat,assay="SCT"){
data <- t(as.matrix(GetAssayData(object = Seurat, assay = "assay", slot = "scale.data")))
tissue_positions_list <- SeuratObject::GetTissueCoordinates(Seurat)

tissue_positions_list <- tissue_positions_list[tissue_positions_list$tissue == 1, c("row", "col")]
distm <- as.matrix(dist(tissue_positions_list, method = "euclidean", upper = TRUE))
diag(distm) <- Inf
neighboursindex <- distm <= 2

if (any(colSums(neighboursindex) == 0)) {
  warning("removing cells without any neighbours in defined distance")
  #data <- data[colSums(neighboursindex) != 0 & rowSums(neighboursindex) != 0, ]
  neighboursindex <- neighboursindex[colSums(neighboursindex) != 0, colSums(neighboursindex) != 0]
  data <- data[colnames(neighboursindex), ]
  tissue_positions_list <- tissue_positions_list[rownames(data),]
}

LociNames <- colnames(neighboursindex)
data <- data[match(LociNames, rownames(data)),]

return(SpaCoObject(neighbours <-  neighboursindex, data <-  as.matrix(data), data_dir <- data_dir, slice <- slice, coordinates<- tissue_positions_list))
}

Spaco_to_Seurat <- function(SpaCoObject,weiteres){

}
