filter.gene_custom <- function (object, thr, assay = NULL)
{
  if(is.null(assay)) assay <- Seurat::DefaultAssay(object)
  a = ifelse((thr >= 0 & thr < 1), thr * nrow(object), thr)
  object@assays[[assay]]@data = object@assays[[assay]]@scale.data[rowSums(object@assays[[assay]]@scale.data !=
                                                              0, na.rm = T) >= a, ]
  object
}
seurat.trans_custom <- function (object, counts = F, sps = F, assay = NULL)
{
  if(is.null(assay)) assay <- Seurat::DefaultAssay(object)
  if (counts == T) {
    c2 = t(as.matrix(object@assays[[assay]]@counts))
  }
  else {
    c2 = t(as.matrix(object@assays[[assay]]@scale.data))
  }
  colnames(c2) <- stringr::str_replace(colnames(c2), "-",
                                       "_")
  counts_frame = data.frame(cell = rownames(c2), c2)
  coor = data.frame(object@images[[Images(object)]]@coordinates,
                    cell = rownames(object@images[[Images(object)]]@coordinates))
  coor = coor[c("cell", "row", "col")]
  data = merge(coor, counts_frame, by = "cell")
  rownames(data) <- data$cell
  data = data[, -which(colnames(data) == "cell")]
  colnames(data) <- stringr::str_replace(colnames(data), "_",
                                         "-")
  if (sps == T) {
    dt = as(as.matrix(data), "sparseMatrix")
    rm(data)
  }
  else {
    dt = data
    rm(data)
  }
  dt
}
