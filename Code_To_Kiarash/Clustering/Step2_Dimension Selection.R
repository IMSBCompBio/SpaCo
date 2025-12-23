library(Seurat)
load(
  "/home/rstudio/SpaCo_paper_code/Clustering/Data/cluster_benchmark_objects.Rdata"
)
dim_reduction <- function(seurat_object,
                          resolution,
                          seed = NULL,
                          algorithm = 4,
                          dim_red = "spaco") {
  if (!is.null(seed))
    set.seed(seed)
  dims <- 1:min(10, ncol(seurat_object@reductions[[dim_red]]))
  seurat_object <-
    FindNeighbors(
      seurat_object,
      reduction = dim_red,
      dims = dims,
      verbose = FALSE
    )
  seurat_object <-
    FindClusters(
      seurat_object,
      resolution = resolution,
      algorithm = algorithm,
      random.seed = seed,
      verbose = FALSE
    )
  return(seurat_object)
}
get_n_clust <- function(seur, dim_red, res) {
  dim_red_obj <- dim_reduction(seur, res, NULL, 4, dim_red)
  cn <- rownames(seur@reductions[[dim_red]])
  length(levels(dim_red_obj@meta.data[cn, "seurat_clusters", drop = TRUE]))
}
n_clust_grid_helper <- function(res, seur) {
  sapply(c("pca", "spatialpca", "spaco", "spaco_svg"),
         function(x)
           get_n_clust(seur, x, res))
}
res_grid <- seq(0.1, 1.5, by = 0.1)
liver_reso_grid <- sapply(res_grid,
                          n_clust_grid_helper, seur = liver_res$seur_1)
colnames(liver_reso_grid) <- res_grid
posterior_reso_grid <- sapply(res_grid,
                              n_clust_grid_helper, seur = posterior_res$seur_1)
colnames(posterior_reso_grid) <- res_grid
anterior_reso_grid <- sapply(res_grid,
                             n_clust_grid_helper, seur = anterior_res$seur_1)
colnames(anterior_reso_grid) <- res_grid

resolutions <- matrix(
  c(1.5, 0.7, 0.8, 0.5,
    0.9, 0.8, 0.7, 0.5,
    1.2, 0.9, 1.1, 0.7),
  nrow = 3,
  ncol = 4,
  byrow = TRUE,
  dimnames = list(
    rows = c("liver", "anterior", "posterior"),
    cols = c("pca", "spaco", "spatialpca", "spaco_svg")
  )
)

save(
  liver_reso_grid,
  anterior_reso_grid,
  posterior_reso_grid,
  resolutions,
  file = "/tmp/cluster_benchmark_nclust_data.Rdata"
)

file.rename(
  from = "/tmp/cluster_benchmark_nclust_data.Rdata",
  to   = "/home/rstudio/SpaCo_paper_code/Clustering/Data/cluster_benchmark_nclust_data.Rdata"
)
