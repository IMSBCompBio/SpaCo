library(SPACO)
library(Seurat)
library(SPARK)
library(HEARTSVG)
source("/home/rstudio/SpaCo_paper_code/SVG/HEARTSVG Helper.R")
dir_name_in <- "/home/rstudio/SpaCo_paper_code/SVG/data/bootstrap_data/"
# dir_name_in <- "/home/rstudio/SpaCo_paper_code/SVG/data/additive_noise_data/"
dir_name <- "/home/rstudio/SpaCo_paper_code/SVG/data/"
fp_names <-
  readRDS(paste0(dir_name, "twin_names.RDS"))
tp_names <-
  readRDS(paste0(dir_name, "totwin.RDS"))
n_rads <- 5
n_bm <- 10
bm_settings <- rep(1:n_rads, each = n_bm)
spots_intersect <-
  readRDS("/home/rstudio/SpaCo_paper_code/SVG/data/intersect.RDS")
moran_i <- function(x, w) {
  # Calculate required components
  n <- length(x)
  S0 <- sum(w)
  xmean <- mean(x)
  deviations <- x - xmean
  S1 <- sum(deviations ^ 2)

  # Calculate spatial lag
  xlag <- w %*% deviations

  # Calculate Moran's I
  I <- (n / S0) * (t(deviations) %*% xlag / S1)

  return(I)
}
n_bootstraps <- 10
dir_name_in <-
  "/home/rstudio/SpaCo_paper_code/SVG/data/bootstrap_data/"
dir_name <- "/home/rstudio/SpaCo_paper_code/SVG/data/"

for (i in seq_along(neighbours_ls)) {
  for (j in 1:n_bootstraps) {
    set.seed((i - 1) * n_bootstraps + j)
    res <- list()
    brain <-
      readRDS(paste0(dir_name_in, "seur_", i, "_run_", j, ".RDS"))
    print("seurat read in")
    SpaCoObject_twin <-
      seurat_to_spaco(
        Seurat = brain,
        # assay = "shuffle",
        assay = "SCT_all",
        n_image = 1,
        layer = "scale.data"
      )
    SpaCoObject_twin@neighbours <-
      SpaCoObject_twin@neighbours[spots_intersect, spots_intersect]
    SpaCoObject_twin@data <-
      SpaCoObject_twin@data[spots_intersect, ]
    SpaCoObject_twin <-
      RunSCA(SpaCoObject_twin,
             set_nspacs = NULL)
    res[[1]] <- SpaCoObject_twin@nSpacs
    res[[2]] <- SVGTest(SpaCoObject_twin)

    coords <- scale(SpaCoObject_twin@coordinates[,])
    data <- brain@assays$SCT_all$scale.data
    data <- data[rowSums(data) != 0,]
    colnames(data) <- gsub("\\.", "-", colnames(data))
    coords <- coords[colnames(data), ]
    tryCatch({
      sparkX <-
        sparkx(as.matrix(data),
               coords,
               # numCores = 2,
               option = "mixture")
    },
    error = function(e) {
      # Handle the error if necessary
      print("Error with SparkX")
      DE_genes_SPARK[[i]] <-
        NULL  # Assign NULL or perform any other desired action
    })
    res[[3]] <- sparkX$res_mtest
    #Run HEARTSVG
    #Use helper function since HEARTSVG doesn't support different assays
    hsvgdat_tmp <- filter.gene_custom(brain, 0.01, "SCT_all")
    hsvgdat <- seurat.trans_custom(hsvgdat_tmp, assay = "SCT_all")
    hsvg_res <- heartsvg(hsvgdat)
    res[[4]] <- hsvg_res
    result_mat <-
      apply(SpaCoObject_twin@data, 2, function(column)
        moran_i(column, SpaCoObject_twin@neighbours))
    res[[5]] <- as.data.frame(result_mat)
    saveRDS(res, paste0(dir_name, "res/res_r_", i, "_run_", j, ".RDS"))
  }
}
