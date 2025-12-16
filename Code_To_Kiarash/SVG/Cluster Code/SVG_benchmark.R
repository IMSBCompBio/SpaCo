rm(list = ls())
library(batchtools)

setwd("/home/koehler/SPACO")
cluster.functions <-
  makeClusterFunctionsSlurm(
    template = "batchtools.slurm.tmpl",
    array.jobs = TRUE,
    scheduler.latency = 1,
    fs.latency = 65
  )
set.seed(1)

n_rads <- 5
n_bm <- 10
bm_settings <- rep(1:n_rads, each = n_bm)
test_bm <- function(i) {
  setwd("/home/koehler/SPACO")
  set.seed(i)
  library(SPACO)
  library(Seurat)
  library(SPARK)
  library(HEARTSVG)
  source("./HEARTSVG Helper.R")
  dir_name_in <- "./data/bootstrap_data/"
  # dir_name_in <- "./data/additive_noise_data/"
  dir_name <- "./data/"
  fp_names <-
    readRDS(paste0(dir_name, "twin_names.RDS"))
  tp_names <-
    readRDS(paste0(dir_name, "totwin.RDS"))
  n_rads <- 5
  n_bm <- 10
  bm_settings <- rep(1:n_rads, each = n_bm)
  rad_i <- (i - 1) %/% 10  + 1
  j <- (i - 1) %% 10 + 1
  spots_intersect <-
    readRDS("./data/intersect.RDS")
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
  res <- list()
  print("Data read in")
  brain <-
    readRDS(paste0(dir_name_in, "seur_",rad_i, "_run_", j, ".RDS"))
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
    SpaCoObject_twin@data[spots_intersect,]
  print("spc created")
  SpaCoObject_twin <-
    RunSCA(SpaCoObject_twin,
           set_nspacs = NULL)
  print("RunSCA done")
  res[[1]] <- SpaCoObject_twin@nSpacs

  res[[2]] <- SVGTest(SpaCoObject_twin)

  coords <- scale(SpaCoObject_twin@coordinates[, ])
  data <- brain@assays$SCT_all$scale.data
  data <- data[rowSums(data) != 0, ]
  colnames(data) <- gsub("\\.", "-", colnames(data))
  coords <- coords[colnames(data),]
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
  print("running HEARTSVG")
  #Use helper function since HEARTSVG doesn't support different assays
  hsvgdat_tmp <- filter.gene_custom(brain, 0.01, "SCT_all")
  hsvgdat <- seurat.trans_custom(hsvgdat_tmp, assay = "SCT_all")
  hsvg_res <- heartsvg(hsvgdat)
  res[[4]] <- hsvg_res
  print("computinmorans")
  result_mat <-
    apply(SpaCoObject_twin@data, 2, function(column)
      moran_i(column, SpaCoObject_twin@neighbours))
  res[[5]] <- as.data.frame(result_mat)
  saveRDS(res, paste0(dir_name, "res/res_r_", rad_i, "_run_", j, ".RDS"))
}
reg <- makeRegistry(file.dir = paste0("/home/koehler/SPACO/results"))

ids <-
  batchMap(fun = test_bm,
           args = list(i = 1:length(bm_settings)),
           #args = list(i = 1:1),
	 reg = reg)

submitJobs(
  ids,
  res = list(
    memory = 100000,
    partition = "batch",
    walltime = 50000
  ),
  reg = reg
)

while (!waitForJobs(ids = ids)) {
  Sys.sleep(1)
}

# results <- reduceResultsList(reg = reg)
# save(results, file = "test.rds")
