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
  library(Matrix)
  library(Rcpp)
  library(SPACO)
  library(Seurat)
  library(SPARK)
  library(HEARTSVG)
  source("./HEARTSVG Helper.R")
  brain <-
    readRDS("./data/brain.RDS")
  neighbours_ls <-
    readRDS("./data/neighmats.RDS")
  orig_data <-
    readRDS("./data/twin_data.RDS")
  orig_names <- rownames(orig_data)

  spots_intersect <-
    intersect(colnames(orig_data), colnames(neighbours_ls[[1]]))
  saveRDS(spots_intersect,
          "./data/intersect.RDS")
  orig_data <- orig_data[rowSums(orig_data) != 0, spots_intersect]
  neighbours_ls <-
    lapply(neighbours_ls, function(mat)
      mat[spots_intersect, spots_intersect])
  tp_names <-
    readRDS("./data/totwin.RDS")
  rads <- names(neighbours_ls)
  fp_names <-
    readRDS("./data/twin_names.RDS")
  tp_names <-
    readRDS("./data/totwin.RDS")

  coverage_adjusted_bootstrap <-
    function(orig_data,
             neighborhood_matrix) {
      resample_rows_fixed_neighborhood <-
        function(data_matrix, neighborhood_matrix) {
          num_cols <- ncol(data_matrix)
          num_rows <- nrow(data_matrix)
          new_data <- matrix(data = NA,
                             nrow = num_rows,
                             ncol = num_cols)
          colnames(new_data) <- colnames(data_matrix)
          rownames(new_data) <- rownames(data_matrix)
          new_spots <- sapply(1:num_rows,
                              function(spot)
                                sample(which(neighborhood_matrix[spot,] == 1),
                                       num_cols,
                                       replace = TRUE))
          new_data <-
            sapply(1:num_cols, function(gene_idx)
              data_matrix[new_spots[gene_idx, ], gene_idx])
          colnames(new_data) <- colnames(data_matrix)
          rownames(new_data) <- rownames(data_matrix)
          return(new_data)
        }
      diag(neighborhood_matrix) <- 1
      coverage_vector <- colSums(orig_data)
      orig_data_abund <-
        as.matrix(sweep(orig_data, 2, coverage_vector, "/"))
      tp_in_data <-
        tp_names[which(tp_names %in% rownames(orig_data_abund))]
      svg_sample_matrix_p <-
        t(resample_rows_fixed_neighborhood(t(orig_data_abund[tp_in_data, ]), neighborhood_matrix))
      orig_data_abund[tp_in_data,] <- svg_sample_matrix_p
      random_coverage <- colSums(orig_data_abund)
      counts_sample <-
        as.matrix(sweep(orig_data_abund, 2, coverage_vector / random_coverage, "*"))
      counts_sample <- round(counts_sample)
      return(counts_sample)
    }
  dir_name <- "./data/bootstrap_data/"
  rad_i <- (i - 1) %/% 10  + 1
  j <- (i - 1) %% 10 + 1
  bootstrapped_data <-
    coverage_adjusted_bootstrap(orig_data = orig_data,
                                neighbours_ls[[rad_i]])
  bootstrapped_data <-
    bootstrapped_data[rowSums(bootstrapped_data) != 0,]
  adt_assay <- CreateAssayObject(counts = bootstrapped_data)
  brain[["shuffle"]] <- adt_assay
  brain <-
    PercentageFeatureSet(brain,
                         pattern = "^mt-" ,
                         col.name = "percent.mt",
                         assay = "shuffle")
  brain <-
    PercentageFeatureSet(brain,
                         pattern = "^Hbb-" ,
                         col.name = "percent.hbb",
                         assay = "shuffle")
  brain <-
    SCTransform(brain,
                assay = "shuffle",
                verbose = F,
                variable.features.n = 3000,
                seed.use = i)
  brain <-
    SCTransform(
      brain,
      assay = "shuffle",
      verbose = F,
      seed.use = i,
      residual.features =
        union(
          rownames(brain@assays$SCT$scale.data),
          c(fp_names, tp_names)
        ),
      new.assay.name = "SCT_all"
    )
  write.table(
    brain@assays$shuffle$data,
    file = paste0(dir_name, "genemat_r_", rad_i, "_run_", j, ".csv")
  )
  brain@assays$Spatial <- NULL
  brain@assays$shuffle <- NULL
  brain@assays$SCT <- NULL
  saveRDS(brain,
          paste0(dir_name, "seur_", rad_i, "_run_", j, ".RDS"))
}
reg <-
  makeRegistry(file.dir = paste0("/home/koehler/SPACO/results"))

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
