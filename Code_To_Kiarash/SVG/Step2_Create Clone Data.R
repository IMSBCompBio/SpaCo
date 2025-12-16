library(dplyr)
library(Seurat)
set.seed(1)
brain <-
  readRDS("./data/brain.RDS")
orig_data <- GetAssayData(brain, assay = "Spatial")
neighbours_ls <-
  readRDS("./data/neighmats.RDS")
spots_intersect <-
  intersect(colnames(orig_data), colnames(neighbours_ls[[1]]))

orig_data <- orig_data[, spots_intersect]
svg_names <-
  readRDS("./data/totwin.RDS")
twin_names <-
  readRDS("./data/twin_names.RDS")
coverage_vector <- colSums(orig_data)
orig_data_abund <-
  as.matrix(sweep(orig_data, 2, coverage_vector, "/"))
data_to_shuffle <-
  orig_data_abund[svg_names,]
new_spots <- t(replicate(length(svg_names),
                         sample(1:ncol(data_to_shuffle), )))
new_data <-
  t(sapply(1:nrow(data_to_shuffle), function(gene_idx)
    data_to_shuffle[gene_idx, new_spots[gene_idx, ]]))
rownames(new_data) <- twin_names
orig_data_abund_twins <-
  rbind(orig_data_abund, new_data)
orig_data_twins <-
  as.matrix(sweep(
    orig_data_abund_twins,
    2,
    coverage_vector / colSums(orig_data_abund_twins),
    "*"
  ))
#probabilistic rounding
orig_data_twins <-
  floor(orig_data_twins) + (runif(length(orig_data_twins)) < (orig_data_twins - floor(orig_data_twins)))
saveRDS(orig_data_twins,
        "./data/twin_data.RDS")
