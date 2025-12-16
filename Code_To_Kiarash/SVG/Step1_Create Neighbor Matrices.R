library(Seurat)
library(distances)

coords <-
  readRDS("/home/rstudio/SpaCo_paper_code/SVG/data/brain_coordinates.RDS")
write.table(coords,
            "/home/rstudio/SpaCo_paper_code/SVG/data/bootstrap_data/coords.csv")
dist_mat <- as.matrix(distances(coords))
diag(dist_mat) <- Inf

create_neighbor_mat <- function(dist_mat, r) {
  dist_mat <= r
}
rads <- c(2, 5, 10, 20, 40)
n_mats <-
  lapply(rads, function(r)
    create_neighbor_mat(dist_mat, r))
n_mats <- lapply(n_mats, function(mat) {
  colnames(mat) <- rownames(coords)
  rownames(mat) <- rownames(coords)
  mat
})
names(n_mats) <- rads
saveRDS(n_mats, file = "/home/rstudio/SpaCo_paper_code/SVG/data/neighmats.RDS")
