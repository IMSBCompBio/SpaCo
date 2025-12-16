library(dplyr)
library(tidyr)
library(mclust)   # for adjustedRandIndex()
library(aricode)  # for NMI()
library(clue)     # for clue::cl_agreement(), if you ever want more measures
library(Seurat)
library(future)
library(future.apply)
load(
  "/home/rstudio/SpaCo_paper_code/Clustering/Data/cluster_benchmark_objects.Rdata"
)
load(
  "/home/rstudio/SpaCo_paper_code/Clustering/Data/cluster_benchmark_nclust_data.Rdata"
)
anterior_map = readRDS('/home/rstudio/SpaCo_paper_code/Clustering/Data/mapped_data_brain.RDS')
posterior_map = readRDS('/home/rstudio/SpaCo_paper_code/Clustering/Data/mapped_data_brain_posterior1.RDS')
liver_map = readRDS('/home/rstudio/SpaCo_paper_code/Clustering/Data/mapped_data_liver.RDS')

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
get_clusters <- function(seur, dim_red, res_c, cn) {
  dim_red_obj <- dim_reduction(seur, res_c[dim_red], NULL, 4, dim_red)
  # cn <- unlist(seur@misc)
  dim_red_obj@meta.data[cn, "seurat_clusters", drop = TRUE]
}
get_clust_wrapper <- function(seur, res_c, cn) {
  clusters <- sapply(c("pca", "spatialpca", "spaco", "spaco_svg"),
                     function(x)
                       get_clusters(seur, x, res_c, cn))
  colnames(clusters) <- c("pca", "spatialpca", "spaco", "spaco_svg")
  clusters
}
get_eval_metrics <- function(cluster_ls) {
  metrics <- lapply(list(NMI, adjustedRandIndex),
                    function(metric)
                      sapply(c("pca", "spatialpca", "spaco", "spaco_svg"),
                             function(dim_red)
                               metric(cluster_ls[[1]][, dim_red],
                                      cluster_ls[[2]][, dim_red])))
  df <- data.frame(matrix(unlist(metrics), ncol = 1))
  rownames(df) <- c(
    "NMI_pca",
    "NMI_spatialpca",
    "NMI_spaco",
    "NMI_spaco_svg",
    "ARI_pca",
    "ARI_spatialpca",
    "ARI_spaco",
    "ARI_spaco_svg"
  )
  df
}
compare_clusterings <- function(res_ls,
                                algorithm = 4,
                                res_c,
                                mapped_data,
                                n_iter = 1000) {
  cn1 <- mapped_data %>% filter(groups == 1) %>% rownames()
  md_sub <- mapped_data[match(cn1, mapped_data$map), , drop = TRUE]
  cn2 <- rownames(md_sub)

  # Setup parallel environment
  options(future.globals.maxSize = 4 * 1024 ^ 3)  # 4GB per worker if needed
  plan(multicore, workers = parallel::detectCores() - 1)

  # Efficient future loop: pass all data as args
  cluster_stats <-
    future_lapply(1:n_iter,
                  function(i, res_ls,
                           res_c,
                           cn1, cn2,
                           algorithm = 4) {
                    # get clusters
                    clusters <-
                      lapply(seq_along(res_ls),
                             function(i)
                               get_clust_wrapper(res_ls[[i]],
                                                 res_c = res_c,
                                                 cn = switch(i, cn1, cn2)))
                    # Metric calculation
                    metrics <- get_eval_metrics(clusters)
                    rbind(metrics, i = i)
                  },
                  # Explicitly pass all data objects as args
                  res_ls = res_ls, res_c, algorithm = algorithm,
                  cn1 = cn1, cn2 = cn2,
                  future.seed = TRUE)

  # Combine results
  return(do.call(rbind, cluster_stats))
}
anterior_cluster_res <-compare_clusterings(anterior_res, 4, resolutions["anterior",], anterior_map) # change to 1000 later
posterior_cluster_res <-compare_clusterings(posterior_res, 4, resolutions["posterior",], posterior_map)
liver_cluster_res <- compare_clusterings(liver_res, 4, resolutions["liver",], liver_map)
save(anterior_cluster_res,
     posterior_cluster_res,
     liver_cluster_res,
     file = "/home/rstudio/SpaCo_paper_code/Clustering/Data/cluster_results.Rdata")

#Example Cluster Figure
library(ggplot2)
set.seed(1)
seur_to_plot <- anterior_res$seur_1
cellnames_to_plot <-
  anterior_map %>% filter(groups == 1) %>% rownames()
clusters <-
  get_clust_wrapper(seur_to_plot, resolutions["anterior", ],
                    cellnames_to_plot)
plt_data <-
  clusters %>% as.data.frame() %>%
  cbind(anterior_map %>% filter(groups == 1) %>% select(x, y)) %>%
  pivot_longer(cols = all_of(colnames(clusters)),
               names_to = "Method",
               values_to = "Cluster") %>%
  filter(Method != "spaco_svg")
# table(plt_data$Method, plt_data$Cluster)
example_cluster_plot <-
  ggplot(plt_data, aes(x = y, y = x, color = Cluster)) +
  geom_point(size = 1.175, pch = 18) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1
    ),
    plot.margin = margin(-0.15,-0.15,-0.25,-0.25, "cm")
  ) +
  xlab("") + ylab("") +
  facet_wrap( ~ Method, ncol = 1)
library(svglite)
ggsave(
  "./Example Clusters.svg",
  plot = example_cluster_plot,
  width = 3,
  height = 8
)
