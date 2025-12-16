library(devtools)
library(Seurat)
library(SPACO)
# mclust for adjustedRandIndex (ARI)
library(mclust)            # >= 6.0.0
# aricode (or infotheo, etc.) for Normalized Mutual Information (NMI)
library(aricode)           # contains NMI() function
# tidyverse for data‐manipulation and plotting
library(tidyverse)
# SpatialPCA for benchmarking
# devtools::install_github(SpatialPCA)
library(SpatialPCA)
library(glmGamPoi)
## reading in the Anterior brain data
anterior_map = readRDS('/home/rstudio/SpaCo_paper_code/Clustering/Data/mapped_data_brain.RDS')
anterior_slice_data = readRDS('/home/rstudio/SpaCo_paper_code/Clustering/Data/brain.RDS')
anterior_slice_seur <- UpdateSeuratObject(anterior_slice_data)

### Read in liver data
liver_map = readRDS('/home/rstudio/SpaCo_paper_code/Clustering/Data/mapped_data_liver.RDS')
liver_data = readRDS('/home/rstudio/SpaCo_paper_code/Clustering/Data/seurat_data_clust_liver.RDS')
liver_seur <- UpdateSeuratObject(liver_data)

### Reading in posterior data
posterior_map = readRDS('/home/rstudio/SpaCo_paper_code/Clustering/Data/mapped_data_brain_posterior1.RDS')
posterior_slice_data = readRDS('/home/rstudio/SpaCo_paper_code/Clustering/Data/seurat_data_clust_brain_posterior1.RDS')
posterior_slice_seur <- UpdateSeuratObject(posterior_slice_data)

set.seed(1)

seur_group_pipeline <-
  function(seur,
           map_data,
           group,
           vars_to_regress = c("^mt-", "^Hbb-"),
           variable_features_n = 3000,
           assay = "Spatial") {
    cn <- map_data %>% filter(groups == group) %>% rownames()
    seur_sub <-
      subset(seur, cells = cn)
    for (var_pattern in vars_to_regress) {
      col_name <- make.names(paste0("percent", var_pattern))
      seur_sub <- Seurat::PercentageFeatureSet(seur_sub,
                                               pattern = var_pattern,
                                               col.name = col_name)
    }
    seur_sub <- SCTransform(
      seur_sub,
      assay = assay,
      variable.features.n = variable_features_n,
      vars.to.regress = make.names(paste0("percent", vars_to_regress))
    )
    seur_sub
  }

seurat_to_spaco_spatialpca <-
  function(seur_spatialpca_ls,
           assay = "SCT",
           layer = "scale.data") {
    data <-
      t(as.matrix(
        Seurat::GetAssayData(
          object = seur_spatialpca_ls$seur,
          assay = Seurat::DefaultAssay(seur_spatialpca_ls$seur),
          layer = layer
        )
      ))
    if (nrow(data) == 0 && ncol(data) == 0) {
      stop("Assay to transform from Seurat is empty.")
    }
    pixel_positions_list <- GetTissueCoordinates(seur_spatialpca_ls$seur)
    tissue_positions_list <-
      as.data.frame(seur_spatialpca_ls$seur@images[[1]]@coordinates)
    SpaCoObject <-
      SpaCoObject(
        neighbours <-
          seur_spatialpca_ls$kernel,
        data <-  as.matrix(data),
        coordinates <- tissue_positions_list
      )
    SpaCoObject
  }

spaco_svg_projs <- function(spaco_obj, genes) {
  gene_idx <- which(colnames(spaco_obj@data) %in% genes)
  spaco_obj@projection <- spaco_obj@data[,gene_idx] %*% spaco_obj@spacs[gene_idx,]
  spaco_obj
}

spacs_to_seurat_svg <-
  function(SpaCoObject, Seurat, nSpacs = SpaCoObject@nSpacs) {
    if (all(colnames(Seurat[[Seurat::DefaultAssay(Seurat)]]) %in% rownames(SpaCoObject@data)) == FALSE) {
      stop(
        "Cells without neighbours in defined distance found in Seurat object. Please subset cells first."
      )
    }
    message("copying significant projections into reduction slot spaco")
    Seurat[["spaco_svg"]] <-
      Seurat::CreateDimReducObject(
        embeddings = SpaCoObject@projection[, 1:max(nSpacs, 2)],
        key = "Spac_svg_",
        assay = Seurat::DefaultAssay(Seurat)
      )

    return(Seurat)
  }

seurat_to_spatialpca <- function(seur) {
  counts <- seur@assays$SCT$counts
  cords <- as.data.frame(seur@images[[1]]@coordinates)
  cords <- cords[cords$tissue == 1, c("row", "col")]
  cords <- as.matrix(cords)
  SpatialPCAObj <-
    CreateSpatialPCAObject(counts, cords)
  SpatialPCAObj = SpatialPCA_buildKernel(
    SpatialPCAObj,
    kerneltype = "gaussian",
    bandwidthtype = "SJ",
    bandwidth.set.by.user = NULL
  )
  SpatialPCAObj = SpatialPCA_EstimateLoading(SpatialPCAObj, fast = FALSE,
                                             SpatialPCnum = 20)
  SpatialPCAObj = SpatialPCA_SpatialPCs(SpatialPCAObj, fast = FALSE)

  # SpatialPCAProjections <- SpatialPCAObj@normalized_expr %*%
  #   t(SpatialPCAObj@SpatialPCs)
  SpatialPCAProjections <- t(SpatialPCAObj@SpatialPCs)
  colnames(SpatialPCAProjections) <-
    paste0("spatialpca_", 1:ncol(SpatialPCAProjections))

  seur[["spatialpca"]] <-
    Seurat::CreateDimReducObject(
      embeddings = SpatialPCAProjections,
      key = "spatialpca_",
      assay = Seurat::DefaultAssay(seur)
    )
  spatialpca_kernel <- SpatialPCAObj@kernelmat
  list(seur = seur,
       kernel = spatialpca_kernel,
       svgs = rownames(SpatialPCAObj@normalized_expr))
}

process_data <-
  function(map_data,
           seur,
           vars_to_regress = c("^mt-", "^Hbb-"),
           variable_features_n = 3000,
           assay = "Spatial") {
    seur_1 <- seur_group_pipeline(seur,
                                  map_data,
                                  1,
                                  vars_to_regress,
                                  variable_features_n,
                                  assay)
    seur_2 <- seur_group_pipeline(seur,
                                  map_data,
                                  2,
                                  vars_to_regress,
                                  variable_features_n,
                                  assay)
    spatialpca_ls_1 <- seurat_to_spatialpca(seur_1)
    spatialpca_ls_2 <- seurat_to_spatialpca(seur_2)

    seur_1 <- spatialpca_ls_1[[1]]
    seur_2 <- spatialpca_ls_2[[1]]

    spaco_obj_1 <- seurat_to_spaco_spatialpca(spatialpca_ls_1, assay)
    spaco_obj_2 <- seurat_to_spaco_spatialpca(spatialpca_ls_2, assay)

    spaco_obj_1 <- RunSCA(spaco_obj_1)
    spaco_obj_2 <- RunSCA(spaco_obj_2)

    spaco_obj_1_svgs <- spaco_svg_projs(spaco_obj_1,
                                        spatialpca_ls_1$svgs)
    spaco_obj_2_svgs <- spaco_svg_projs(spaco_obj_2,
                                        spatialpca_ls_2$svgs)

    seur_1 <- spacs_to_seurat(spaco_obj_1, seur_1)
    seur_2 <- spacs_to_seurat(spaco_obj_2, seur_2)

    seur_1 <- spacs_to_seurat_svg(spaco_obj_1_svgs, seur_1)
    seur_2 <- spacs_to_seurat_svg(spaco_obj_2_svgs, seur_2)

    return(list(
      seur_1 = seur_1,
      seur_2 = seur_2
    ))
  }

liver_res <- process_data(liver_map, liver_seur, assay = "RNA")
anterior_res <- process_data(anterior_map, anterior_slice_seur)
posterior_res <- process_data(posterior_map, posterior_slice_seur)

anterior_res$seur_1 <- RunPCA(anterior_res$seur_1, verbose = FALSE)
anterior_res$seur_2 <- RunPCA(anterior_res$seur_2, verbose = FALSE)

save(
  liver_res,
  anterior_res,
  posterior_res,
  file = "/home/rstudio/SpaCo_paper_code/Clustering/Data/cluster_benchmark_objects.Rdata"
)
