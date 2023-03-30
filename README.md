SPaCo analysis and visualisation of spatial sequencing data
================

# SPaCo guided tutorial

## Setup the SPaCo object and normalize the data.

### Setup a SPaCo object from 10x genomics Visium data.

``` r
#Specify data paths
library(SPACO)

data_dir="~/Liver_Organoid_project/scRNA_spatial_tissue/VX03_A006200136/cellranger/141895/outs"
slice = "~/Liver_Organoid_project/scRNA_spatial_tissue/VX03_A006200136/cellranger/141895/outs/spatial"
filename ="filtered_feature_bc_matrix.h5" 
spatial_file = "tissue_positions_list.csv"

#Initialize the object from raw (non-normalized data and use the Seurat library for pre-processing)
SpaCoObject <- read_10x_for_spaco(data_dir = data_dir,slice = slice,filename = filename,
                                  variable_features_n = 3000,
                                  spatial_file = spatial_file,
                                  vars_to_regress = c("^mt-","^Hbb-"))
```

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   4%  |                                                                              |=====                                                                 |   7%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  18%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  25%  |                                                                              |====================                                                  |  29%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  75%  |                                                                              |=======================================================               |  79%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  96%  |                                                                              |======================================================================| 100%
    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   4%  |                                                                              |=====                                                                 |   7%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  18%  |                                                                              |===============                                                       |  21%  |                                                                              |==================                                                    |  25%  |                                                                              |====================                                                  |  29%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  39%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  46%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================                                |  54%  |                                                                              |========================================                              |  57%  |                                                                              |==========================================                            |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  71%  |                                                                              |====================================================                  |  75%  |                                                                              |=======================================================               |  79%  |                                                                              |==========================================================            |  82%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |=================================================================     |  93%  |                                                                              |====================================================================  |  96%  |                                                                              |======================================================================| 100%

### Setup a SPaCo object from existing Seurat object

``` r
library(Seurat)
library(SeuratData)
```

    ## ── Installed datasets ───────────────────────────────────── SeuratData v0.2.2 ──

    ## ✔ stxBrain 0.1.1

    ## ────────────────────────────────────── Key ─────────────────────────────────────

    ## ✔ Dataset loaded successfully
    ## ❯ Dataset built with a newer version of Seurat than installed
    ## ❓ Unknown version of Seurat installed

``` r
library(ggplot2)
library(patchwork)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
#install data from SeuratData
InstallData("stxBrain")
```

    ## Warning: The following packages are already installed and will not be
    ## reinstalled: stxBrain

``` r
brain <- LoadData("stxBrain", type = "anterior1")
brain <- PercentageFeatureSet(brain, pattern = "^mt-" ,col.name = "percent.mt")
brain <- PercentageFeatureSet(brain, pattern = "^Hbb-" ,col.name = "percent.hbb")
brain <- SCTransform(brain, assay = "Spatial", variable.features.n = 3000, 
                     vars.to.regress = c("percent.mt","percent.hbb"))
```

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 17668 by 2696

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 2000 genes, 2696 cells

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%

    ## Found 90 outliers - those will be ignored in fitting/regularization step

    ## Second step: Get residuals using fitted parameters for 17668 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |====                                                                  |   6%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  17%  |                                                                              |==============                                                        |  19%  |                                                                              |================                                                      |  22%  |                                                                              |==================                                                    |  25%  |                                                                              |===================                                                   |  28%  |                                                                              |=====================                                                 |  31%  |                                                                              |=======================                                               |  33%  |                                                                              |=========================                                             |  36%  |                                                                              |===========================                                           |  39%  |                                                                              |=============================                                         |  42%  |                                                                              |===============================                                       |  44%  |                                                                              |=================================                                     |  47%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  53%  |                                                                              |=======================================                               |  56%  |                                                                              |=========================================                             |  58%  |                                                                              |===========================================                           |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |===============================================                       |  67%  |                                                                              |=================================================                     |  69%  |                                                                              |===================================================                   |  72%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================                |  78%  |                                                                              |========================================================              |  81%  |                                                                              |==========================================================            |  83%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |================================================================      |  92%  |                                                                              |==================================================================    |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%

    ## Computing corrected count matrix for 17668 genes

    ##   |                                                                              |                                                                      |   0%  |                                                                              |==                                                                    |   3%  |                                                                              |====                                                                  |   6%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  11%  |                                                                              |==========                                                            |  14%  |                                                                              |============                                                          |  17%  |                                                                              |==============                                                        |  19%  |                                                                              |================                                                      |  22%  |                                                                              |==================                                                    |  25%  |                                                                              |===================                                                   |  28%  |                                                                              |=====================                                                 |  31%  |                                                                              |=======================                                               |  33%  |                                                                              |=========================                                             |  36%  |                                                                              |===========================                                           |  39%  |                                                                              |=============================                                         |  42%  |                                                                              |===============================                                       |  44%  |                                                                              |=================================                                     |  47%  |                                                                              |===================================                                   |  50%  |                                                                              |=====================================                                 |  53%  |                                                                              |=======================================                               |  56%  |                                                                              |=========================================                             |  58%  |                                                                              |===========================================                           |  61%  |                                                                              |=============================================                         |  64%  |                                                                              |===============================================                       |  67%  |                                                                              |=================================================                     |  69%  |                                                                              |===================================================                   |  72%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================                |  78%  |                                                                              |========================================================              |  81%  |                                                                              |==========================================================            |  83%  |                                                                              |============================================================          |  86%  |                                                                              |==============================================================        |  89%  |                                                                              |================================================================      |  92%  |                                                                              |==================================================================    |  94%  |                                                                              |====================================================================  |  97%  |                                                                              |======================================================================| 100%

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 51.03828 secs

    ## Determine variable features

    ## Place corrected count matrix in counts slot

    ## Regressing out percent.mt, percent.hbb

    ## Centering data matrix

    ## Set default assay to SCT

``` r
SpaCoObject <- seurat_to_spaco(Seurat = brain, assay = "SCT", n_image= 1, slot = "scale.data")
```

### Run the Spatial component analysis as dimensionality reduction and computer the number of informative spatial components.

Here we use the new implementation to compute the numer of relevant
spac’s. If you want to use the standard implementation just use RunSCA.

``` r
SpaCoObject<- RunSCA2(SpaCoObject, PC_criterion = "percent",
                      PC_value = .8, compute_nSpacs = T,
                      compute_projections = TRUE, nSpacQuantile = 0.05)
```

    ## computing number of releveant spacs

    ## computing projections

``` r
SpaCoObject@nSpacs
```

    ## [1] 52

we can plot the computed meta gene projections directly from the
SPaCoObject or transfer the projections right back into an existing
Seurat Object.

``` r
spacplot <- Spaco_plot(SpaCoObject, spac = 1:4, ncol = NULL, combine = T)
spacplot
```

![](README_files/figure-gfm/plot%20the%20first%204%20meta%20genes-1.png)<!-- -->

The number of relevant spatial components (spacs) is stored in the
@nSpacs slot in the SPaCoObject

``` r
SpaCoObject@nSpacs
```

    ## [1] 52

### Compute spatial variable genes (SVG’s).

``` r
DE_genes<- FindSVG(SpaCoObject, nSim = 1000, nSpacs = SpaCoObject@nSpacs)
```

    ## computing emprirical p-values this may take a while.

``` r
head(DE_genes)
```

    ##            Gene       score  pVal pAdjust
    ## 1         Rgs20 0.016793129 0.001   0.834
    ## 2         Oprk1 0.016852438 0.001   0.834
    ## 3        Npbwr1 0.019321797 0.015   0.998
    ## 4          St18 0.018823283 0.001   0.834
    ## 5 3110035E14Rik 0.007950325 0.001   0.834
    ## 6 1700034P13Rik 0.020306009 0.568   0.998

###Add dimension reduction information to an existing Seurat object.
First we remove all spots from the Seurat object which have no direct
neighbors as they violate SPaCo assumptions. The computed SPaCo
projections are stored in the object slot “object\[\[”spaco”\]\]”.

``` r
brain <- SPACO::subset_non_neighbour_cells(SpaCoObject,brain)

brain <- SPACO::spacs_to_seurat(SpaCoObject,brain)
```

    ## copying projections into reduction slot spaco

To compare the spatial component analysis to the analog in the Seurat
pipeline we do all the downstream processing right in the Seurat object.

``` r
cc <- SpatialFeaturePlot(brain,features = c("Spac_1","Spac_2","Spac_3","Spac_4"),combine = T)
brain_2 <- RunPCA(brain, verbose = F)
dd <- SpatialFeaturePlot(brain_2,features = c("PC_1","PC_2","PC_3","PC_4"),combine = T)

brain <- RunUMAP(brain,reduction = "spaco",dims = 1:SpaCoObject@nSpacs,n.neighbors = 45, verbose = F)
```

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

``` r
brain_2 <- RunUMAP(brain_2,reduction = "pca",dims = 1:30, verbose = F)

brain <- FindNeighbors(brain,reduction = "spaco" , dims = 1:5, verbose = F)
brain_2 <- FindNeighbors(brain_2,reduction = "pca" , dims = 1:30, verbose = F)

brain <- FindClusters(brain,resolution = 0.24, verbose = F)
brain_2 <- FindClusters(brain_2, verbose = F)

aa <- DimPlot(brain,group.by="seurat_clusters")+ggtitle("SpaCo")
bb <- DimPlot(brain_2,group.by="seurat_clusters")+ggtitle("Pca")
a <- DimPlot(brain,reduction = "spaco")+ggtitle("SpaCo")
b <- DimPlot(brain_2,reduction = "pca")+ggtitle("Pca")

rr <- SpatialDimPlot(brain)+ggtitle("SpaCo")
qq <- SpatialDimPlot(brain_2)+ggtitle("Pca")
```

``` r
patchwork::wrap_plots(a,b,aa,bb,rr,qq,ncol=2)
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
a+b
```

![](README_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
aa+bb
```

![](README_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
rr+qq
```

![](README_files/figure-gfm/unnamed-chunk-5-4.png)<!-- -->
