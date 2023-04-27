SPaCo analysis and visualisation of spatial sequencing data
================

# SPaCo guided tutorial

## Setup the SPaCo object and normalize the data.

### Setup a SPaCo object from 10x genomics Visium data.

### Setup a SPaCo object from existing Seurat object

``` r
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

#install data from SeuratData
InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")
brain <- PercentageFeatureSet(brain, pattern = "^mt-" ,col.name = "percent.mt")
brain <- PercentageFeatureSet(brain, pattern = "^Hbb-" ,col.name = "percent.hbb")
brain <- SCTransform(brain, assay = "Spatial", variable.features.n = 3000, verbose = F)

SpaCoObject <- seurat_to_spaco(Seurat = brain, assay = "SCT", n_image= 1, slot = "scale.data")
```

### Run the Spatial component analysis as dimensionality reduction and computer the number of informative spatial components.

Here we use the new implementation to compute the numer of relevant
spac’s that you can find in the @nSpacs slot. If you want to use the
standard implementation just use RunSCA and include the nSim parameter.
It’s recommended to use at least 1000 iterations. The RunSCA2 function
uses a statistical test to compute the significant spac’s but might be
less precise than the sampling approach but will have a much shorter
runtime.

``` r
SpaCoObject<- RunSCA2(SpaCoObject, PC_criterion = "percent",
                      PC_value = .8, compute_nSpacs = T,
                      compute_projections = TRUE, nSpacQuantile = 0.05)
```

    ## computing number of releveant spacs

    ## Using quantile level of 0.05

    ## computing projections

``` r
SpaCoObject@nSpacs
```

    ## [1] 52

``` r
SpaCoObject_sample<- RunSCA(SpaCoObject, PC_criterion = "percent",
                      PC_value = .8, compute_nSpacs = T,
                      compute_projections = TRUE, nSpacQuantile = 0.05, nSim = 1000)
```

    ## computing number of releveant spacs

    ## computing projections this may take a while

``` r
SpaCoObject_sample@nSpacs
```

    ## [1] 35

we can plot the computed meta gene projections directly from the
SPaCoObject or transfer the projections right back into an existing
Seurat Object.

``` r
spacplot <- Spaco_plot(SpaCoObject, spac = 1:4, ncol = NULL, combine = T)
spacplot
```

![](README_files/figure-gfm/plot%20the%20first%204%20meta%20genes-1.png)<!-- -->

``` r
SpaCoObject <-  smooth_profiles(SpaCoObject)
smoothed_projection_plot(SpaCoObject,features = "Ppp1r1b")
```

![](README_files/figure-gfm/plot%20the%20first%204%20meta%20genes-2.png)<!-- -->

``` r
feature_plot(SpaCoObject, features = "Ppp1r1b", ncol = NULL, combine = TRUE)
```

![](README_files/figure-gfm/plot%20the%20first%204%20meta%20genes-3.png)<!-- -->

###Add dimension reduction information to an existing Seurat object.
First we remove all spots from the Seurat object which have no direct
neighbors as they violate SPaCo assumptions. The computed SPaCo
projections are stored in the object slot “object\[\[”spaco”\]\]”.

``` r
brain <- SPACO::subset_non_neighbour_cells(SpaCoObject,brain)

brain <- SPACO::spacs_to_seurat(SpaCoObject,brain)
```

    ## copying significan projections into reduction slot spaco

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

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
a+b
```

![](README_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
aa+bb
```

![](README_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

``` r
rr+qq
```

![](README_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->

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
