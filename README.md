SPaCo analysis and visualisation of spatial sequencing data
================

# SPaCo guided tutorial

## Installation

This beta of SPACO is currently not yet available on Bioconductor or
CRAN. Therefore it must be installed from github via devtools.

``` r
library(devtools)
devtools::install_github("IMSBCompBio/SpaCo")
```

## Setup the SPaCo object and normalize the data.

### Setup a SPaCo object from 10x genomics Visium data.

SPACO uses the importer function of the Seurat library to import 10X
Visium spatial data. The imported data is preprocessed using Variance
Stabilizing Transformations by the sctransform package (URAL) The number
of variable features to keep can be defined and any variables that
should be regressed out. In out tutorial we correct for mitochondrial
and hemoglobin genes detected in the spots

``` r
#Specify data paths
library(SPACO)

data_dir="~/path/to/directory"
slice = "~/path/to/directory/with/the/H&E/stain"
filename ="name of the feature matrix" 
spatial_file = "name of the tissue positions list"

#Initialize the object from raw (non-normalized data and use the Seurat library for pre-processing)
SpaCoObject <- read_10x_for_spaco(data_dir = data_dir,slice = slice,filename = filename,
                                  variable_features_n = 3000,
                                  spatial_file = spatial_file,
                                  vars_to_regress = c("^mt-","^Hbb-"))
```

### Setup a SPaCo object from existing Seurat object

``` r
library(Seurat)
library(SPACO)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

#install data from SeuratData
#InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")
brain <- PercentageFeatureSet(brain, pattern = "^mt-" ,col.name = "percent.mt")
brain <- PercentageFeatureSet(brain, pattern = "^Hbb-" ,col.name = "percent.hbb")
brain <- SCTransform(brain, assay = "Spatial", variable.features.n = 3000, verbose = F)

SpaCoObject <- seurat_to_spaco(Seurat = brain, assay = "SCT", n_image= 1, slot = "scale.data")
```

### Run the Spatial component analysis as dimensionality reduction and computer the number of informative spatial components.

In this step we actually do the spatial component analysis. Further, we
determine the number of relevant spatial components via sampling. These
can be found in the @nSpacs slot in the object.

``` r
SpaCoObject <- RunSCA(SpaCoObject, PC_criterion = "percent",
                      PC_value = .8, compute_nSpacs = T, nSpacQuantile = 0.05, nSim = 110 )
```

    ## computing number of releveant spacs

``` r
SpaCoObject@nSpacs
```

    ## [1] 47

### Visualisation and denoising

we can plot the computed meta gene projections directly from the
SPaCoObject or transfer the projections right back into an existing
Seurat Object.

``` r
spacplot <- Spaco_plot(SpaCoObject, spac = 1:4, ncol = NULL, combine = T)
```

    ## Loading required package: scales

``` r
spacplot
```

![](README_files/figure-gfm/plot%20the%20first%204%20meta%20genes-1.png)<!-- -->

``` r
SpaCoObject <-  denoise_profiles(SpaCoObject)
denoised <- denoised_projection_plot(SpaCoObject,features = "Ybx1")+ggtitle("denoised")
original <- feature_plot(SpaCoObject, features = "Ybx1", ncol = NULL, combine = TRUE)+ggtitle("original")
original+denoised
```

![](README_files/figure-gfm/plot%20the%20first%204%20meta%20genes-2.png)<!-- -->

### Add dimension reduction information to an existing Seurat object.

First we remove all spots from the Seurat object which have no direct
neighbors as they violate SPaCo assumptions. The computed SPaCo
projections are stored in the object slot “object\[\[”spaco”\]\]”.

``` r
brain <- SPACO::subset_non_neighbour_cells(SpaCoObject,brain)

brain <- SPACO::spacs_to_seurat(SpaCoObject,brain)
```

    ## copying significan projections into reduction slot spaco

### Use Seurat for visualisation

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

brain <- FindNeighbors(brain,reduction = "spaco" , dims = 1:SpaCoObject@nSpacs, verbose = F)
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

### Compute spatial variable genes (SVG’s).

``` r
DE_genes<- SVGTest(SpaCoObject)
```

    ## Loading required package: mgcv

    ## Loading required package: nlme

    ## 
    ## Attaching package: 'nlme'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse

    ## This is mgcv 1.9-0. For overview type 'help("mgcv-package")'.

    ## Loading required package: CompQuadForm

``` r
head(DE_genes)
```

    ##                    score     pVal p.adjust
    ## Rgs20         0.21917877 0.000000        0
    ## Oprk1         0.24688536 0.000000        0
    ## Npbwr1        0.04620468 0.000000        0
    ## St18          0.07229480 0.000000        0
    ## 3110035E14Rik 0.74777639 0.000000        0
    ## 1700034P13Rik 0.01280162 0.781657        1

``` r
DE_genes["Ppp1r1b",]
```

    ##             score pVal p.adjust
    ## Ppp1r1b 0.9162679    0        0

### Test for Spatially variable genes

### GO-Term and KEGG pathway enrichment of spatially variable genes
