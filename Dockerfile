FROM rocker/rstudio:4.4

# Set environment variables 
# non interactive mode is for dealing with arguements that require response 
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies for RStudio Server
RUN apt-get update && apt-get install -y \
    psmisc \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libfftw3-dev \
    libgsl-dev \
    libclang-dev \
    libharfbuzz-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libfontconfig1-dev \
    libpng-dev \
    libfribidi-dev \
    libglpk40 \
    git \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install R packages

# Install devtools package in R
RUN R -e "install.packages('remotes')"
RUN R -e "install.packages('devtools')"

# Installing Necessary R packages
RUN R -e "library(remotes); remotes::install_version('ggplot2', version='3.5.2')"
#RUN R -e "library(remotes); remotes::install_version('Seurat', version='5.3.1')"
RUN R -e "remotes::install_github('satijalab/seurat@v5.3.1')"
RUN R -e "library(remotes); remotes::install_version('Matrix', version='1.7.3')"
RUN R -e "library(remotes); remotes::install_github('IMSBCompBio/SpaCo')"
RUN R -e "library(remotes); remotes::install_github('xzhoulab/SPARK')" 
RUN R -e "library(devtools); devtools::install_github('https://github.com/cz0316/HEARTSVG',force = T)" 
RUN R -e "library(remotes); remotes::install_version('tidyverse', version='2.0.0')"
RUN R -e "library(remotes); remotes::install_version('aricode', version='1.0.3')"
RUN R -e "library(remotes); remotes::install_version('clue', version='0.3.66')"
RUN R -e "library(remotes); remotes::install_version('mclust', version='6.1.1')"
RUN R -e "library(remotes); remotes::install_version('proxy', version='0.4.27')"
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install(c( \
    'S4Vectors', \
    'IRanges', \
    'XVector', \
    'S4Arrays', \
    'GenomeInfoDb', \
    'SparseArray', \
    'GenomicRanges', \
    'DelayedArray', \
    'DelayedMatrixStats', \
    'SummarizedExperiment', \
    'SingleCellExperiment', \
    'scuttle', \
    'splatter', \
    'glmGamPoi' \
    ))"

RUN R -e "library(remotes); remotes::install_github('shangll123/SpatialPCA')"

# 0.4.27
# Home Directory
WORKDIR /home/rstudio

# Make directory for figure generation scripts
# RUN mkdir -p /home/rstudio/figure_generation_scripts

# load data 
COPY Code_To_Kiarash /home/rstudio/SpaCo_paper_code
#COPY R_scripts/fig2_scripts /home/rstudio/figure_generation_scripts
#COPY R_scripts/brain_seuratobject_raw.RData /home/rstudio
#COPY R_scripts/fig2_scripts/Kia*/ home/rstudio/figure_generation_scripts/
#COPY R_scripts/SpaCo_R_kia /home/rstudio/SpaCo_R_kia

# Expose RStudio Server port
EXPOSE 8787

# run the server-daemonize in the foreground
# and disable authentication
RUN echo 'server-daemonize=0' >> /etc/rstudio/rserver.conf && \
    echo 'auth-none=0' >> /etc/rstudio/rserver.conf

# Set default user and password
ENV USERNAME=rstudio

# Keep the container running
CMD ["/init"]