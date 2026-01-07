# pull base image

FROM rocker/tidyverse:4.1.3
#FROM tschaffter/rstudio



# who maintains this image
LABEL maintainer Luis Montano "luis.montano@ccri.at"
LABEL version 4.1.3-v1

# how to build this image
#docker build -t cancerbits/dockr:ML2CellR -f ml2cellR.Dockerfile .

# change some permissions
RUN chmod -R a+rw ${R_HOME}/site-library # so that everyone can dynamically install more libraries within container
RUN chmod -R a+rw ${R_HOME}/library

# add custom options for rstudio sessions
# make sure sessions stay alive forever
RUN echo "session-timeout-minutes=0" >> /etc/rstudio/rsession.conf
# make sure authentication is not needed so often
RUN echo "auth-timeout-minutes=0" >> /etc/rstudio/rserver.conf
RUN echo "auth-stay-signed-in-days=365" >> /etc/rstudio/rserver.conf

RUN apt-get update && \
    apt-get install -y --no-install-recommends gnupg
    
#RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com --recv-keys 467B942D3A79BD29

##important line needed to upgrade from 4.1.3 to 4.3.2
#RUN apt-get update && \
#    apt-get install -y --no-install-recommends gnupg ca-certificates ubuntu-keyring && \
#    rm -rf /var/lib/apt/lists/*

##updating certificate...to fix several keyring bugs    
#RUN apt-key adv --keyserver pgp.mit.edu --recv-keys 467B942D3A79BD29

#install a million os level dependencies before bioconductor software

RUN apt-get update -qq \
  && apt-get -y --no-install-recommends install \
  htop \
  nano \
  libigraph-dev \
  libcairo2-dev \
  libxt-dev \
  libcurl4-openssl-dev \
  libcurl4 \
  libxml2-dev \
  libxt-dev \
  openssl \
  libssl-dev \
  wget \
  curl \
  bzip2 \
  libbz2-dev \
  libpng-dev \
  libhdf5-dev \
  pigz \
  libudunits2-dev \
  libgdal-dev \
  libgeos-dev \
  libboost-dev \
  libboost-iostreams-dev \
  libboost-log-dev \
  libboost-system-dev \
  libboost-test-dev \
  libz-dev \
  libarmadillo-dev \
  libglpk-dev \
  jags \
  libgsl-dev \
   libgsl0-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  cmake \
 libxt6 \
libglpk40 \
imagemagick \
  libmagick++-dev \
  libmagickcore-dev \
  default-jdk \
  liblzma-dev \
  python3 \
  python3-pip \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
##



# install conda (mambaforge):
RUN wget -O conda.sh https://github.com/conda-forge/miniforge/releases/download/25.3.0-3/Miniforge3-25.3.0-3-Linux-x86_64.sh && \
/bin/sh conda.sh -b -p /opt/conda 
RUN chmod -R a+rw /opt/conda
ENV PATH=/opt/conda/bin:$PATH


RUN R -e "install.packages(c('reticulate'))"

# create a conda environment for scVI:
RUN R -e "reticulate::conda_create(envname='/opt/conda/envs/scvi', python_version='3.11.5', pip=TRUE)" && \
  R -e "reticulate::conda_install(envname='/opt/conda/envs/scvi', packages=c('scvi_tools==0.20.3', 'scanpy==1.9.3'), pip=TRUE)"

# create a conda environment for arboreto/GRNboost2:
RUN R -e "reticulate::conda_create(envname='/opt/conda/envs/arboreto', python_version='3.8.17', pip=TRUE)" && \
  R -e "reticulate::conda_install(envname='/opt/conda/envs/arboreto', packages=c('arboreto==0.1.6'), pip=TRUE)"


# create a conda environment for deeptools:
#RUN R -e "reticulate::conda_create(envname='/opt/conda/envs/py37', python_version='3.7', pip=TRUE, channel=c('bioconda','conda-forge'))" && \
#  R -e "reticulate::conda_install(envname='/opt/conda/envs/py37', packages=c('deeptools==3.5.4'), pip=TRUE)"

# Install deepTools and ucsc-bigwigmerge in a Conda environment
RUN conda create -n py37 python=3.7 -y && \
    conda install -n py37 -c bioconda -c conda-forge deeptools=3.5.4 ucsc-bigwigmerge -y

RUN conda create -n nfcore python=3.12.1 -y && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict && \
    conda install -n nfcore nextflow=23.10.1 nf-core -y

RUN echo "\
channels:\n\
  - conda-forge\n\
  - bioconda\n\
  - defaults\n\
channel_priority: strict\n" > ~/.condarc

# create a conda environment for nfcore:
#RUN R -e "reticulate::conda_create(envname='/opt/conda/envs/nfcore', python_version='3.12.1', pip=TRUE, channel=c('bioconda', 'conda-forge'))" && \
#  R -e "reticulate::conda_install(envname='/opt/conda/envs/nfcore', packages=c('nexftlow==23.10.1', 'nf-core'))"


#RUN conda create -n py37 python=3.7 

#RUN echo "source activate env" > ~/.bashrc
#ENV PATH /opt/conda/envs/env/bin:$PATH

#environment location: /opt/miniconda/envs/py37

# make RUN commands use the new environment
#SHELL ["conda", "run", "-n", "py37", "/bin/bash", "-c"]


#RUN conda activate py37
#RUN conda install -c bioconda samtools
#RUN conda install -c bioconda deeptools
#RUN conda install -c bioconda ucsc-bigwigmerge
 

# The default CRAN mirror is set to a snapshot from 2022-04-21
# Bioconductor is at version 3.14

# install a couple of R packages that we commonly use
#RUN apt-get -y update && apt-get -y install \
#  libglpk40 \
#  libxt6 \
#  && apt-get clean \
#  && rm -rf /tmp/* /var/tmp/*

RUN R -e "install.packages(c('remotes', 'ggVenn', 'ggVennDiagram'))"
RUN R -e "install.packages('BiocManager', version='3.14')"
RUN R -e "BiocManager::install(c('markdown', 'sparseMatrixStats', 'edgeR', 'apeglm', 'DESeq2', 'fgsea', 'hypeR', 'patchwork','data.table', 'dplyr', 'tidyr','Rhtslib', 'biovisBase', 'Rsamtools','GenomicAlignments', 'BSgenome', 'VariantAnnotation','rtracklayer', 'GenomicFeatures', 'OrganismDbi', 'ensembldb', 'ggplot2', 'ggrastr','ggsignif','ggrepel', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'ggrepel','TFBSTools', 'ggwordcloud','GSVA', 'motifmatchr', 'BSgenome.Hsapiens.UCSC.hg38', 'signatureSearch','motifmatchr', 'qusage'))"
RUN R -e "install.packages(c('devtools',  'ggwordcloud', 'pROC', 'simpleCache', 'randomForest', 'ggVenDiagram' ))"
RUN R -e "devtools::install_github(repo = 'cancerbits/canceRbits', ref = '1f3dd76')"
RUN R -e "remotes::install_github(repo = 'jokergoo/complexHeatmap')"
#RUN R -e "install.packages('tinytex', repos='https://cloud.r-project.org')"
#RUN R -e "tinytex::install_tinytex()"


# more package installation below this line

#ENTRYPOINT ["conda", "run", "-n", "myenv", "py37", "src/server.py"]
