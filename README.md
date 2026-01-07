# Supplementary code repository for: Directing stem cell differentiation by chromatin state approximation

Authors: Luis F. Montano-Gutierrez<sup>1,\*</sup>, Sophie Müller<sup>1,2,\*</sup>, Ana P. Kutschat<sup>1,2</sup>, Igor Adameyko<sup>3,4</sup>, Davide Seruggia<sup>1,2,#</sup>, Florian Halbritter<sup>1,#</sup>

1 St. Anna Children’s Cancer Research Institute (CCRI), Vienna, Austria

2 CeMM Research Center for Molecular Medicine of the Austrian Academy of Sciences, Vienna, Austria

3 Department of Physiology and Pharmacology, Karolinska Institute, Solna, Sweden

4 Department of Neuroimmunology, Center for Brain Research, Medical University of Vienna, Vienna, Austria

<sup>\*,#</sup>, These authors contributed equally to this work.


## Abstract:

A prime goal of regenerative medicine is to replace dysfunctional cells in the body. To design protocols for producing target cells in the laboratory, one may need to consider exponentially large combinations of culture components. Here, we tested the potential of iteratively approximating the target phenotype by quantifying the distance between chromatin profiles (ATAC-seq) of differentiating cells in vitro and their in-vivo counterparts. We tested this approach on the well-studied generation of erythroblasts from haematopoietic stem cells, evaluating a fixed number of components over two sequential differentiation rounds (8x8 protocols). We found that the most erythroblast-like cells upon the first round yielded the most erythroblast-like cells at the second round, suggesting that greedy selection by chromatin approximation can be a viable optimisation strategy. Furthermore, by targeting transcriptional regulators linked to chromatin regions that were incompletely reprogrammed even after two rounds of differentiation, we could make a data-driven selection of additions to the protocol that further improved erythropoiesis. In future, our methodology can help craft notoriously difficult cells in vitro, such as B cells.

## Repository structure:

* `ml2cellR.Dockerfile` defines the environment used to carry out R analyses 
* `config.yaml` is used to set paths 
* `R/` holds R function definitions and misc utility scripts
* `Rmd/` holds R markdown documents for the individual steps of the project
* `nextflow/` has custom nextflow pipelines
* `bash/` holds shell scripts to build and run the docker image, and to parse the config file
* `metadata/` holds custom geneset definitions

## Links:

* Gene Expression Omnibus (GEO): [GSE291386](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291386)
* Preprint (bioRxiv): [doi:10.1101/2025.04.24.650451][https://doi.org/10.1101/2025.04.24.650451]

# ML2Cell

## Table of Contents

1.  [Using ML2Cell](#using-ml2cell)
    -   [1.1 Purpose](#11-purpose)\
    -   [1.2 Overview of Steps](#12-overview-of-the-steps)\
    -   [1.3 Required Datasets](#13-collecting-datasets-for-ml2cell)\
    -   [1.4 Pre-analysing ATAC-seq
        Data](#14-pre-analysing-your-datasets-with-an-atac-seq-pipeline)
2.  [Running ML2Cell](#running-ml2cell--step-by-step-instructions)
    -   [2.0.1 Manual
        Clone](#201-cloning-the-git-repository-manually-no-command-line-experience-required)\
    -   [2.0.2 Cloning via terminal](#202-cloning-the-git-repository-via-terminal-command-line-experience-required)\
    -   [2.1 Setup](#21-setup)
        -   [2.1.1
            dataset_info.csv](#211-findcreate-folders-in-the-computer-and-fill-in-the-dataset_infocsv)\
        -   [2.1.2
            checking and preparing the config yaml ](#212-check-all-the-input-options-in-the-config-file)\
        -   [2.1.3 Container
            Setup](#213-setting-up-the-container-environment-to-run-the-code)\
        -   [2.1.4
            Samplesheets](#214-preparing-a-samplesheet-for-each-dataset)
3.  [Adjustments After First
    Run](#3-potential-adjustments-after-the-first-run)

------------------------------------------------------------------------

# Using ML2Cell

## 1. Overview

### 1.1 Purpose

1.  ML2Cell takes your ATAC-seq samples plus a cell-type reference
    (including a *target* cell type). It outputs how well each sample
    matches the target cell type's defining chromatin features.

2.  ML2Cell identifies chromatin discrepancies and suggests treatments
    that may help shift your samples toward the target cell type.

------------------------------------------------------------------------

### 1.2 Overview of the Steps

1.  Run nf-core/ATAC-seq on all datasets *(future built-in step)*\
2.  Map sample reads to reference consensus peaks\
3.  Identify cell-type-informative peaks using ROC analysis and
    differential accessibility\
4.  Compute each sample's distance to the target in a
    cell-type-informative PCA\
5.  Perform signature analysis to identify residual cell-type
    contributions\
6.  Identify chromatin peaks discordant with the reference\
7.  Detect enriched transcription factor binding sites in discordant
    regions\
8.  Suggest treatments targeting priority transcription factors

------------------------------------------------------------------------

### 1.3 Collecting Datasets for ML2Cell

ML2Cell accepts three dataset categories:

-   **Reference** --- ATAC-seq of your reference cell types\
-   **Test** --- your samples\
-   **External** --- optional published datasets to analyze separately

------------------------------------------------------------------------

### 1.4 Pre-analysing Your Datasets with an ATAC-seq Pipeline

ML2Cell expects your datasets to be preprocessed.\
We **strongly recommend** using the **nf-core/atacseq** Nextflow
pipeline.

------------------------------------------------------------------------

# Running ML2Cell --- Step-by-Step Instructions

## 2.0.1 Cloning the Git Repository Manually (No Command Line Experience Required)

Download the repository:\
**https://github.com/cancerbits/ml2cell_code**\
Click **"Download this repository"**.

------------------------------------------------------------------------

## 2.0.2 Cloning the Git Repository via Terminal (Command Line Experience Required)

Install Git: https://git-scm.com/\
Then:

``` bash
git clone https://github.com/cancerbits/ml2cell_code
```

------------------------------------------------------------------------

# 2.1 Setup

## 2.1.1 Find/Create Folders in the Computer and Fill in the `dataset_info.csv`

For each dataset, provide:

-   Path to folder containing ATAC-seq processed data\
-   Path to metadata (samplesheets)\
-   Path to nf-core/atacseq outputs\
-   Path to a clean ML2Cell output folder
    -   *If not provided, outputs go into the current folder*

Use `dataset_info_example.csv` as a guide.

------------------------------------------------------------------------

## 2.1.2 Check ALL Input Options in the Config File

**IMPORTANT:**\
Complete the entire configuration file, especially **target cell type
parameters**.

Recommendation:\
Run with default settings for the first pass.\
Later, modify parameters such as QC thresholds or Random Forest
classification settings.

------------------------------------------------------------------------

## 2.1.3 Setting Up the Container Environment to Run the Code

### Docker, Podman, and Containerisation Overview

ML2Cell runs inside a Docker container.\
More information: https://www.docker.com\
Podman is an open-source alternative.

------------------------------------------------------------------------

### Using the ML2Cell Container

#### Strategy 1 (Recommended): Preconstructed Docker Image

1.  Install Docker Desktop or Docker Engine\
2.  Pull the public ML2Cell image:

``` bash
docker pull cancerbits:ML2cellR
```

------------------------------------------------------------------------

#### Strategy 2 (Advanced): Local Docker + Repository Setup

*(Requires completed `config.yaml`)*

Run:

``` bash
bash run_docker_rstudio.sh <PORT> <PASSWORD>
```

Using Podman:

``` bash
bash run_podman_rstudio.sh <PORT> <PASSWORD>
```

------------------------------------------------------------------------

## 2.1.4 Preparing a Samplesheet for Each Dataset

Samplesheets must be comma-separated `.csv` files with the columns:

    fastq_1, fastq_2, Experiment, celltype, Condition

For test/validation samples with treatments, add:

    treatment_id.<round>

Example:

    hydrocortisone.1

Treatment value definitions:

-   `NA` → not part of treatment group (e.g. reference)
-   `0` → control group, no treatment
-   `1` → treatment applied at 1×
-   any number → treatment applied at that fold (e.g., `3` = 3× dose)

------------------------------------------------------------------------

# 3. Potential Adjustments After the First Run

Consider modifying the following in `config.yaml`:

### `randomforest_gini_threshold`

-   Controls number of PCA components used in classification\
-   Should stay **\> 1** but can be permissive

### `randomforest_ntrees`

-   Higher values may refine classification\
-   increases computation time

------------------------------------------------------------------------

ML2Cell depends on **nf-core/ATAC-seq** processed data and correct
configuration via `config.yaml`.
