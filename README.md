# The cervicovaginal microbiome impacts spatially restricted host transcriptional signatures across the human ectocervical epithelium and submucosa

Vilde Kaldhusdal, Mathias Franzen Boger, Adam D. Burgener, Julie Lajoie, Kenneth Omollo, Joshua Kimani, Annelie Tjernlund, Keith Fowke, Douglas S. Kwon, Gabriella Edfeldt and Kristina Broliden

![](./resources/Graphical%20abstract.png)

## Table of contents

-   [General info](#general-info)
-   [Workflow](#Workflow)
-   [Dependencies](#dependencies)
-   [Data Availability Statment](#data-availability-statment)
-   [Repo description](#repo-description)
-   [Setup](#setup)

## General info

doi:[10.1038/s41598-024-83775-9](https://doi.org/10.1038/s41598-024-83775-9)

This repository contains the code related to the article "The cervicovaginal microbiome impacts spatially restricted host transcriptional signatures across the human ectocervical epithelium and submucosa"

This project used multiple datasets:
-   Spatial transcriptomics data (10x Visium)
-   Transcriptomics data (bulk mRNA-SEQ)

## Workflow

### Analysis scripts

0.  [00_load_st_data](https://vildeka.github.io/Spatial_Microbiota/00_load_st_data)
1.  [01_QC_st_data](https://vildeka.github.io/Spatial_Microbiota/01_QC_st_data)
2.  [02_integrate_st_data](https://vildeka.github.io/Spatial_Microbiota/02_integrate_st_data)
3.  [03_clustering_st_data](https://vildeka.github.io/Spatial_Microbiota/03_clustering_st_data)
4.  [04_deconvolute_st_data](https://vildeka.github.io/Spatial_Microbiota/04_deconvolute_st_data)
5.  [05_DGE_clusters_st_data](https://vildeka.github.io/Spatial_Microbiota/05_DGE_clusters_st_data)
6.  [06_DGE_condition_st_data](https://vildeka.github.io/Spatial_Microbiota/06_DGE_condition_st_data)
7.  [07_hdWGCNA_analysis](https://vildeka.github.io/Spatial_Microbiota/07_hdWGCNA_analysis)
8.  [08_spatial_distance](https://vildeka.github.io/Spatial_Microbiota/08_spatial_distance)

### Manuscript figures
1. [Figure_1.Rmd](https://vildeka.github.io/Spatial_Microbiota/Figure_1)
2. [Figure_2.Rmd](https://vildeka.github.io/Spatial_Microbiota/Figure_2)
3. [Figure_3.Rmd](https://vildeka.github.io/Spatial_Microbiota/Figure_3)
4. [Figure_4.Rmd](https://vildeka.github.io/Spatial_Microbiota/Figure_4)
5. [Figure_5.Rmd](https://vildeka.github.io/Spatial_Microbiota/Figure_5)
6. [Figure_6.Rmd](https://vildeka.github.io/Spatial_Microbiota/Figure_6)

## Dependencies

Project is created with:

-   R version: 4.1.2
-   RStudio version: 1.4.1106 (use version seperate from conda env.)
-   renv version: 0.15.2
-   Seurat version:

## Data Availability Statement

**Spattial transcriptomics count data** files and RDS object can be accessed in the Gene Expression Omnibus public repository, SuperSeries ID GSE217237. The raw transcriptomic sequencing data cannot be held in a public repository due to the sensitive nature of such personal data. Request for data access can be made to the Karolinska Institutet Research Data Office (contact via rdo\@ki.se), and access will be granted if the request meets the requirements of the data policy.

**Bulk transcriptomics count data** files can be accessed in the Gene Expression Omnibus public repository, SuperSeries ID GSE217237. The raw transcriptomic sequencing data cannot be held in a public repository due to the sensitive nature of such personal data. Request for data access can be made to the Karolinska Institutet Research Data Office (contact via rdo\@ki.se), and access will be granted if the request meets the requirements of the data policy.

## Repo description

-   **src**\
    contains all the analysis script including all preprocessing steps
-   **manuscript**\
    reproducible code for figures included in the manuscript
-   **md_files**\
    rendered versions of the Rmds in src and manuscript folder respectively
-   **bin**\
    R scripts with various helper functions, mostly visualization functions for the ST data
-   **data**\
    Empty folder to store data downloaded from GEO

<!-- -->

```         
project
│   README.md
│   renv.loc    
└───src
│   │   00_Preprocessing.Rmd
│   │   02_Analysis.Rmd
│   │   ...
|   └───md files (rendered versions of the Rmds)
│   |       │   00_Preprocessing.md
│   |       │   02_Analysis.md
│   |       │   ...
|   |
│   └───manuscript
│       │   Figure01.Rmd
│       │   Figure02.Rmd
│       │   ...
|       └───md files
│           │   Figure01.md
│           │   Figure02.md
│           │   ...
└───bin
│   │   file01.txt
│   │   file02.txt
└───data
│   │   file01.txt
│   │   file02.txt
│
```

## Setup

Recommended setup to run this project:

Conda + renv

1.  Clone the repo
2.  If not already installed download mini conda/conda
3.  In the terminal navigate to the project directory
4.  create a new enviroment:<br/> `conda env create -n Spatial_Microbiota -f environment.yml`
5.  Activate the enviroment:<br/> `conda activate Spatial_Microbiota`
6.  Open Rstudio:<br/> `rstudio& Spatial DMPA`
7.  Install all packages specified by the lockfile:<br/> `renv::restore()`
