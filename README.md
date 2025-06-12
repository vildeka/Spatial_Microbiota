# The cervicovaginal microbiome impacts spatially restricted host transcriptional signatures throughout the human ectocervical epithelium and submucosa

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

doi:[...](https://doi.org/)

This repository contains the code related to the article "The cervicovaginal microbiome impacts spatially restricted host transcriptional signatures throughout the human ectocervical epithelium and submucosa"

This project used:
-   21 samples of Spatial transcriptomics data (10x Visium)

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
1. [Figure_1.Rmd](https://vildeka.github.io/Spatial_Microbiota/Figure1)
2. [Figure_2&3.Rmd](https://vildeka.github.io/Spatial_Microbiota/Figure2&3)
3. [Figure_4.Rmd](https://vildeka.github.io/Spatial_Microbiota/Figure4)
4. [Figure_5.Rmd](https://vildeka.github.io/Spatial_Microbiota/Figure5)
5. [Figure_6.Rmd](https://vildeka.github.io/Spatial_Microbiota/Figure6)

## Dependencies

Project is created with:

-   R version: 4.3.3
-   RStudio version: 2025.05.0 (use version seperate from conda env.)
-   renv version: 0.15.2
-   Seurat version: 4.4.0

## Data Availability Statement

**Spattial transcriptomics count data** files and RDS object can be accessed in the Gene Expression Omnibus public repository, SuperSeries ID GSE217237. The raw transcriptomic sequencing data cannot be held in a public repository due to the sensitive nature of such personal data. Request for data access can be made to the Karolinska Institutet Research Data Office (contact via rdo\@ki.se), and access will be granted if the request meets the requirements of the data policy.

The 21 samples have been published in three separate GEO records as follows:
GSE (12 new sample): P020,P045,P050,P057,P001,P014,P018,P087,P021,P024,P081,P117
GSE (4 samples): P118, P105, P080, P031
GSE290350 (5 samples): P004, P008, P026, P044, P067

## Repo description

-   **src**\
    contains all the analysis scripts
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
