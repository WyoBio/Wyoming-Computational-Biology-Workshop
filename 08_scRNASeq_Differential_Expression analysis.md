# Differential expression analysis

In this section, we'll utilize the previously generated Seurat object that has undergone various preprocessing steps, clustering, and cell typing. We'll use this object for gene expression and differential expression analyses. Two key differences in the differential expression analyses conducted in this module are:

- **Identification of Epithelial Cells:** We'll start by identifying epithelial cells in our data using Epcam expression as a marker, considering that tumor cells should exhibit epithelial characteristics. Subsequently, we'll perform a differential expression analysis within the Epcam-positive population(s).

- **Comparison of T Cell Phenotypes:** Another analysis will compare T cell populations in mice treated with ICB therapy against those with T cells depleted, which underwent ICB therapy. This comparison aims to identify differences in T-cell phenotypes related to the treatment.

### Setup

If the saved Seurat object from the previous step isn't already loaded in your current R session, read it into the session.

```R
library(Seurat)
library(dplyr)
library(EnhancedVolcano)
library(presto)

datadir =
setwd(datadir)
merged <- readRDS('preprocessed_object.rds')
```

### Analysis of gene expression in epithelial cells.


