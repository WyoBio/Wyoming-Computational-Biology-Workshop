# Quality assessment

We are going to begin our single-cell analysis by loading in the output from CellRanger. We will load in our different samples, create a Seurat object with them, and take a look at the quality of the cells.

The general steps to preprocessing your single-cell data with Seurat:

1. **Create a Seurat object**
2. **Filter low-quality cells**
3. **Merge samples**
4. **Normalize counts**
5. **Find variable features**
6. **Scale data**
7. **Determine PCs for Clustering**
8. **Clustering** - FindNeighbors, FindClusters, RunUMAP

### Setup
Begin by initiating a fresh R script and loading the essential packages. Among these, Seurat takes the spotlight. Renowned for its distinctive data structure and robust tools, Seurat facilitates quality control, analysis, and in-depth exploration of single-cell RNA sequencing data. Widely embraced, it stands as a staple in the realm of single-cell RNA analysis, setting a benchmark for its comprehensive functionality and reliability.

```R
library("Seurat") 
library("ggplot2") # creating plots
library("cowplot") # add-on to ggplot, we use the plot_grid function to put multiple violin plots in one image file
library("dplyr")   # a set of functions for data frame manipulation -- a core package of Tidyverse (THE R data manipulation package)
library("Matrix")  # a set of functions for operating on matrices
library("viridis") # color maps for graphs that are more readable than default colors
library("gprofiler2")

datadir = "/project/biocompworkshop/Data_Vault/scRNAseq"
outdir = "/project/biocompworkshop/rshukla/scRNASeq_Results"

setwd(datadir)


