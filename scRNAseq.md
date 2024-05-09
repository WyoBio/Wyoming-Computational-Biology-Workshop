### Quality assessment

We are going to begin our single-cell analysis by loading in the output from CellRanger. We will load in our different samples, create a Seurat object with them, and take a look at the quality of the cells.

The general steps to preprocessing your single-cell data with Seurat:

-   Create a Seurat object
-   Filter low-quality cells
-   Merge samples
-   Normalize counts
-   Find variable features
-   Scale data
-   Determine PCs for Clustering
-   Clustering -\> FindNeighbors, FindClusters, RunUMAP

#### Load in Data

Seurat is the pivotal package for this stage, offering a distinct data structure and a suite of tools tailored for quality control, analysis, and exploration of single-cell RNA sequencing data. Widely acclaimed and recognized as a benchmark tool, Seurat stands out as the standard choice for single-cell RNA analysis.

```{R}
library("Seurat") 
library("ggplot2") # creating plots
library("cowplot") # add-on to ggplot, we use the plot_grid function to put multiple violin plots in one image file
library("dplyr")   # a set of functions for data frame manipulation -
library("Matrix")  # a set of functions for operating on matrices
library("viridis") # color maps for graphs that are more readable than default colors
library("gprofiler2") 
```

Let's ensure that our workspace is set up according to our preferences. We'll designate our working directory for saving files. It's crucial to confirm your current directory and ensure that you've selected the correct directory for saving your files.

```{R}
getwd()                         
dir()                           
setwd("/project/biocompworkshop/rshukla/scRNAseq")        
getwd()                        
outdir="/project/biocompworkshop/rshukla/scRNAseq/outdir_single_cell_rna" 
dir.create(outdir)
```

We've already loaded our single-cell data files into the session. Let's now examine them:

```{R}
list.files("/project/biocompworkshop/rshukla/scRNAseq")
```

You should see six h5 files:\
- `Rep1_ICB-sample_filtered_feature_bc_matrix.h5`\
- `Rep1_ICBdT-sample_filtered_feature_bc_matrix.h5`\
- `Rep3_ICB-sample_filtered_feature_bc_matrix.h5`\
- `Rep3_ICBdT-sample_filtered_feature_bc_matrix.h5`\
- `Rep5_ICB-sample_filtered_feature_bc_matrix.h5`\
- `Rep5_ICBdT-sample_filtered_feature_bc_matrix.h5`

#### Loading in one sample
