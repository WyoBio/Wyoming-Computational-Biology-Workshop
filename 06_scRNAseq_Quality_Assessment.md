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

