# Cell Type Annotation

### What's the importance of labeling cell types in single-cell RNAseq data?
Assigning cell type annotations in single-cell gene expression data is a dynamic area of research. This task is challenging due to two main factors:

- Gene expression levels often exist on a continuum rather than in discrete categories.
- Variations in gene expression don't always equate to functional differences between cells (Pasquini et al.).

The field is rapidly evolving with the introduction of various annotation tools and the development of multiple scRNA-seq databases (Wang et al.).

One effective approach for cell type annotation is utilizing web resources. These tools are user-friendly and do not require advanced scripting or programming skills. Notable examples include Azimuth, Tabula Sapiens, and SciBet. Researchers often employ multiple tools for cell type annotation to ensure accuracy. Additionally, it is advisable to validate annotations through experiments, statistical analyses, or consultations with subject matter experts.

Methods for assigning cell types typically fall into one of two categories:

- **Databases with Marker Genes for Manual Annotation**: These databases contain marker genes for various cell types (primarily human and mouse). By using these marker genes, you can manually annotate the cells in your datasets. For example, in the 10x Genomics Loupe Browser, you can identify top differentially expressed genes for each cluster. You can then search these genes in the database to determine if they are markers for specific cell types. To ensure accuracy, you may need to search multiple top genes in each cluster. This process can be laborious and time-consuming, depending on the dataset's complexity and your prior knowledge of the cell types.

- **Resources for Automated Cell Type Annotation (Reference-Based)**: These community-developed tools automatically annotate cells by comparing new data with existing references. This article highlights several web tools that do not require programming skills. A major limitation of these tools is that the quality of the results heavily depends on the quality of the pre-annotated reference datasets.

### Setup

```R
library(SingleR)
library(celldex)
library(Seurat)
library(cowplot)

merged <- readRDS("outdir_single_cell_rna/rep135_clustered.rds")
```


