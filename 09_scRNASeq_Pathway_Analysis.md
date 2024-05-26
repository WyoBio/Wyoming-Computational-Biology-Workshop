# Gene set enrichment and pathway analysis

After performing differential expression analysis and identifying a list of interesting genes, the next step is often enrichment or pathway analysis. Broadly, enrichment analyses can be categorized into two types: overrepresentation analysis and gene set enrichment analysis (GSEA).

**Overrepresentation Analysis:** This method takes a list of significantly differentially expressed (DE) genes and determines if these genes are known to be differentially regulated within a specific pathway or gene set. It is particularly useful when we have a set of highly differentially expressed genes and want to identify the processes they may be involved in. Mathematically, it calculates a p-value using a hypergeometric distribution to see if a gene set from a database is significantly over-represented among our DE genes. Key points to note are:
- We can determine the list of genes used as inputs by setting thresholds for p-values and log2FC, which in turn define the gene list.

**Gene Set Enrichment Analysis (GSEA):** GSEA uses a list of genes and their corresponding fold change values as inputs. Unlike overrepresentation analysis, GSEA includes all genes without applying filters based on log2FC or p-values. GSEA is useful for detecting incremental changes in gene expression that collectively impact a specific pathway. It ranks genes based on their 'enrichment scores' (ES), which measure the degree to which a set of genes is over-represented at the top or bottom of a list of genes ordered by their log2FC values.

**Databases:** An important aspect of any enrichment analysis is the choice of databases. The main pitfall to avoid is using multiple or broad databases, as this can lead to many spurious results. It is advisable to select reference databases based on their biological relevance to obtain meaningful and accurate results.

### Enrichment Analysis Using clusterProfiler and Enrichr

There are several tools available for enrichment analysis, and in this case, we will use clusterProfiler. This tool supports both overrepresentation and GSEA analyses, is widely used in the field, and has numerous helpful tutorials and resources.

Additionally, we will utilize the web tool Enrichr for some of our analyses.

We will begin by examining the Epcam positive clusters identified in the Differential Expression section. First, load the necessary R libraries and read in the DE file generated previously. Recall that this file was created using the `FindMarkers` function in Seurat, with `ident.1` set to  `cluster 9` and `ident.2` set to `cluster 12`. Thus, we are comparing `cluster 9` to `cluster 12`, where positive `log2FC` values indicate genes upregulated in `cluster 9` or downregulated in `cluster 12`, and negative `log2FC` values indicate the opposite.

```R
# Load R libraries
library("Seurat")
library("ggplot2")
library("cowplot")
library("dplyr")
library("clusterProfiler")
library("org.Mm.eg.db")
library("msigdbr")
library("DOSE")
library("stringr")
library("enrichplot")

# Read in the epithelial DE file
de_gsea_df <- read.csv('outdir_single_cell_rna/epithelial_de_gsea.tsv', sep = '\t')

head(de_gsea_df)
# Open this file in Rstudio and get a sense for the distribution of foldchange values and see if their p values are significant
# Alternatively try making a histogram of log2FC values using ggplot and the geom_histogram() function. 
ggplot(de_gsea_df, aes(avg_log2FC)) + geom_histogram()
# You can also go one step further and impose a p-value cutoff (say 0.01) and plot the distribution.
```
