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

datadir = "/project/biocompworkshop/rshukla/scRNASeq_Results"

setwd("/project/biocompworkshop/rshukla/scRNASeq_Results")

de_gsea_df <- read.csv('epithelial_de_gsea.tsv', sep = '\t')

head(de_gsea_df)

# Open this file in Rstudio and get a sense for the distribution of foldchange values and see if their p values are significant
# Alternatively try making a histogram of log2FC values using ggplot and the geom_histogram() function. 

ggplot(de_gsea_df, aes(avg_log2FC)) + geom_histogram()

```

### Overrepresentation analysis

You may notice that we have several genes with significant fold change values. While fold change does not affect overrepresentation analysis, it can help set thresholds for selecting genes. Since many genes have fold changes greater than +2 or less than -2, we can use this as our cutoff. Additionally, we will set an adjusted p-value cutoff of 0.01. For the overrepresentation analysis, we will filter de_gsea_df based on these log2FC and p-value thresholds to get the list of genes for our analysis.

```R
# Filter de_gsea_df by subsetting it to only include genes that are significantly DE (pval<0.01) and their absolute log2FC is > 2.
# The abs(de_gsea_df$avg_log2FC) ensures that we keep both the up and downregulated genes

overrep_df <- de_gsea_df[de_gsea_df$p_val_adj < 0.01 & abs(de_gsea_df$avg_log2FC) > 2,] 
overrep_gene_list <- rownames(overrep_df)
```

Next, we will establish our reference. Although `clusterProfiler` typically uses the msigdb reference by default, we will demonstrate how to download a mouse-specific cell type signature reference geneset from msigdb for your analysis. Specifically, we will utilize the M8 geneset from the msigdb mouse collections. By clicking on the `Gene Symbols` link on the right, you can download the dataset and upload it to your workspace. These files are in a `gmt` (gene matrix transposed) format and can be read in using the built-in R function, `read.gmt`.

Once we have the reference data loaded, we will employ the `enricher` function from the `clusterProfiler` library for the overrepresentation analysis. This function requires inputs such as the DE gene list, the reference database, the statistical method for p-value adjustment, and a p-value cutoff threshold. The enricher function generates an overrepresentation R object that can be used in visualization functions like `barplot()` and `dotplot()` to create typical pathway analysis figures. Additionally, we can utilize the web tool Enrichr for a rapid analysis across multiple databases. For this purpose, we will save the gene list used for the overrepresentation analysis to a TSV file.

```R
# Accessing the Tabula Muris GMT File
# You can download it from: https://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp?targetSpeciesDB=Mouse#M8
# Note that you'll need to create a user ID and password to access the file. For convenience, I've uploaded the files to the data vault.

msigdb_m8 <- read.gmt('/project/biocompworkshop/Data_Vault/m8.all.v2023.2.Mm.symbols.gmt')

# Click on the dataframe in RStudio to see how it's formatted- we have 2 columns, #the first with the genesets, and the other with genes that are in that geneset.
# Try to determine how many different pathways are in this database

overrep_msigdb_m8 <- enricher(gene = overrep_gene_list, TERM2GENE = msigdb_m8, pAdjustMethod = "BH", pvalueCutoff = 0.05)

# Visualize data using the barplot and dotplot functions

barplot(overrep_msigdb_m8, showCategory = 10)
dotplot(overrep_msigdb_m8, showCategory = 10)

# Save overrep_gene_list to a tsv file (overrep_gene_list is our list of genes and 
# File is the name we want the file to have when it's saved. 
# The remaining arguments are optional- row.names=FALSE stops R from adding numbers (effectively an S.No column), 
# col.names gives our single column TSV a column name, 
# and quote=FALSE ensures the genes don't have quotes around them which is the default way R saves string values to a TSV)

write.table(x = overrep_gene_list, file = 'outdir_single_cell_rna/epithelial_overrep_gene_list.tsv', row.names = FALSE, col.names = 'overrep_genes', quote=FALSE)
```

For the Enrichr webtool analysis, we will begin by opening the TSV file in our RStudio session. We'll then copy the genes from the file and paste them directly into the text box provided on the right side of the webtool interface. Upon submission, the webtool will generate multiple barplots depicting enriched pathways. Feel free to explore the results by clicking around the plots. To compare these results with those generated in R, navigate to the `Cell Types` tab at the top and look for `Tabula Muris`.

A crucial aspect of a robust overrepresentation analysis involves leveraging biological expertise alongside the identified pathways to generate hypotheses. Not every pathway displayed in the plots may be relevant, but considering our knowledge of bladder cancer (for this dataset), we can infer that basal and luminal bladder cancers share similar expression profiles with basal and luminal breast cancers. Therefore, the overrepresentation analysis indicating genesets such as 'Tabula Muris senis mammary gland basal cell ageing' and 'Tabula muris senis mammary gland luminal epithelial cell of mammary gland ageing' could suggest that the differences observed in unsupervised clusters 9 and 12 may originate from basal and luminal cells.

To delve deeper into this hypothesis, we can compile a list of basal and luminal markers from literature, calculate a combined score for those genes using Seurat's `AddModuleScore` function, and ascertain if the clusters segregate into basal and luminal categories. For now, we will utilize the same markers defined in the original manuscript of this dataset.

```R
# Define lists of marker genes

basal_markers <- c('Cd44', 'Krt14', 'Krt5', 'Krt16', 'Krt6a')
luminal_markers <- c('Cd24a', 'Erbb2', 'Erbb3', 'Foxa1', 'Gata3', 'Gpx2', 'Krt18', 'Krt19', 'Krt7', 'Krt8', 'Upk1a')

# Read in the seurat object if it isn't loaded in your R session

merged <- readRDS('preprocessed_object.rds')

# Use AddModuleScore to calculate a single score that summarizes the gene expression for each list of markers

merged <- AddModuleScore(merged, features=list(basal_markers), name='basal_markers_score')
merged <- AddModuleScore(merged, features=list(luminal_markers), name='luminal_markers_score')

# Visualize these scores using FeaturePlot and VlnPlots
FeaturePlot(merged, features=c('basal_markers_score1', 'luminal_markers_score1'))
VlnPlot(merged, features=c('basal_markers_score1', 'luminal_markers_score1'), group.by = 'seurat_clusters_res0.8', pt.size=0)
```

Interesting! This analysis suggests that cluster 12 consists of basal epithelial cells, while cluster 9 consists of luminal epithelial cells. Next, let's use GSEA to identify distinct biological processes between these clusters.

### GSEA Analysis

For GSEA, we need to start by creating a named vector where the values are the log fold change values, and the names are the gene names. GSEA analysis relies on identifying incremental gene expression changes (not just those that are statistically significant), so we will use our original unfiltered dataframe to get these values. This will be used as input to the `gseGO` function in the clusterProfiler library, which uses gene ontology for GSEA analysis.

The other parameters for the function include:
- `OrgDb = org.Mm.eg.db`, the organism database from where all the pathways’ genesets will be determined.
- `ont = "ALL"`, specifies the subontologies, with possible options being BP (Biological Process), MF (Molecular Function), CC (Cellular Compartment), or ALL.
- `keyType = "SYMBOL"` tells `gseGO` that the genes in our named vector are gene symbols as opposed to Entrez IDs or Ensembl IDs.
- `pAdjustMethod="BH"` and `pvalueCutoff=0.05` specify the p-value adjustment statistical method to use and the corresponding cutoff.

```R
# Read in the epithelial de df we generated previously
de_gsea_df <- read.csv('outdir_single_cell_rna/epithelial_de_gsea.tsv', sep = '\t')

# Get all the foldchange values from the original dataframe

gene_list <- de_gsea_df$avg_log2FC

# Set names for this vector to gene names

names(gene_list) <- rownames(de_gsea_df)

# Sort list in descending order of log2FC values as that is required by the gseGO function

gene_list = sort(gene_list, decreasing = TRUE)

# Now we can run the gseGO function

gse <- gseGO(geneList=gene_list, 
             OrgDb = org.Mm.eg.db,
             ont ="ALL",              
             keyType = "SYMBOL", 
             pAdjustMethod = "BH",             
             pvalueCutoff = 0.05)

# Explore the gse object by opening it in RStudio. 
# It basically has a record of all the parameters and inputs used for the function, 
# along with a results dataframe.
# We can pull this result dataframe out to view it in more detail

gse_result <- gse@result
```

Upon examining the dataframe, you may observe approximately 900 rows. Similar to the overrepresentation analysis, not all of these rows are likely to hold significant biological relevance. While you can manually sift through the results to identify meaningful pathways, for efficiency, we will focus on subsets within the `gse` object that contain pathways with the term "epithelial" in their names for plotting purposes.

To achieve this, we will first determine the indices of rows containing these values using R's `which` and `grepl` functions. Subsequently, we will subset the results dataframe in the gse object to include only those rows. Notably, all these genesets exhibit negative enrichment scores (indicating downregulation in our putative luminal cell cluster). To facilitate plotting, we will add one index with a positive enrichment score. The plotting functions we will employ include `dotplot`, `cnetplot`, and `heatplot`.

```R
# Start by grabbing the indices we'll need to subset the `gse` object
# Epithelial indices from gse object using which and grepl

epithelial_indices <- which(grepl("epithelial", gse@result$Description))

# Index for the most positive enrichment score using which.max()
positive_index <- which.max(gse@result$enrichmentScore)

# Concatenate these indices to get the list of indices that will be used to subset the gsea object

subset_indices <- c(positive_index,epithelial_indices)

# Now create a new gse_epithelial object and subset the gse object to these indices

gse_epithelial <- gse
gse_epithelial@result <- gse_epithelial@result[subset_indices,]

# Plot!
# dotplot - splitting by 'sign' and facet_grid together allow us to separate activated and suppressed pathways

dotplot(gse_epithelial, showCategory=20, split=".sign") + facet_grid(.~.sign) 

# Heatplot - allows us to see the genes that are being considered 
# for each of the pathways/genesets and their corresponding fold change

heatplot(gse_epithelial, foldChange=gene_list)

# cnetplot - allows us to see the genes along with the various 
# Pathways/genesets and how they related to each other

cnetplot(gse_epithelial, foldChange=gene_list)
```

From these findings, we can infer that cluster 9 (putative luminal cells) exhibit reduced expression levels in several pathways associated with epithelial cell proliferation compared to cluster 12 (putative basal cells).

As an additional exercise, let's conduct an overrepresentation and/or GSEA analysis for the DE analysis conducted on CD8 T cells. Since we didn't save the DE results file earlier, you can download the DE file from the provided link first and utilize it for your analysis.

```R
download.file(url = 'http://genomedata.org/cri-workshop/reference_files/cd8tcells_de_gsea.tsv',
              destfile = 'outdir_single_cell_rna/cd8tcells_de_gsea.tsv')
```

Hint:

The comparison between ICB and ICBdT in CD8 T cells might present subtler differences compared to the luminal vs. basal epithelial cell comparison. Consider adjusting the fold-change cutoff for differentially expressed genes, for example, to 0.5.

The M8 cell type signature gene sets may not be the most pertinent choice. Consider downloading a different set of gene sets, such as the MH hallmark gene sets. You can visit the msigdb mouse collections, download the 'Gene Symbols' GMT file for the desired gene sets, upload it to your posit cloud environment, and adjust the relevant R commands to load the gmt file accordingly.

