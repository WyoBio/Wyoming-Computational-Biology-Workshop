# Pathway analysis
### Setup
In this section, we'll utilize the GAGE (Generally Applicable Gene-set Enrichment tool) tool in R to assess significantly enriched gene sets among those exhibiting significant "up" and "down" regulation in our UHR vs HBR gene expression analysis. The question we aim to answer is: Do genes associated with brain cell types and processes show enrichment among the DE genes with notable expression differences between UHR and HBR samples?

GAGE, a widely used bioconductor tool, conducts gene-set enrichment and pathway analysis regardless of sample sizes, experimental designs, or assay platforms, suitable for microarray and RNAseq datasets. Here, we'll apply GAGE using gene sets from the "Gene Ontology" (GO) and MSigDB databases for pathway analysis.

```R
# Load R libraries 
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
library(gage)
```
- Importing Differential Expression Results for Pathway Analysis
```R
# Define and set working dir paths
datadir = "/Users/ramshukla/NPRG Dropbox/Computaional Biology Workshop/Day2/Differential Expression/DE_Results"
setwd(datadir)

# Load in the DE results file with only significant genes
DE_genes <-read.table("DE_sig_genes_DESeq2.tsv", sep="\t", header=T, stringsAsFactors = F)
```
### Preparing Gene Set Databases
- In order to conduct our pathway analysis effectively, we require a comprehensive list of pathways and their corresponding genes. Numerous databases house collections of genes or gene sets that can help determine the functional relationships among mutated or differentially expressed genes. Some key resources include GO, KEGG, MSigDB, and WikiPathways. For this task, our focus will be on investigating GO and MSigDB.

- The GAGE package offers a real-time querying function for GO called go.gsets(). This function, when provided with a species argument, generates a list of gene sets along with useful metadata for easy subsetting. Understanding GO's three gene ontologies—Biological Process, Molecular Function, and Cellular Component—will be beneficial for our analysis.

- While GAGE lacks a similar tool for exploring MSigDB gene sets, MSigDB conveniently provides a downloadable .gmt file containing all gene sets. This file format can be effortlessly read into GAGE using the readList() function. Upon reviewing MSigDB, you'll find eight distinct gene set collections, each with unique features. For this exercise, we'll utilize the C8 - cell type signature gene sets collection, which comprises gene sets containing cluster markers from single-cell sequencing studies of human tissue.

```R
# Set up go database
go.hs <- go.gsets(species="human")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF]
go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC]

# Here we will read in an MSigDB gene set that was selected for this exercise. 
c8 <-"http://genomedata.org/rnaseq-tutorial/c8.all.v7.2.entrez.gmt"
all_cell_types <-readList(c8)
```
### Annotating genes
Alright, we have our set of differentially expressed genes and our gene sets. Yet, if you examine one of the objects containing these gene sets, you'll find a sequence of integers within each set. These integers represent Entrez gene identifiers. However, do we have similar information in our list of DE genes? Currently, no. Our prior results utilize Ensembl IDs as gene identifiers. To proceed with pathway analysis, we must convert our gene identifiers to the format used in GO and MSigDB gene sets.

Thankfully, Bioconductor maintains comprehensive genome-wide annotation data for various species, accessible via the OrgDb bioc view command. This simplifies the conversion of gene identifiers. Below, we employ the mapIds() function to query the OrganismDb object for Entrez IDs based on Ensembl IDs. Given potential one-to-many mappings, we use multiVals="first" to retrieve only the first identifier. Alternatively, multiVals="asNA" can be used to disregard one-to-many mappings.

```R
DE_genes$entrez <- mapIds(org.Hs.eg.db, keys=DE_genes$ensemblID, column="ENTREZID", keytype="ENSEMBL", multiVals="first")
```
### Refining and Mapping Identifiers
Upon completing the annotation process, you may observe that some of our Ensembl gene IDs were not mapped to an Entrez gene ID. Why does this occur? This issue delves into the complex realm of gene definition and annotation, reflecting differences in how resources annotate the human genome. Consequently, some discrepancies are anticipated. In the subsequent steps, we will address this by initially removing the ERCC spike-in genes and then utilizing an alternate identifier for further mapping.

```R
# Remove spike-in
DE_genes_clean <- DE_genes[!grepl("ERCC", DE_genes$ensemblID),]

# Just so we know what we have removed 
ERCC_gene_count <-nrow(DE_genes[grepl("ERCC", DE_genes$ensemblID),])
ERCC_gene_count

# Deal with genes that we do not have an Entrez ID for 
missing_ensembl_key<-DE_genes_clean[is.na(DE_genes_clean$entrez),]
DE_genes_clean <-DE_genes_clean[!DE_genes_clean$ensemblID %in% missing_ensembl_key$ensemblID,]

# Try mapping using a different key
missing_ensembl_key$entrez <- mapIds(org.Hs.eg.db, keys=missing_ensembl_key$Symbol, column="ENTREZID", keytype="SYMBOL", multiVal='first')

# Remove remaining genes 
missing_ensembl_key_update <- missing_ensembl_key[!is.na(missing_ensembl_key$entrez),]

# Create a Final Gene list of all genes where we were able to find an Entrez ID (using two approaches)
DE_genes_clean <-rbind(DE_genes_clean, missing_ensembl_key_update)
```
### Preparing DESeq2 Results for GAGE Analysis
Last step! Let’s format the differential expression results into a format suitable for the GAGE package. Basically this means obtaining the log2 fold change values and assigning entrez gene identifiers to these values.

```R
# grab the log fold changes for everything
De_gene.fc <- DE_genes_clean$log2FoldChange

# set the name for each row to be the Entrez Gene ID
names(De_gene.fc) <- DE_genes_clean$entrez
```
### Running pathway analysis
Now, we can utilize the `gage()` function to extract significantly perturbed pathways from our differential expression experiment.

Abbreviations:
- "bp": Biological Process
- "mf": Molecular Function
- "cc": Cellular Component

These abbreviations represent the three primary categories of gene ontology terms/annotations mentioned earlier.

```R
# Run GAGE
# go 
fc.go.bp.p <- gage(De_gene.fc, gsets = go.bp.gs)
fc.go.mf.p <- gage(De_gene.fc, gsets = go.mf.gs)
fc.go.cc.p <- gage(De_gene.fc, gsets = go.cc.gs)

# msigdb
fc.c8.p <- gage(De_gene.fc, gsets =all_cell_types)


###  Convert to dataframes 
# Results for testing for GO terms which are up-regulated
fc.go.bp.p.up <- as.data.frame(fc.go.bp.p$greater)
fc.go.mf.p.up <- as.data.frame(fc.go.mf.p$greater)
fc.go.cc.p.up <- as.data.frame(fc.go.cc.p$greater)

# Results for testing for GO terms which are down-regulated
fc.go.bp.p.down <- as.data.frame(fc.go.bp.p$less)
fc.go.mf.p.down <- as.data.frame(fc.go.mf.p$less)
fc.go.cc.p.down <- as.data.frame(fc.go.cc.p$less)

# Results for testing for MSigDB C8 gene sets which are up-regulated
fc.c8.p.up <- as.data.frame(fc.c8.p$greater)

# Results for testing for MSigDB C8 gene sets which are down-regulated
fc.c8.p.down <- as.data.frame(fc.c8.p$less)
```
### Explore significant results
We now have results with corresponding p-values.

Let's dive into understanding what "up-regulated" or "down-regulated" signifies in the context of our UHR vs HBR comparison. It might be beneficial to open and review the data in your DE_genes_DESeq2.tsv file.

Take a look at the cellular process results from our GO analysis. Do these results align with your expectations?

```R
# Try doing something like this to find some significant results:
# View the top 20 significantly up- or down-regulated GO terms from the Cellular Component Ontology
head(fc.go.cc.p.up[order(fc.go.cc.p.up$p.val),], n=20)
head(fc.go.cc.p.down[order(fc.go.cc.p.down$p.val),], n=20)

# You can do the same thing with your results from MSigDB
head(fc.c8.p.up)
head(fc.c8.p.down)
```
### More exploration
Now, let's transition away from R and delve deeper into exploring our results locally. For the rest of the exercise, our focus will be on the results obtained from GO analysis. We'll utilize an online tool to visualize the interrelationships among the GO terms we uncovered.

```R
write.table(fc.go.cc.p.up, "fc.go.cc.p.up.tsv", quote = F, sep = "\t", col.names = T, row.names = T)
write.table(fc.go.cc.p.down, "fc.go.cc.p.down.tsv", quote = F, sep = "\t", col.names = T, row.names = T)
```
