# Visualization
### Setup
We will establish working directories, load the necessary R libraries, and import the DESeqDataSet object and the results table generated in the previous section into the R environment.

```R
# Set working directories
datadir = "/project/biocompworkshop/rshukla/DE_Results"
outdir = "/Users/ramshukla/NPRG Dropbox/Computaional Biology Workshop/Day2/Differential Expression/DE_Results"

# load R libraries
library(DESeq2)
library(data.table)
library(pheatmap)
library(gplots)
library(ggplot2)
library(ggrepel)

# Load in the DESeqDataSet object
dds <- readRDS('dds.rds')

# Load in the results object before shrinkage
res <- readRDS('res.rds')

# Load in the results object after shrinkage
resLFC <- readRDS('resLFC.rds')

# Load in the final results file with all sorted DE results
deGeneResultSorted <- fread('DE_all_genes_DESeq2.tsv')
```
### MA-plot before LFC shrinkage
Originally designed for assessing microarray expression data, MA-plots utilize the log ratio (M) and mean average (A) derived from scanned intensity measurements. Even in RNAseq differential expression (DE) experiments with two conditions, these plots retain their utility. They swiftly reveal the count of significantly differentially expressed genes, the balance between up- and down-regulated genes, and any outliers.

To interpret MA-plots effectively, remember that the Y axis (M) represents the log2 fold change between tested conditions, where a larger fold change signifies a more pronounced difference. On the X axis (A), we gauge gene read alignment, with higher positions indicating greater total aligned reads—an indicator of overall gene expression (albeit without considering gene length in raw read counts).

Utilizing DESeq2's built-in `plotMA` function, genes are color-coded based on significance thresholds. Notably, genes with elevated expression values and larger fold changes more frequently attain significance, aligning with expectations.

```R
# use DESeq2 built in MA-plot function
plotMA(res, ylim=c(-2, 2))
```
### MA-plot after LFC shrinkage
Upon executing DESeq2, we derived two outcomes—one incorporating log-fold change shrinkage and one without. In scenarios with genes registering low hits, we often encounter disproportionately large fold changes. Consider a scenario where one gene has 1 hit compared to another with 6 hits; this results in a 6x fold change. However, this substantial variance likely reflects noise rather than genuine biological signals.

By employing `plotMA` on our findings post log-fold change shrinkage algorithm implementation, we observe a more controlled depiction of this phenomenon.

```R
plotMA(resLFC, ylim=c(-2,2))
```
Given the focused nature of our dataset (specifically chr22 genes), the observed effect is quite nuanced. However, a comparison between the two plots reveals subtle shifts, particularly noticeable in the upper left and bottom left corners, where certain fold change values are converging toward 0.

### Inspecting gene counts across conditions
It's often beneficial to examine the normalized counts for a gene across our samples. DESeq2 offers a convenient built-in function for this purpose, operating on the dds object. Let's take SEPT3, identified as significantly higher in the UHR cohort in our DE analysis. This view allows us to assess the per-sample distribution of our corrected counts, promptly identifying any outliers within each group and facilitating further investigation if necessary.

```R
# view SEPT3 normalized counts
plotCounts(dds, gene='ENSG00000100167', intgroup = 'Condition')

# view PRAME normalized counts
plotCounts(dds, gene='ENSG00000185686', intgroup = 'Condition')
```
### Examining pairwise sample clustering
Examining the relatedness between samples can provide valuable insights into their overall similarity or dissimilarity. Although not included in DESeq2, a convenient library allows for the creation of a hierarchically clustered heatmap from our DESeq2 data. It's important to note that before conducting any distance calculations, the count data should be transformed using `vst()` or `rlog()`, which can be directly applied to the dds object.

```R
# We opt for rlog due to our relatively smaller gene set. For DE experiments with thousands of genes, it's advisable to utilize the vst() function.
rld <- rlog(dds, blind=F)
# # view the structure of this object
rld

# Compute sample distances using the `dist` function (defaulting to Euclidean distance metric).
# Extract the rlog-transformed data using `assay`.
# Transpose the data using `t`.
# Calculate distance values using `dist`.
# The distance is computed for each vector of sample gene values, comparing all samples in a pairwise manner.

# view the first few lines of raw data
head(assay(dds))

# see the rlog transformed data
head(assay(rld))

# see the impact of transposing the matrix
t(assay(rld))[1:6,1:5]

# see the distance values
dist(t(assay(rld)))

# put it all together and store the result
sampleDists <- dist(t(assay(rld)))

# convert the distance result to a matrix
sampleDistMatrix <- as.matrix(sampleDists)

# view the distance numbers directly in the pairwise distance matrix
head(sampleDistMatrix)

# construct clustered heatmap, important to use the computed sample distances for clustering
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists)
```
We can opt for a similarity metric like Pearson correlation instead of a distance metric. There are various correlation and distance options available:

- Correlation metrics: "pearson", "kendall", "spearman"
- Distance metrics: "euclidean", "maximum", "manhattan", "canberra", "binary", or "minkowski"

```R
sampleCorrs <- cor(assay(rld), method="pearson")
sampleCorrMatrix <- as.matrix(sampleCorrs)
head(sampleCorrMatrix)

pheatmap(sampleCorrMatrix)
```

Instead of summarizing all gene count data into a distance metric for each sample, you can gain a similar understanding of the pattern by visualizing all genes simultaneously.

```R
pheatmap(mat=t(assay(rld)), show_colnames = FALSE)
```

### Additional R Visualization for Differential Expression
At times, you might want to customize and manipulate expression estimates in R in a more flexible manner. This section covers an in-depth tutorial on visualizing your results in R and conducting "old-school" (non-DESeq2) visualization of your data.

We will begin by generating normalized data from the `dds` object and then proceed with our analysis. Working directly on count data is not feasible because it lacks the appropriate scaling and normalization required for accurate comparisons. The steps for creating the dds object will be reiterated.

```R
setwd(datadir)
rawdata=read.table("featurecounts.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1)
colnames(rawdata)
rawdata <- rawdata[,-c(1,2,3,4,5,6,10)]
# Extract and edit column names
# Sample column names
column_names <- colnames(rawdata)

# Extract "HBR_RepX" and "UHR_RepX" from column names
extracted_names <- gsub(".*\\.([A-Z]+)_Rep([0-9]+)\\.bam", "\\1_Rep\\2", column_names)

colnames(rawdata) <- extracted_names

# Check dimensions
dim(rawdata)

# Check dimensions again to see effect of filtering
dim(rawdata)

metaData <- data.frame('Condition'=c('UHR', 'UHR', 'UHR', 'HBR', 'HBR', 'HBR'))

metaData$Condition <- factor(metaData$Condition, levels=c('HBR', 'UHR'))

rownames(metaData) <- colnames(rawdata)

head(metaData)

all(rownames(metaData) == colnames(rawdata))

dds <- DESeqDataSetFromMatrix(countData = rawdata, colData = metaData, design = ~Condition)

dds <- estimateSizeFactors(dds)

sizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)

normalized_counts <- as.data.frame(normalized_counts)

gene_expression <- normalized_counts
```
Load other files and the results generated from the previous section.

```R
setwd(outdir)

# Load gene name mapping file
gene_names=read.table("ENSG_ID2Name.txt", header=TRUE, stringsAsFactors=FALSE)
colnames(gene_names)=c("gene_id","gene_name")

# Load DE results from the DESeq2 pipeline
results_genes <-read.table("DE_all_genes_DESeq2.tsv", sep="\t", header=T, stringsAsFactors = F)
```
Let's take a quick look at the imported data to get a sense of its content.

```R
#### Working with 'dataframes'
# View the first five rows of data (all columns).
head(gene_expression)

# View the column names
colnames(gene_expression)

# View the row names
row.names(gene_expression)

# Determine the dimensions of the dataframe. 
dim(gene_expression)

# Get the first 3 rows of data and a selection of columns
gene_expression[1:3,c(1:3,6)]

# Do the same thing, but using the column names instead of numbers
gene_expression[1:3, c("HBR_Rep1","HBR_Rep2","HBR_Rep3","UHR_Rep3")]

# Now, exlore the differential expression (DESeq2 results) 
head(results_genes)
dim(results_genes)

# Assign some colors for use later.  You can specify color by RGB, Hex code, or name
# To get a list of color names:
colours()
data_colors=c("tomato1","tomato2","tomato3","royalblue1","royalblue2","royalblue3")
```
The following code blocks generate various plots using the dataset.

### Plot #1: Visualizing Normalized Count (NC) Value Range and Distribution across 6 Libraries

```R
# Create boxplots for this purpose
# Display on a log2 scale and set a minimum non-zero value to avoid log2(0)
min_nonzero=1

# Set the columns for finding TPM and create shorter names for figures
data_columns=c(1:6)
short_names=c("HBR_1","HBR_2","HBR_3","UHR_1","UHR_2","UHR_3")

boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(TPM)", main="Distribution of NC for all 6 libraries")
# Note that the bold horizontal line on each boxplot is the median
```
### Plot #2: Assessing Technical Replicates' Reproducibility

```R
# Transform the data by adding a small arbitrary value and then converting it to the log2 scale to avoid issues with log2(0).
x = gene_expression[,"UHR_Rep1"]
y = gene_expression[,"UHR_Rep2"]

plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="TPM (UHR, Replicate 1)", ylab="TPM (UHR, Replicate 2)", main="Comparison of expression values for a pair of replicates")

# Add a straight line of slope 1, and intercept 0
abline(a=0,b=1)

# Calculate the correlation coefficient and display in a legend
rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")
```
### Plot #3: Transforming Scatter Plots into Density Scatter Plots

```R
colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
smoothScatter(x=log2(x+min_nonzero), y=log2(y+min_nonzero), xlab="TPM (UHR, Replicate 1)", ylab="TPM (UHR, Replicate 2)", main="Comparison of expression values for a pair of replicates", colramp=colors, nbin=200)
```

### Plot #4: Consolidated Scatter Plots of Replicate Sets

```R
# Generate an R plot using a function that accepts two libraries for comparison and a plot name as input.
plotCor = function(lib1, lib2, name){
  x=gene_expression[,lib1]
  y=gene_expression[,lib2]
  zero_count = length(which(x==0)) + length(which(y==0))
  colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  smoothScatter(x=log2(x+min_nonzero), y=log2(y+min_nonzero), xlab=lib1, ylab=lib2, main=name, colramp=colors, nbin=275)
  abline(a=0,b=1)
  rs=cor(x,y, method="pearson")^2
  legend_text = c(paste("R squared = ", round(rs, digits=3), sep=""), paste("Zero count = ", zero_count, sep=""))
  legend("topleft", legend_text, lwd=c(1,NA), col="black", bg="white", cex=0.8)
}

# Now make a call to our custom function created above, once for each library comparison
par(mfrow=c(1,3))
plotCor("UHR_Rep1", "UHR_Rep2", "UHR_1 vs UHR_2")
plotCor("UHR_Rep2", "UHR_Rep3", "UHR_2 vs UHR_3")
plotCor("UHR_Rep1", "UHR_Rep3", "UHR_1 vs UHR_3")

#### Compare the correlation between all replicates
# Is the observed pattern consistent across all eight libraries, with replicates showing the most similarity followed by tumor versus normal samples?

# Calculate the NC sum for all 6 libraries
gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)

# Identify genes with a total NC sum of at least 5 - we'll exclude genes with extremely low expression levels.
i = which(gene_expression[,"sum"] > 5)

# Calculate the correlation between all pairs of data
r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")
r
```

### Plot #5: Converting Correlation to Distance for Multi-Dimensional Scaling

```R
# Generate 2D coordinates for plotting points representing each library.
# Libraries with similar expression patterns, indicated by high correlation, should cluster together.
# What grouping pattern do we anticipate, considering the library types (technical replicates, biological replicates, tumor/normal)?
d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
par(mfrow=c(1,1))
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes)", xlim=c(-0.12,0.12), ylim=c(-0.12,0.12))
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], short_names, col=data_colors)
```

### Plot #6: Distribution of differential expression values as a histogram

```R
# Display only those results that are significant according to DESeq2 (loaded above)
sig=which(results_genes$pvalue<0.05)
hist(results_genes[sig,"log2FoldChange"], breaks=50, col="seagreen", xlab="log2(Fold change) UHR vs HBR", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topleft", "Fold-change > 4", lwd=2, lty=2)
```

### Plot #7: Mean Expression Values Comparison with Significant Differential Expression

```R
gene_expression[,"HBR_mean"]=apply(gene_expression[,c(1:3)], 1, mean)
gene_expression[,"UHR_mean"]=apply(gene_expression[,c(4:6)], 1, mean)

x=log2(gene_expression[,"UHR_mean"]+min_nonzero)
y=log2(gene_expression[,"HBR_mean"]+min_nonzero)
plot(x=x, y=y, pch=16, cex=0.25, xlab="UHR TPM (log2)", ylab="HBR TPM (log2)", main="UHR vs HBR TPMs")
abline(a=0, b=1)
xsig=x[sig]
ysig=y[sig]
points(x=xsig, y=ysig, col="magenta", pch=16, cex=0.5)
legend("topleft", "Significant", col="magenta", pch=16)

# Get the gene symbols for the top N (according to corrected p-value) and display them on the plot
# topn = order(abs(results_genes[sig,"log2FoldChange"]), decreasing=TRUE)[1:25]
topn = order(results_genes[sig,"padj"])[1:25]
text(x[topn], y[topn], results_genes[topn,"Symbol"], col="black", cex=0.75, srt=45)
```

### Plot #8: Heatmap for Expression Differences Across Six Samples

```R
# Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

# Create a subset of significant genes with p-value<0.05 and log2 fold-change >= 2
sigpi = which(results_genes[,"pvalue"]<0.05)
sigp = results_genes[sigpi,]
sigfc = which(abs(sigp[,"log2FoldChange"]) >= 2)
sigDE = sigp[sigfc,]

main_title="sig DE Genes"
par(cex.main=0.8)
sigDE_genes=sigDE[,"ensemblID"]
sigDE_genenames=sigDE[,"Symbol"]

data=log2(as.matrix(gene_expression[as.vector(sigDE_genes),data_columns])+1)
heatmap.2(data, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(10,4), Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, cexRow=0.3, cexCol=1, labRow=sigDE_genenames,col=rev(heat.colors(75)))
```

### Plot #9: Volcano plot

```R
# Set differential expression status for each gene - default all genes to "no change"
results_genes$diffexpressed <- "No"

# if log2Foldchange > 2 and pvalue < 0.05, set as "Up regulated"
results_genes$diffexpressed[results_genes$log2FoldChange >= 2 & results_genes$pvalue < 0.05] <- "Higher in UHR"

# if log2Foldchange < -2 and pvalue < 0.05, set as "Down regulated"
results_genes$diffexpressed[results_genes$log2FoldChange <= -2 & results_genes$pvalue < 0.05] <- "Higher in HBR"

# write the gene names of those significantly upregulated/downregulated to a new column
results_genes$gene_label <- NA
results_genes$gene_label[results_genes$diffexpressed != "No"] <- results_genes$Symbol[results_genes$diffexpressed != "No"]

ggplot(data=results_genes[results_genes$diffexpressed != "No",], aes(x=log2FoldChange, y=-log10(pvalue), label=gene_label, color = diffexpressed)) +
  xlab("log2Foldchange") +
  scale_color_manual(name = "Differentially expressed", values=c("blue", "red")) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  geom_point(data = results_genes[results_genes$diffexpressed == "No",], aes(x=log2FoldChange, y=-log10(pvalue)), colour = "black")
```

