## Setup

### Load R libraries we will use in this section

```R
library(DESeq2)
library(data.table)
library(apeglm)
```
### Define working dir paths

```R
datadir = "/Users/ramshukla/NPRG Dropbox/Computaional Biology Workshop/Day2/Differential Expression"
outdir = "/Users/ramshukla/NPRG Dropbox/Computaional Biology Workshop/Day2/Differential Expression/DE_Results"
```
### Set working directory to datadir

```R
setwd(datadir)
```
### Read in the RNAseq read counts for each gene (produced by featurecounts)

```R
rawdata=read.table("featurecounts.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1)
colnames(rawdata)
rawdata <- rawdata[,-c(1,2,3,4,5,6,10)] # Remove columns which are not required
```
### Extract and edit column names
By default, featurecount outputs the sample name as the file name along with the file path. We will modify this to ensure the correct sample names are used.

```R
column_names <- colnames(rawdata)
extracted_names <- gsub(".*\\.([A-Z]+)_Rep([0-9]+)\\.bam", "\\1_Rep\\2", column_names) # Extract "HBR_RepX" and "UHR_RepX" from column names
colnames(rawdata) <- extracted_names
dim(rawdata) # Check dimensions
```
### Filter raw counts
Before running DESeq2 or any differential expression analysis, it's beneficial to pre-filter the data. This process offers computational advantages, such as reducing the memory size of R objects, which in turn speeds up DESeq2's processing. Removing "low-quality" data not only streamlines the statistical tests but also minimizes the need for multiple test corrections, potentially highlighting more significant genes.

```R
* # Require at least 1/6 of samples to have expressed count >= 10*
sample_cutoff <- (1/6)
count_cutoff <- 10

# Define a function to calculate the fraction of values expressed above the count cutoff
getFE <- function(data,count_cutoff){
  FE <- (sum(data >= count_cutoff)/length(data))
  return(FE)
}
```
Apply the function to all genes, and filter out genes not meeting the sample cutoff

```R
fraction_expressed <- apply(rawdata, 1, getFE, count_cutoff)

keep <- which(fraction_expressed >= sample_cutoff)

rawdata <- rawdata[keep,]

dim(rawdata) # Check dimensions again to see effect of filtering
```
### Specifying the experimental design
DESeq2 relies on a clear understanding of the experimental design, specifically how samples are grouped into different conditions for testing. For this dataset, the experimental design is straightforward: we have six samples representing one condition. This simplicity allows us to construct the experimental design directly within R. However, it's crucial to remember that DESeq2 doesn't verify sample names; it assumes that the column names in our counts matrix precisely match the specified row names in the experimental design.

```R
# construct a mapping of the meta data for our experiment (comparing UHR cell lines to HBR brain tissues)
# in simple terms this is defining the biological condition/label for each experimental replicate
# create a simple one column dataframe to start
metaData <- data.frame('Condition'=c('UHR', 'UHR', 'UHR', 'HBR', 'HBR', 'HBR'))

# convert the "Condition" column to a factor data type, this will determine the direction of log2 fold-changes for the genes (i.e. up or down regulated)
metaData$Condition <- factor(metaData$Condition, levels=c('HBR', 'UHR'))

# set the row names of the metaData dataframe to be the names of our sample replicates from the read counts matrix
rownames(metaData) <- colnames(rawdata)

# view the metadata dataframe
head(metaData)

# check that names of the rawdat columns match the names of the meta data rows
# use the "all" function which tests whether an entire logical vector is TRUE
all(rownames(metaData) == colnames(rawdata))
```
### Construct the DESeq2 object piecing all the data together
Now that the data is correctly formatted, we can consolidate all necessary information for running differential expression into a single object. This object serves as a container for the input data, intermediate computations, and specifies the condition under examination.

```R
# make deseq2 data sets
# here we are setting up our experiment by supplying: (1) the gene counts matrix, (2) the sample/replicate for each column, and (3) the biological conditions we wish to compare.
# this is a simple example that works for many experiments but these can also get more complex
# for example, including designs with multiple variables such as "~ group + condition",
# and designs with interactions such as "~ genotype + treatment + genotype:treatment".
dds <- DESeqDataSetFromMatrix(countData = rawdata, colData = metaData, design = ~Condition)
```
### Running DESeq2
Now that everything is set up, running DESeq2 initiates several crucial steps:

-1. Estimating size factors  
-2. Estimating dispersion  
-3. Implementing "independent filtering" to minimize the number of statistical tests conducted (refer to ?results).    
-4. Fitting Negative Binomial GLM and executing the Wald statistical test.    
-5. Adjusting p-values for multiple testing using the Benjamini and Hochberg method.    

```R
# run the DESeq2 analysis on the "dds" object
dds <- DESeq(dds)

# view the DE results
res <- results(dds)
summary(res)
```

resultsNames(dds)

resLFC <- lfcShrink(dds, coef="Condition_UHR_vs_HBR", type="apeglm")

deGeneResult <- resLFC

summary(res)
summary(deGeneResult)

head(res)
head(deGeneResult)

mapping <- fread("ENSG_ID2Name.txt", header=F)

setnames(mapping, c('ensemblID', 'Symbol'))

head(mapping)

deGeneResult$ensemblID <- rownames(deGeneResult)
deGeneResult <- as.data.table(deGeneResult)
deGeneResult <- merge(deGeneResult, mapping, by='ensemblID', all.x=T)


original_counts = as.data.frame(rawdata)
original_counts[,"ensemblID"] = rownames(rawdata)
deGeneResult = merge(deGeneResult, original_counts, by='ensemblID', all.x=T)
head(deGeneResult)


deGeneResult[order(deGeneResult$padj),]

deGeneResult[order(deGeneResult$log2FoldChange),]

dim(deGeneResult)[1] # number of genes tested
dim(deGeneResult[deGeneResult$padj < 0.05])[1] #number of significant genes

deGeneResultSorted = deGeneResult[order(deGeneResult$padj),]

deGeneResultSignificant = deGeneResultSorted[deGeneResultSorted$padj < 0.05]


setwd(outdir)

# save the final DE result (all genes)  to an output file
fwrite(deGeneResultSorted, file='DE_all_genes_DESeq2.tsv', sep="\t")

normalizedData <- as.data.table(counts(dds, normalized = TRUE),keep.rownames="Gene_ID")
fwrite(normalizedData, file = 'gene_tpm_all_samples.tsv', sep = "\t")

# save the final DE result (significant genes only)  to an output file
fwrite(deGeneResultSignificant, file='DE_sig_genes_DESeq2.tsv', sep="\t")

# save the DESeq2 objects for the data visualization section
saveRDS(dds, 'dds.rds')
saveRDS(res, 'res.rds')
saveRDS(resLFC, 'resLFC.rds')




