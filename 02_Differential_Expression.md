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

- Estimating size factors  
- Estimating dispersion  
- Implementing "independent filtering" to minimize the number of statistical tests conducted (refer to ?results).    
- Fitting Negative Binomial GLM and executing the Wald statistical test.    
- Adjusting p-values for multiple testing using the Benjamini and Hochberg method.    

```R
# run the DESeq2 analysis on the "dds" object
dds <- DESeq(dds)

# view the DE results
res <- results(dds)
summary(res)
```
### Log-fold change shrinkage
Shrinking the log-fold change values is a common practice aimed at minimizing exaggerated differences that arise from genes with low counts, leading to significant variability that may not reflect true biological signals. Take, for instance, a gene observed in two samples: one sample with 1 read and another with 6 reads. This scenario would yield a 6-fold change, which might not accurately represent the underlying biology. Several algorithms exist for shrinking log2 fold change values, and in this case, we'll utilize the apeglm algorithm. Please note that the apeglm package must be installed to implement this approach.

```R
# shrink the log2 fold change estimates
# shrinkage of effect size (log fold change estimates) is useful for visualization and ranking of genes
# In simplistic terms, the goal of calculating "dispersion estimates" and "shrinkage" is also to account for the problem that
# genes with low mean counts across replicates tend of have higher variability than those with higher mean counts.
# Shrinkage attempts to correct for this. For a more detailed discussion of shrinkage refer to the DESeq2 vignette

# first get the name of the coefficient (log fold change) to shrink
resultsNames(dds)

# now apply the shrinkage approach
resLFC <- lfcShrink(dds, coef="Condition_UHR_vs_HBR", type="apeglm")

# make a copy of the shrinkage results to manipulate
deGeneResult <- resLFC

#contrast the values for a few genes before and after shrinkage
head(res)
head(deGeneResult)
```
### Annotate gene symbols onto the DE results
The initial run of DESeq2 utilized Ensembl gene IDs as identifiers, which can be less intuitive for human interpretation. To enhance the interpretability of results, gene symbols were incorporated into the list of differentially expressed genes, making the outcomes more user-friendly and easier to understand.

##### - ENSG_ID2Name.txt is located in here: 

```R
# read in gene ID to name mappings (using "fread" an alternative to "read.table")
mapping <- fread("ENSG_ID2Name.txt", header=F)

# add names to the columns in the "mapping" dataframe
setnames(mapping, c('ensemblID', 'Symbol'))

# view the first few lines of the gene ID to name mappings
head(mapping)

# merge on gene names
deGeneResult$ensemblID <- rownames(deGeneResult)
deGeneResult <- as.data.table(deGeneResult)
deGeneResult <- merge(deGeneResult, mapping, by='ensemblID', all.x=T)

# merge the original raw count values onto this final dataframe to aid interpretation
original_counts = as.data.frame(rawdata)
original_counts[,"ensemblID"] = rownames(rawdata)
deGeneResult = merge(deGeneResult, original_counts, by='ensemblID', all.x=T)

# view the final merged dataframe
# based on the raw counts and fold change values what does a negative fold change mean with respect to our two conditions HBR and UHR?
head(deGeneResult)
```

### Data filitering
After completing the DE analysis, it's valuable to review and refine the data frames, focusing on genes that are relevant. Basic data manipulation techniques are applied to filter for significant genes based on specific thresholds, streamlining the analysis to key insights.

```R
# view the top genes according to adjusted p-value
deGeneResult[order(deGeneResult$padj),]

# view the top genes according to fold change
deGeneResult[order(deGeneResult$log2FoldChange),]

# determine the number of up/down significant genes at FDR < 0.05 significance level
dim(deGeneResult)[1] # number of genes tested
dim(deGeneResult[deGeneResult$padj < 0.05])[1] #number of significant genes

# order the DE results by adjusted p-value
deGeneResultSorted = deGeneResult[order(deGeneResult$padj),]

# create a filtered data frame that limits to only the significant DE genes (adjusted p.value < 0.05)
deGeneResultSignificant = deGeneResultSorted[deGeneResultSorted$padj < 0.05]
```

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




