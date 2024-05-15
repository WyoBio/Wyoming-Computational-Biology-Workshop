## Setup

#### Load R libraries we will use in this section
```
library(DESeq2)
library(data.table)
library(apeglm)
```
#### Define working dir paths
```
datadir = "/Users/ramshukla/NPRG Dropbox/Computaional Biology Workshop/Day2/Differential Expression"
outdir = "/Users/ramshukla/NPRG Dropbox/Computaional Biology Workshop/Day2/Differential Expression/DE_Results"
```
#### Set working directory to datadir
```
setwd(datadir)
```
#### read in the RNAseq read counts for each gene (produced by featurecounts)
```
rawdata=read.table("featurecounts.txt", header=TRUE, stringsAsFactors=FALSE, row.names=1)
colnames(rawdata)
rawdata <- rawdata[,-c(1,2,3,4,5,6,10)] # Remove columns which are not required
```
#### Extract and edit column names
By default, featurecount outputs the sample name as the file name along with the file path. We will modify this to ensure the correct sample names are used.
```
column_names <- colnames(rawdata)
extracted_names <- gsub(".*\\.([A-Z]+)_Rep([0-9]+)\\.bam", "\\1_Rep\\2", column_names) # Extract "HBR_RepX" and "UHR_RepX" from column names
colnames(rawdata) <- extracted_names
dim(rawdata) # Check dimensions
```


# Require at least 1/6 of samples to have expressed count >= 10
sample_cutoff <- (1/6)
count_cutoff <- 10

#Define a function to calculate the fraction of values expressed above the count cutoff
getFE <- function(data,count_cutoff){
  FE <- (sum(data >= count_cutoff)/length(data))
  return(FE)
}

#Apply the function to all genes, and filter out genes not meeting the sample cutoff
fraction_expressed <- apply(rawdata, 1, getFE, count_cutoff)
keep <- which(fraction_expressed >= sample_cutoff)
rawdata <- rawdata[keep,]

# Check dimensions again to see effect of filtering
dim(rawdata)

metaData <- data.frame('Condition'=c('UHR', 'UHR', 'UHR', 'HBR', 'HBR', 'HBR'))

metaData$Condition <- factor(metaData$Condition, levels=c('HBR', 'UHR'))

rownames(metaData) <- colnames(rawdata)


head(metaData)

all(rownames(metaData) == colnames(rawdata))

dds <- DESeqDataSetFromMatrix(countData = rawdata, colData = metaData, design = ~Condition)

dds <- DESeq(dds)
res <- results(dds)
summary(res)

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




