# Visualization
### Setup
We will establish working directories, load the necessary R libraries, and import the DESeqDataSet object and the results table generated in the previous section into the R environment.

```R
# Set working directories
datadir = "/Users/ramshukla/NPRG Dropbox/Computaional Biology Workshop/Day2/Differential Expression"
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

plotMA(res, ylim=c(-2, 2))
plotMA(resLFC, ylim=c(-2,2))


plotCounts(dds, gene='ENSG00000100167', intgroup = 'Condition')
plotCounts(dds, gene='ENSG00000185686', intgroup = 'Condition')


rld <- rlog(dds, blind=F)

rld
head(assay(dds))

head(assay(rld))

t(assay(rld))[1:6,1:5]

dist(t(assay(rld)))

sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix(sampleDists)

head(sampleDistMatrix)

pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists)

sampleCorrs <- cor(assay(rld), method="pearson")
sampleCorrMatrix <- as.matrix(sampleCorrs)
head(sampleCorrMatrix)

pheatmap(sampleCorrMatrix)

pheatmap(mat=t(assay(rld)), show_colnames = FALSE)


setwd(datadir)
dir()

gene_expression=read.table("gene_tpm_all_samples.tsv", header=TRUE, stringsAsFactors=FALSE, row.names=1)

gene_names=read.table("ENSG_ID2Name.txt", header=TRUE, stringsAsFactors=FALSE)
colnames(gene_names)=c("gene_id","gene_name")

setwd(outdir)

results_genes <-read.table("DE_all_genes_DESeq2.tsv", sep="\t", header=T, stringsAsFactors = F)

head(gene_expression)
colnames(gene_expression)

row.names(gene_expression)

dim(gene_expression)

gene_expression[1:3,c(1:3,6)]

gene_expression[1:3, c("HBR_Rep1","HBR_Rep2","HBR_Rep3","UHR_Rep3")]

head(results_genes)
dim(results_genes)

colours()
data_colors=c("tomato1","tomato2","tomato3","royalblue1","royalblue2","royalblue3")

min_nonzero=1

# Set the columns for finding TPM and create shorter names for figures
data_columns=c(1:6)
short_names=c("HBR_1","HBR_2","HBR_3","UHR_1","UHR_2","UHR_3")

boxplot(log2(gene_expression[,data_columns]+min_nonzero), col=data_colors, names=short_names, las=2, ylab="log2(TPM)", main="Distribution of TPMs for all 6 libraries")

x = gene_expression[,"UHR_Rep1"]
y = gene_expression[,"UHR_Rep2"]


plot(x=log2(x+min_nonzero), y=log2(y+min_nonzero), pch=16, col="blue", cex=0.25, xlab="TPM (UHR, Replicate 1)", ylab="TPM (UHR, Replicate 2)", main="Comparison of expression values for a pair of replicates")

abline(a=0,b=1)

rs=cor(x,y)^2
legend("topleft", paste("R squared = ", round(rs, digits=3), sep=""), lwd=1, col="black")


colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
smoothScatter(x=log2(x+min_nonzero), y=log2(y+min_nonzero), xlab="TPM (UHR, Replicate 1)", ylab="TPM (UHR, Replicate 2)", main="Comparison of expression values for a pair of replicates", colramp=colors, nbin=200)


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

par(mfrow=c(1,3))
plotCor("UHR_Rep1", "UHR_Rep2", "UHR_1 vs UHR_2")
plotCor("UHR_Rep2", "UHR_Rep3", "UHR_2 vs UHR_3")
plotCor("UHR_Rep1", "UHR_Rep3", "UHR_1 vs UHR_3")

gene_expression[,"sum"]=apply(gene_expression[,data_columns], 1, sum)

i = which(gene_expression[,"sum"] > 5)
r=cor(gene_expression[i,data_columns], use="pairwise.complete.obs", method="pearson")
r


d=1-r
mds=cmdscale(d, k=2, eig=TRUE)
par(mfrow=c(1,1))
plot(mds$points, type="n", xlab="", ylab="", main="MDS distance plot (all non-zero genes)", xlim=c(-0.12,0.12), ylim=c(-0.12,0.12))
points(mds$points[,1], mds$points[,2], col="grey", cex=2, pch=16)
text(mds$points[,1], mds$points[,2], short_names, col=data_colors)

sig=which(results_genes$pvalue<0.05)
hist(results_genes[sig,"log2FoldChange"], breaks=50, col="seagreen", xlab="log2(Fold change) UHR vs HBR", main="Distribution of differential expression values")
abline(v=-2, col="black", lwd=2, lty=2)
abline(v=2, col="black", lwd=2, lty=2)
legend("topleft", "Fold-change > 4", lwd=2, lty=2)


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
topn = order(results_genes[sig,"padj"])[1:25]
text(x[topn], y[topn], results_genes[topn,"Symbol"], col="black", cex=0.75, srt=45)


mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

#Create a subset of significant genes with p-value<0.05 and log2 fold-change >= 2
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


