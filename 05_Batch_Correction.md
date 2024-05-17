# Batch Correction

### Introduction to batch correction
- Excerpt from: *ComBat-seq: batch effect adjustment for RNA-seq count data by Yuqing Zhang, Giovanni Parmigiani, and W Evan Johnson*   

***Genomic data are often produced in batches due to logistical or practical restrictions, but technical variation and differences across batches, often called batch effects, can cause significant heterogeneity across batches of data. Batch effects often result in discrepancies in the statistical distributions across data from different technical processing batches, and can have unfavorable impact on downstream biological analysis … Batch effects often cannot be fully addressed by normalization methods and procedures. The differences in the overall expression distribution of each sample across batch may be corrected by normalization methods, such as transforming the raw counts to (logarithms of) CPM, TPM or RPKM/FPKM, the trimmed mean of M values (TMM), or relative log expression (RLE). However, batch effects in composition, i.e. the level of expression of genes scaled by the total expression (coverage) in each sample, cannot be fully corrected with normalization. … While the overall distribution of samples may be normalized to the same level across batches, individual genes may still be affected by batch-level bias.***


### Setup
In this section, we will demonstrate the principles and application of batch correction using the ComBat-Seq tool in R (Bioconductor). Although we do not expect batch effects in our test data due to its generation at a single center, at one time, and with consistent methodology, we will use a different dataset that is highly related to showcase the impact of batch correction in this module.

The raw count data we will use is from a publicly available RNA-seq dataset obtained from a comprehensive multi-platform comparison of sequencing platforms. This study also explored the impact of generating data across multiple sites, using polyA vs. ribo-reduction for enrichment, and examining the effects of RNA degradation (PMID: 25150835): "Multi-platform and cross-methodological reproducibility of transcriptome profiling by RNA-seq in the ABRF Next-Generation Sequencing Study."

The UHR (cancer cell lines) and HBR (brain tissue) samples used in this course were also utilized in this publication. To assess a significant batch effect, we will perform a differential expression analysis comparing Ribo-depleted ("Ribo") and polyA-enriched ("Poly") samples in the UHR vs. HBR comparison.

You can find the data at this location:




datadir = "~/NPRG Dropbox/Computaional Biology Workshop/Day2/BatchCorrection"
outdir = "~/NPRG Dropbox/Computaional Biology Workshop/Day2/BatchCorrection/Results"
library("sva") 
library("ggplot2")
library("gridExtra")
library("edgeR")
library("UpSetR")

#load in the uncorrected data as raw counts
setwd(datadir)
uncorrected_data = read.table("GSE48035_ILMN.Counts.SampleSubset.ProteinCodingGenes.tsv", header=TRUE, sep="\t", as.is=c(1,2))
setwd(outdir)

#simplify the names of the data columns
# (A = Universal Human Reference RNA and B = Human Brain Reference RNA)
# RNA = polyA enrichment and RIBO = ribosomal RNA depletion
# 1, 2, 3, 4 are replicates
names(uncorrected_data) = c("Gene", "Chr", "UHR_Ribo_1", "UHR_Ribo_2", "UHR_Ribo_3", "UHR_Ribo_4", "HBR_Ribo_1", "HBR_Ribo_2", "HBR_Ribo_3", "HBR_Ribo_4", 
                            "UHR_Poly_1", "UHR_Poly_2", "UHR_Poly_3", "UHR_Poly_4", "HBR_Poly_1", "HBR_Poly_2", "HBR_Poly_3", "HBR_Poly_4")
sample_names = names(uncorrected_data)[3:length(names(uncorrected_data))]

#review data structure
head(uncorrected_data)
dim(uncorrected_data)

#define conditions, library methods, and replicates
conditions = c("UHR", "UHR", "UHR", "UHR", "HBR", "HBR", "HBR", "HBR", "UHR", "UHR", "UHR", "UHR", "HBR", "HBR", "HBR", "HBR")
library_methods = c("Ribo", "Ribo", "Ribo", "Ribo", "Ribo", "Ribo", "Ribo", "Ribo", "Poly", "Poly", "Poly", "Poly", "Poly", "Poly", "Poly", "Poly")
replicates = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4)

#calculate principal components for the uncorrected data
pca_uncorrected_obj = prcomp(uncorrected_data[,sample_names])

#pull PCA values out of the PCA object
pca_uncorrected = as.data.frame(pca_uncorrected_obj[2]$rotation)

#assign labels to the data frame
pca_uncorrected[,"condition"] = conditions
pca_uncorrected[,"library_method"] = library_methods
pca_uncorrected[,"replicate"] = replicates

#plot the PCA
#create a classic 2-dimension PCA plot (first two principal components) with conditions and library methods indicated
cols <- c("UHR" = "#481567FF", "HBR" = "#1F968BFF")
p1 = ggplot(data=pca_uncorrected, aes(x=PC1, y=PC2, color=condition, shape=library_method))
p1 = p1 + geom_point(size=3)
p1 = p1 + stat_ellipse(type="norm", linetype=2)
p1 = p1 + labs(title="PCA, RNA-seq counts for 16 HBR/UHR and Ribo/PolyA samples (uncorrected data)", color="Condition", shape="Library Method")
p1 = p1 + scale_colour_manual(values = cols)
p1



#perform the batch correction

#first we need to transform the format of our groups and batches from names (e.g. "UHR", "HBR", etc.) to numbers (e.g. 1, 2, etc.)
#in the command below "sapply" is used to apply the "switch" command to each element and convert names to numbers as we define
groups = sapply(as.character(conditions), switch, "UHR" = 1, "HBR" = 2, USE.NAMES = F)
batches = sapply(as.character(library_methods), switch, "Ribo" = 1, "Poly" = 2, USE.NAMES = F)

#now run ComBat_seq
corrected_data = ComBat_seq(counts = as.matrix(uncorrected_data[,sample_names]), batch = batches, group = groups)

#join the gene and chromosome names onto the now corrected counts from ComBat_seq
corrected_data = cbind(uncorrected_data[,c("Gene","Chr")], corrected_data)

#compare dimensions of corrected and uncorrected data sets
dim(uncorrected_data)
dim(corrected_data)

#visually compare values of corrected and uncorrected data sets
head(uncorrected_data)
head(corrected_data)


pca_corrected_obj = prcomp(corrected_data[,sample_names])

#pull PCA values out of the PCA object
pca_corrected = as.data.frame(pca_corrected_obj[2]$rotation)

#assign labels to the data frame
pca_corrected[,"condition"] = conditions
pca_corrected[,"library_method"] = library_methods
pca_corrected[,"replicate"] = replicates

#as above, create a PCA plot for comparison to the uncorrected data
cols <- c("UHR" = "#481567FF", "HBR" = "#1F968BFF")
p2 = ggplot(data=pca_corrected, aes(x=PC1, y=PC2, color=condition, shape=library_method))
p2 = p2 + geom_point(size=3)
p2 = p2 + stat_ellipse(type="norm", linetype=2)
p2 = p2 + labs(title="PCA, RNA-seq counts for 16 HBR/UHR and Ribo/PolyA samples (batch corrected data)", color="Condition", shape="Library Method")
p2 = p2 + scale_colour_manual(values = cols)
p2

grid.arrange(p1, p2, nrow = 2)


#perform differential expression analysis on the uncorrected data and batch corrected data sets

#first define the sets of samples to be compared to each other
uhr_ribo_samples = c("UHR_Ribo_1", "UHR_Ribo_2", "UHR_Ribo_3", "UHR_Ribo_4")
uhr_poly_samples = c("UHR_Poly_1", "UHR_Poly_2", "UHR_Poly_3", "UHR_Poly_4")
hbr_ribo_samples = c("HBR_Ribo_1", "HBR_Ribo_2", "HBR_Ribo_3", "HBR_Ribo_4")
hbr_poly_samples = c("HBR_Poly_1", "HBR_Poly_2", "HBR_Poly_3", "HBR_Poly_4")
uhr_samples = c(uhr_ribo_samples, uhr_poly_samples)
hbr_samples = c(hbr_ribo_samples, hbr_poly_samples)

#create a function that will run edgeR (DE analysis) for a particular pair of sample sets
run_edgeR = function(data, group_a_name, group_a_samples, group_b_samples, group_b_name){
  #create a list of all samples for this current comparison
  samples_for_comparison = c(group_a_samples, group_b_samples)
  
  #define the class factor for this pair of sample sets
  class = factor(c(rep(group_a_name,length(group_a_samples)), rep(group_b_name,length(group_b_samples))))
  
  #create a simplified data matrix for only these samples
  rawdata = data[,samples_for_comparison]
  
  #store gene names for later
  genes = rownames(data)
  gene_names = data[,"Gene"]
  
  #make DGElist object
  y = DGEList(counts=rawdata, genes=genes, group=class)
  
  #perform TMM normalization
  y <- calcNormFactors(y)
  
  #estimate dispersion
  y <- estimateCommonDisp(y, verbose=TRUE)
  y <- estimateTagwiseDisp(y)
  
  #perform the differential expression test
  et <- exactTest(y)
  
  #print number of up/down significant genes at FDR = 0.05 significance level and store the DE status in a new variable (de)
  summary(de <- decideTests(et, adjust.method="fdr", p=.05))
  
  #create a matrix of the DE results
  mat <- cbind(
    genes, 
    gene_names,
    sprintf('%0.3f', log10(et$table$PValue)),
    sprintf('%0.3f', et$table$logFC)
  )
  
  #create a version of this matrix that is limited to only the *significant* results
  mat <- mat[as.logical(de),]
  
  #add name to the columns of the final matrix
  colnames(mat) <- c("Gene", "Gene_Name", "Log10_Pvalue", "Log_fold_change")
  
  return(mat)
}

uhr_ribo_vs_hbr_ribo_uncorrected = run_edgeR(data=uncorrected_data, group_a_name="UHR", group_a_samples=uhr_ribo_samples, group_b_name="HBR", group_b_samples=hbr_ribo_samples)
uhr_poly_vs_hbr_poly_uncorrected = run_edgeR(data=uncorrected_data, group_a_name="UHR", group_a_samples=uhr_poly_samples, group_b_name="HBR", group_b_samples=hbr_poly_samples)
uhr_ribo_vs_hbr_poly_uncorrected = run_edgeR(data=uncorrected_data, group_a_name="UHR", group_a_samples=uhr_ribo_samples, group_b_name="HBR", group_b_samples=hbr_poly_samples)
uhr_poly_vs_hbr_ribo_uncorrected = run_edgeR(data=uncorrected_data, group_a_name="UHR", group_a_samples=uhr_poly_samples, group_b_name="HBR", group_b_samples=hbr_ribo_samples)
uhr_vs_hbr_uncorrected = run_edgeR(data=uncorrected_data, group_a_name="UHR", group_a_samples=uhr_samples, group_b_name="HBR", group_b_samples=hbr_samples)

#run the same five comparisons through edgeR using the *batch corrected data*
uhr_ribo_vs_hbr_ribo_corrected = run_edgeR(data=corrected_data, group_a_name="UHR", group_a_samples=uhr_ribo_samples, group_b_name="HBR", group_b_samples=hbr_ribo_samples)
uhr_poly_vs_hbr_poly_corrected = run_edgeR(data=corrected_data, group_a_name="UHR", group_a_samples=uhr_poly_samples, group_b_name="HBR", group_b_samples=hbr_poly_samples)
uhr_ribo_vs_hbr_poly_corrected = run_edgeR(data=corrected_data, group_a_name="UHR", group_a_samples=uhr_ribo_samples, group_b_name="HBR", group_b_samples=hbr_poly_samples)
uhr_poly_vs_hbr_ribo_corrected = run_edgeR(data=corrected_data, group_a_name="UHR", group_a_samples=uhr_poly_samples, group_b_name="HBR", group_b_samples=hbr_ribo_samples)
uhr_vs_hbr_corrected = run_edgeR(data=corrected_data, group_a_name="UHR", group_a_samples=uhr_samples, group_b_name="HBR", group_b_samples=hbr_samples)

#how much of a difference does batch correction make when doing the comparison of all UHR vs all HBR samples?
dim(uhr_vs_hbr_uncorrected)
dim(uhr_vs_hbr_corrected)


listInput1 = list("4 UHR Ribo vs 4 HBR Ribo" = uhr_ribo_vs_hbr_ribo_uncorrected[,"Gene"], 
                  "4 UHR Poly vs 4HBR Poly" = uhr_poly_vs_hbr_poly_uncorrected[,"Gene"],
                  "4 UHR Ribo vs 4 HBR Poly" = uhr_ribo_vs_hbr_poly_uncorrected[,"Gene"],
                  "4 UHR Poly vs 4 HBR Ribo" = uhr_poly_vs_hbr_ribo_uncorrected[,"Gene"],
                  "8 UHR vs 8 HBR" = uhr_vs_hbr_uncorrected[,"Gene"])
upset(fromList(listInput1), order.by = "freq", number.angles=45, point.size=3)


listInput2 = list("4 UHR Ribo vs 4 HBR Ribo" = uhr_ribo_vs_hbr_ribo_corrected[,"Gene"], 
                  "4 UHR Poly vs 4 HBR Poly" = uhr_poly_vs_hbr_poly_corrected[,"Gene"],
                  "4 UHR Ribo vs 4 HBR Poly" = uhr_ribo_vs_hbr_poly_corrected[,"Gene"],
                  "4 UHR Poly vs 4 HBR Ribo" = uhr_poly_vs_hbr_ribo_corrected[,"Gene"],
                  "8 UHR vs 8 HBR" = uhr_vs_hbr_corrected[,"Gene"])
upset(fromList(listInput2), order.by = "freq", number.angles=45, point.size=3)

write.table(uhr_vs_hbr_corrected, file="DE_genes_uhr_vs_hbr_corrected.tsv", quote=FALSE, row.names=FALSE, sep="\t")


