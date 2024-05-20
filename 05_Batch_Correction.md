# Batch Correction

### Introduction to batch correction

Excerpt from: *ComBat-seq: batch effect adjustment for RNA-seq count data by Yuqing Zhang, Giovanni Parmigiani, and W Evan Johnson*   

***Genomic data are often produced in batches due to logistical or practical restrictions, but technical variation and differences across batches, often called batch effects, can cause significant heterogeneity across batches of data. Batch effects often result in discrepancies in the statistical distributions across data from different technical processing batches, and can have unfavorable impact on downstream biological analysis … Batch effects often cannot be fully addressed by normalization methods and procedures. The differences in the overall expression distribution of each sample across batch may be corrected by normalization methods, such as transforming the raw counts to (logarithms of) CPM, TPM or RPKM/FPKM, the trimmed mean of M values (TMM), or relative log expression (RLE). However, batch effects in composition, i.e. the level of expression of genes scaled by the total expression (coverage) in each sample, cannot be fully corrected with normalization. … While the overall distribution of samples may be normalized to the same level across batches, individual genes may still be affected by batch-level bias.***

### Data Overview
In this section, we will demonstrate the principles and application of batch correction using the ComBat-Seq tool in R (Bioconductor). Although we do not expect batch effects in our test data due to its generation at a single center, at one time, and with consistent methodology, we will use a different dataset that is highly related to showcase the impact of batch correction in this module.

The raw count data we will use is from a publicly available RNA-seq dataset obtained from a comprehensive multi-platform comparison of sequencing platforms. This study also explored the impact of generating data across multiple sites, using polyA vs. ribo-reduction for enrichment, and examining the effects of RNA degradation (PMID: 25150835): "Multi-platform and cross-methodological reproducibility of transcriptome profiling by RNA-seq in the ABRF Next-Generation Sequencing Study."

The UHR (cancer cell lines) and HBR (brain tissue) samples used in this course were also utilized in this publication. To assess a significant batch effect, we will perform a differential expression analysis comparing Ribo-depleted ("Ribo") and polyA-enriched ("Poly") samples in the UHR vs. HBR comparison.

**Important Note: Considerations for Batch Correction**

Batch correction requires representation of each condition of interest in each batch. For instance, if all HBR samples were processed using Riboreduction and all UHR samples with PolyA enrichment, we would be unable to effectively model the batch effect versus the condition effect.

You can find the data at this location:

### Setup

```R
# Load R libraries 
library("sva") 
library("ggplot2")
library("gridExtra")
library("edgeR")
library("UpSetR")

# Define working directories 
datadir = "/project/biocompworkshop/Data_Vault"
outdir = "/project/biocompworkshop/rshukla/BatchCorrection_Results"

#load in the uncorrected data as raw counts
setwd(datadir)

uncorrected_data = read.table("GSE48035_ILMN.Counts.SampleSubset.ProteinCodingGenes.tsv", header=TRUE, sep="\t", as.is=c(1,2))

setwd(outdir)

# Simplify the names of the data columns
# (A = Universal Human Reference RNA and B = Human Brain Reference RNA)
# RNA = polyA enrichment and RIBO = ribosomal RNA depletion
# 1, 2, 3, 4 are replicates
names(uncorrected_data) = c("Gene", "Chr", "UHR_Ribo_1", "UHR_Ribo_2", "UHR_Ribo_3", "UHR_Ribo_4", "HBR_Ribo_1", "HBR_Ribo_2", "HBR_Ribo_3", "HBR_Ribo_4", 
                            "UHR_Poly_1", "UHR_Poly_2", "UHR_Poly_3", "UHR_Poly_4", "HBR_Poly_1", "HBR_Poly_2", "HBR_Poly_3", "HBR_Poly_4")

sample_names = names(uncorrected_data)[3:length(names(uncorrected_data))]

# Review data structure
head(uncorrected_data)
dim(uncorrected_data)

# Define conditions, library methods, and replicates
conditions = c("UHR", "UHR", "UHR", "UHR", "HBR", "HBR", "HBR", "HBR", "UHR", "UHR", "UHR", "UHR", "HBR", "HBR", "HBR", "HBR")
library_methods = c("Ribo", "Ribo", "Ribo", "Ribo", "Ribo", "Ribo", "Ribo", "Ribo", "Poly", "Poly", "Poly", "Poly", "Poly", "Poly", "Poly", "Poly")
replicates = c(1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4)"
```

### Conduct PCA on Uncorrected Counts 
PCA Analysis and Batch Effect Identification

Principal component analysis (PCA) serves as a valuable tool for spotting potential batch effects within your data. The primary strategy involves using PCA to uncover patterns of similarity or difference in expression signatures among your samples and determining if these patterns align with expected biological conditions. The PCA plot can be annotated with biological conditions as well as potential sources of batch effects, such as sequencing source, data generation date, lab technician, library construction kit batches, matrigel batches, mouse litters, software or instrumentation versions, and so on.

PCA is a method for reducing the dimensionality of large datasets, such as thousands of gene expression values across multiple samples. It aims to represent a vast set of variables with a smaller set that captures the essential information from the larger set. While PCA is a versatile exploratory data analysis tool with numerous applications and complexities, the specifics are beyond the scope of this demonstration focusing on batch effect correction.

In the context of this module, PCA allows for visualizing samples as clusters based on their overall gene expression patterns. These clusters, typically visualized in 2D or interactive 3D plots, can aid in interpreting high-level differences between samples and validating expectations regarding similarity between conditions and replicates.

We will conduct PCA analysis both before and after batch correction, labeling samples according to biological condition (UHR vs HBR) and library preparation type (Ribo vs PolyA).

```R
# Calculate principal components for the uncorrected data
pca_uncorrected_obj = prcomp(uncorrected_data[,sample_names])

# Pull PCA values out of the PCA object
pca_uncorrected = as.data.frame(pca_uncorrected_obj[2]$rotation)

# Assign labels to the data frame
pca_uncorrected[,"condition"] = conditions
pca_uncorrected[,"library_method"] = library_methods
pca_uncorrected[,"replicate"] = replicates

# Plot the PCA
# Create a classic 2-dimension PCA plot (first two principal components) with conditions and library methods indicated
cols <- c("UHR" = "#481567FF", "HBR" = "#1F968BFF")
p1 = ggplot(data=pca_uncorrected, aes(x=PC1, y=PC2, color=condition, shape=library_method))
p1 = p1 + geom_point(size=3)
p1 = p1 + stat_ellipse(type="norm", linetype=2)
p1 = p1 + labs(title="PCA, RNA-seq counts for 16 HBR/UHR and Ribo/PolyA samples (uncorrected data)", color="Condition", shape="Library Method")
p1 = p1 + scale_colour_manual(values = cols)
p1
```
### The ComBat-Seq Package for Batch Correction
The ComBat-Seq package is part of the SVA package, designed for Surrogate Variable Analysis. This collection of methods aims to remove batch effects and unwanted variations in large datasets. ComBat-Seq, a modification of the ComBat method, is specifically tailored for count-based data in bulk RNA-seq datasets.

- Key advantages of ComBat-Seq include:
  - 1. Utilization of a negative binomial regression model, fitting the characteristics of bulk RNA-seq count data.
  - 2. Generation of corrected data retaining the count nature, suitable for various differential expression analysis methods like EdgeR and DESeq2.

ComBat-Seq offers a straightforward interface with default settings for most arguments. Basic documentation for these arguments can be found here: https://github.com/zhangyuqing/ComBat-seq

Each of the ComBat-Seq arguments is briefly explained below:

**counts:** This is your matrix of gene expression read counts (raw counts). Each row represents a gene, each column represents a sample, and each cell contains an integer count indicating the RNA-seq counts observed for that gene/sample combination. In R, ensure that this data is in matrix format (use `as.matrix()` if necessary) before passing it into ComBat-Seq.

**batch:** This is a vector describing the batches you are concerned about. For instance, if you have eight samples where the first four were processed using library kit (A) and the last four using library kit (B), your batch vector would be defined as: `c(1,1,1,1,2,2,2,2)`.

**group = NULL:** This is a vector describing your biological condition of interest. For example, if your experiment compares drug-treated and untreated cells in pairs with 4 biological replicates, your group vector would be: `c(1,2,1,2,1,2,1,2)`.

**covar_mod = NULL:** If you have multiple biological conditions of interest, define them using covar_mod (covariates) instead of group. For instance, consider the same experiment as above but with alternating male and female cells for each replicate. You would define a covariate matrix for covar_mod as follows:

```R
# Treatment_group = c(1,2,1,2,1,2,1,2)
# Sex_group = c(1,1,2,2,1,1,2,2)
# Covariate_matrix = cbind(treatment_group, sex_group)
```

**full_mod = TRUE:** If set to TRUE, include the condition of interest in the model. This is typically recommended as the default setting. We haven't encountered a compelling reason to set this to FALSE.

**shrink = FALSE:** Determines whether to apply shrinkage on parameter estimation.

**shrink.disp = FALSE:** Specifies whether to apply shrinkage on dispersion.

**gene.subset.n = NULL:** Indicates the number of genes to use in empirical Bayes estimation, only applicable when shrink = TRUE.

A detailed discussion of shrinkage (related to the `shrink`, `shrink.disp`, and `gene_subset.n` arguments) is beyond the scope of this tutorial. Briefly, shrinkage refers to a set of methods that attempt to correct for gene-specific variability in the counts observed in RNA-seq datasets. More specifically, it relates to the dispersion parameter of the negative binomial distribution used to model RNA-seq count data that can suffer from overdispersion. The dispersion parameter describes how much variance deviates from the mean. In simple terms, shrinkage methods are an attempt to correct for problematic dispersion. A more detailed discussion of these statistical concepts can be found in the DESeq2 paper. However, for our purposes here, the bottom line is that the ComBat-Seq authors state that “We have shown that applying empirical Bayes shrinkage is not necessary for ComBat-seq because the approach is already sufficiently robust due to the distributional assumption.” So we will leave these arguments at their default FALSE settings.

### Perform the batch correction

```R
# First we need to transform the format of our groups and batches from names (e.g. "UHR", "HBR", etc.) to numbers (e.g. 1, 2, etc.)
# in the command below "sapply" is used to apply the "switch" command to each element and convert names to numbers as we define
groups = sapply(as.character(conditions), switch, "UHR" = 1, "HBR" = 2, USE.NAMES = F)
batches = sapply(as.character(library_methods), switch, "Ribo" = 1, "Poly" = 2, USE.NAMES = F)

# Now run ComBat_seq
corrected_data = ComBat_seq(counts = as.matrix(uncorrected_data[,sample_names]), batch = batches, group = groups)

# Join the gene and chromosome names onto the now corrected counts from ComBat_seq
corrected_data = cbind(uncorrected_data[,c("Gene","Chr")], corrected_data)

# Compare dimensions of corrected and uncorrected data sets
dim(uncorrected_data)
dim(corrected_data)

# Visually compare values of corrected and uncorrected data sets
head(uncorrected_data)
head(corrected_data)

pca_corrected_obj = prcomp(corrected_data[,sample_names])

# Pull PCA values out of the PCA object
pca_corrected = as.data.frame(pca_corrected_obj[2]$rotation)

# Assign labels to the data frame
pca_corrected[,"condition"] = conditions
pca_corrected[,"library_method"] = library_methods
pca_corrected[,"replicate"] = replicates

# As above, create a PCA plot for comparison to the uncorrected data
cols <- c("UHR" = "#481567FF", "HBR" = "#1F968BFF")
p2 = ggplot(data=pca_corrected, aes(x=PC1, y=PC2, color=condition, shape=library_method))
p2 = p2 + geom_point(size=3)
p2 = p2 + stat_ellipse(type="norm", linetype=2)
p2 = p2 + labs(title="PCA, RNA-seq counts for 16 HBR/UHR and Ribo/PolyA samples (batch corrected data)", color="Condition", shape="Library Method")
p2 = p2 + scale_colour_manual(values = cols)
p2

grid.arrange(p1, p2, nrow = 2)
```
How does batch correction influence differential gene expression results? We will utilize UpSet plots to explore the overlap of significant DE genes identified in the following comparisons:

- UHR-Ribo vs HBR-Ribo (same library type, 4 vs 4 replicates)
- UHR-Poly vs HBR-Poly (same library type, 4 vs 4 replicates)
- UHR-Ribo vs HBR-Poly (different library types, 4 vs 4 replicates)
- UHR-Poly vs HBR-Ribo (different library types, 4 vs 4 replicates)
- UHR-Comb vs HBR-Comb (combined library types, 8 vs 8 replicates)

These five sets of differential expression analyses will be conducted using both the uncorrected and corrected data sets.

```R
# Perform differential expression analysis on the uncorrected data and batch corrected data sets

# First define the sets of samples to be compared to each other
uhr_ribo_samples = c("UHR_Ribo_1", "UHR_Ribo_2", "UHR_Ribo_3", "UHR_Ribo_4")
uhr_poly_samples = c("UHR_Poly_1", "UHR_Poly_2", "UHR_Poly_3", "UHR_Poly_4")
hbr_ribo_samples = c("HBR_Ribo_1", "HBR_Ribo_2", "HBR_Ribo_3", "HBR_Ribo_4")
hbr_poly_samples = c("HBR_Poly_1", "HBR_Poly_2", "HBR_Poly_3", "HBR_Poly_4")
uhr_samples = c(uhr_ribo_samples, uhr_poly_samples)
hbr_samples = c(hbr_ribo_samples, hbr_poly_samples)

# Create a function that will run edgeR (DE analysis) for a particular pair of sample sets
run_edgeR = function(data, group_a_name, group_a_samples, group_b_samples, group_b_name){
  # Create a list of all samples for this current comparison
  samples_for_comparison = c(group_a_samples, group_b_samples)
  
  # Define the class factor for this pair of sample sets
  class = factor(c(rep(group_a_name,length(group_a_samples)), rep(group_b_name,length(group_b_samples))))
  
  # Create a simplified data matrix for only these samples
  rawdata = data[,samples_for_comparison]
  
  # Store gene names for later
  genes = rownames(data)
  gene_names = data[,"Gene"]
  
  # Make DGElist object
  y = DGEList(counts=rawdata, genes=genes, group=class)
  
  # Perform TMM normalization
  y <- calcNormFactors(y)
  
  # Estimate dispersion
  y <- estimateCommonDisp(y, verbose=TRUE)
  y <- estimateTagwiseDisp(y)
  
  # Perform the differential expression test
  et <- exactTest(y)
  
  # Print number of up/down significant genes at FDR = 0.05 significance level and store the DE status in a new variable (de)
  summary(de <- decideTests(et, adjust.method="fdr", p=.05))
  
  # Create a matrix of the DE results
  mat <- cbind(
    genes, 
    gene_names,
    sprintf('%0.3f', log10(et$table$PValue)),
    sprintf('%0.3f', et$table$logFC)
  )
  
  # Create a version of this matrix that is limited to only the *significant* results
  mat <- mat[as.logical(de),]
  
  # Add name to the columns of the final matrix
  colnames(mat) <- c("Gene", "Gene_Name", "Log10_Pvalue", "Log_fold_change")
  
  return(mat)
}

uhr_ribo_vs_hbr_ribo_uncorrected = run_edgeR(data=uncorrected_data, group_a_name="UHR", group_a_samples=uhr_ribo_samples, group_b_name="HBR", group_b_samples=hbr_ribo_samples)
uhr_poly_vs_hbr_poly_uncorrected = run_edgeR(data=uncorrected_data, group_a_name="UHR", group_a_samples=uhr_poly_samples, group_b_name="HBR", group_b_samples=hbr_poly_samples)
uhr_ribo_vs_hbr_poly_uncorrected = run_edgeR(data=uncorrected_data, group_a_name="UHR", group_a_samples=uhr_ribo_samples, group_b_name="HBR", group_b_samples=hbr_poly_samples)
uhr_poly_vs_hbr_ribo_uncorrected = run_edgeR(data=uncorrected_data, group_a_name="UHR", group_a_samples=uhr_poly_samples, group_b_name="HBR", group_b_samples=hbr_ribo_samples)
uhr_vs_hbr_uncorrected = run_edgeR(data=uncorrected_data, group_a_name="UHR", group_a_samples=uhr_samples, group_b_name="HBR", group_b_samples=hbr_samples)

# Run the same five comparisons through edgeR using the *batch corrected data*
uhr_ribo_vs_hbr_ribo_corrected = run_edgeR(data=corrected_data, group_a_name="UHR", group_a_samples=uhr_ribo_samples, group_b_name="HBR", group_b_samples=hbr_ribo_samples)
uhr_poly_vs_hbr_poly_corrected = run_edgeR(data=corrected_data, group_a_name="UHR", group_a_samples=uhr_poly_samples, group_b_name="HBR", group_b_samples=hbr_poly_samples)
uhr_ribo_vs_hbr_poly_corrected = run_edgeR(data=corrected_data, group_a_name="UHR", group_a_samples=uhr_ribo_samples, group_b_name="HBR", group_b_samples=hbr_poly_samples)
uhr_poly_vs_hbr_ribo_corrected = run_edgeR(data=corrected_data, group_a_name="UHR", group_a_samples=uhr_poly_samples, group_b_name="HBR", group_b_samples=hbr_ribo_samples)
uhr_vs_hbr_corrected = run_edgeR(data=corrected_data, group_a_name="UHR", group_a_samples=uhr_samples, group_b_name="HBR", group_b_samples=hbr_samples)

# How much of a difference does batch correction make when doing the comparison of all UHR vs all HBR samples?
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
```

An UpSet plot serves as an alternative to a Venn Diagram, displaying the overlap (intersection) among multiple sets of values. In this instance, we are assessing the intersection of genes identified as significantly differentially expressed (DE) across five distinct comparisons. The connected black circles represent each combination of sets under consideration. The bar graph above each column visualizes the number of shared genes within those sets. For instance, the initial column contains all five black circles, and the corresponding bar indicates the count of genes found in all five DE comparisons conducted.



