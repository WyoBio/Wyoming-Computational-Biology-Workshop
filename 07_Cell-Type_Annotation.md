# Cell Type Annotation

### What's the importance of labeling cell types in single-cell RNAseq data?
Assigning cell type annotations in single-cell gene expression data is a dynamic area of research. This task is challenging due to two main factors:

- Gene expression levels often exist on a continuum rather than in discrete categories.
- Variations in gene expression don't always equate to functional differences between cells (Pasquini et al.).

The field is rapidly evolving with the introduction of various annotation tools and the development of multiple scRNA-seq databases (Wang et al.).

One effective approach for cell type annotation is utilizing web resources. These tools are user-friendly and do not require advanced scripting or programming skills. Notable examples include Azimuth, Tabula Sapiens, and SciBet. Researchers often employ multiple tools for cell type annotation to ensure accuracy. Additionally, it is advisable to validate annotations through experiments, statistical analyses, or consultations with subject matter experts.

Methods for assigning cell types typically fall into one of two categories:

- **Databases with Marker Genes for Manual Annotation**: These databases contain marker genes for various cell types (primarily human and mouse). By using these marker genes, you can manually annotate the cells in your datasets. For example, in the 10x Genomics Loupe Browser, you can identify top differentially expressed genes for each cluster. You can then search these genes in the database to determine if they are markers for specific cell types. To ensure accuracy, you may need to search multiple top genes in each cluster. This process can be laborious and time-consuming, depending on the dataset's complexity and your prior knowledge of the cell types.

- **Resources for Automated Cell Type Annotation (Reference-Based)**: These community-developed tools automatically annotate cells by comparing new data with existing references. This article highlights several web tools that do not require programming skills. A major limitation of these tools is that the quality of the results heavily depends on the quality of the pre-annotated reference datasets.

### Setup

```R
library(SingleR)
library(celldex)
library(Seurat)
library(cowplot)

datadir = "/project/biocompworkshop/rshukla/scRNASeq_Results"
setwd(datadir)
list.files()

merged <- readRDS("outdir_single_cell_rna/rep135_clustered.rds")
```

### Choosing the appropriate reference dataset
The function below presents normalized expression values from 830 microarray samples created by ImmGen using pure populations of murine immune cells. The data processing and normalization were conducted according to the methods outlined in Aran, Looney, and Liu et al. (2019). Specifically, CEL files from the Gene Expression Omnibus (GEO; GSE15907 and GSE37448) were acquired, processed, and normalized using the robust multi-array average (RMA) technique on probe-level data. This dataset encompasses 20 broad cell types ("label.main") and 253 finely resolved cell subtypes ("label.fine"). Furthermore, the subtypes have been linked to the Cell Ontology ("label.ont," unless specified as "none" in cell.ont), facilitating additional programmatic inquiries.

```R
#cell typing with single R
#load singleR immgen reference
ref_immgen <- celldex::ImmGenData()
```
Invoking the ImmGenData() function returns a SummarizedExperiment object that includes a matrix of log-expression values along with sample-level labels.

```R
ref_immgen
```
```
class: SummarizedExperiment 
dim: 22134 830 
metadata(0):
assays(1): logcounts
rownames(22134): Zglp1 Vmn2r65 ... Tiparp Kdm1a
rowData names(0):
colnames(830): GSM1136119_EA07068_260297_MOGENE-1_0-ST-V1_MF.11C-11B+.LU_1.CEL
  GSM1136120_EA07068_260298_MOGENE-1_0-ST-V1_MF.11C-11B+.LU_2.CEL ...
  GSM920654_EA07068_201214_MOGENE-1_0-ST-V1_TGD.VG4+24ALO.E17.TH_1.CEL
  GSM920655_EA07068_201215_MOGENE-1_0-ST-V1_TGD.VG4+24ALO.E17.TH_2.CEL
colData names(3): label.main label.fine label.ont
```
Let's now examine the appearance of each label:

```R
head(ref_immgen$label.main, n=10)
```
Analyzing the primary labels reveals general cell types like Macrophages and Monocytes.

```R
[1] "Macrophages" "Macrophages" "Macrophages"
[4] "Macrophages" "Macrophages" "Macrophages"
[7] "Monocytes"   "Monocytes"   "Monocytes"  
[10] "Monocytes"
```
```
head(ref_immgen$label.fine, n=10)
```
Looking at the detailed labels, we observe a refinement of the broader cell types mentioned earlier. Instead of just 6 labels for Macrophages, we now encounter specific subtypes like Macrophages (MF.11C-11B+).

```R
[1] "Macrophages (MF.11C-11B+)" "Macrophages (MF.11C-11B+)"
[3] "Macrophages (MF.11C-11B+)" "Macrophages (MF.ALV)"     
[5] "Macrophages (MF.ALV)"      "Macrophages (MF.ALV)"     
[7] "Monocytes (MO.6+I-)"       "Monocytes (MO.6+2+)"      
[9] "Monocytes (MO.6+2+)"       "Monocytes (MO.6+2+)"
```
```R
head(ref_immgen$label.ont, n=10)
```
By examining the ont labels, we observe that the subtypes are now linked to Cell Ontology IDs.

```R
[1] "CL:0000235" "CL:0000235" "CL:0000235" "CL:0000583" "CL:0000583"
[6] "CL:0000583" "CL:0000576" "CL:0000576" "CL:0000576" "CL:0000576"
```
### Utilizing the ImmGen cell reference with our dataset

```R
#generate predictions for our seurat object
predictions_main = SingleR(test = GetAssayData(merged), 
                      ref = ref_immgen,
                      labels = ref_immgen$label.main)

predictions_fine = SingleR(test = GetAssayData(merged), 
                           ref = ref_immgen,
                           labels = ref_immgen$label.fine)
```














