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
library(ggplot2)
library(viridis)

datadir = "/project/biocompworkshop/rshukla/scRNASeq_Results"
setwd(datadir)
list.files()

merged <- readRDS("rep135_clustered.rds")
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
# Generate predictions for our seurat object
predictions_main = SingleR(test = GetAssayData(merged), 
                      ref = ref_immgen,
                      labels = ref_immgen$label.main)

predictions_fine = SingleR(test = GetAssayData(merged), 
                           ref = ref_immgen,
                           labels = ref_immgen$label.fine)
```
What is the structure of these objects?

```R
head(predictions_main)
```
You should see a format similar to the one below, where each row represents a barcode and the columns correspond to scores, labels (assigned cell type), delta, and pruned.labels.

```
DataFrame with 6 rows and 4 columns
                                                         scores      labels delta.next pruned.labels
                                                       <matrix> <character>  <numeric>   <character>
Rep1_ICBdT_AAACCTGAGCCAACAG-1 0.4156037:0.4067582:0.2845856:...         NKT  0.0124615           NKT
Rep1_ICBdT_AAACCTGAGCCTTGAT-1 0.4551058:0.3195934:0.2282272:...     B cells  0.1355124       B cells
Rep1_ICBdT_AAACCTGAGTACCGGA-1 0.0717647:0.0621878:0.0710026:... Fibroblasts  0.1981683   Fibroblasts
Rep1_ICBdT_AAACCTGCACGGCCAT-1 0.2774994:0.2569566:0.2483387:...    NK cells  0.0577608      NK cells
Rep1_ICBdT_AAACCTGCACGGTAAG-1 0.3486259:0.3135662:0.3145100:...     T cells  0.1038542       T cells
Rep1_ICBdT_AAACCTGCATGCCACG-1 0.0399733:0.0229926:0.0669236:... Fibroblasts  0.2443470   Fibroblasts
```
#### What information do these columns provide us?

- The **scores** column includes a matrix for each barcode, indicating SingleR's confidence in assigning each cell type to the corresponding barcode row.
  - Each cell within this matrix reflects the level of confidence for a specific cell type assignment.
- The **labels** column represents SingleR's most confident assignment for each barcode, highlighting the predominant cell type identification.
- The **delta** column encompasses the "delta" value per cell, indicating the difference between the score for the assigned label and the median score across all labels.
  - A small delta suggests uniform confidence across labels, possibly indicating less meaningful assigned labels.
- SingleR can eliminate cells with low delta values caused by:
  - Ambiguous assignments with closely related reference labels.
  - Incorrect assignments that poorly match all reference labels.
- The **pruned.labels** column contains "cleaner" or more reliable labels after SingleR has discarded cells with low delta values, ensuring more accurate cell type assignments.

How many cells have labels with low confidence?

```R
unique(predictions_main$pruned.labels)
```
```
 [1] "NKT"               "B cells"           "Fibroblasts"       "NK cells"          "T cells"           "Neutrophils"      
 [7] "DC"                "Monocytes"         "ILC"               "Epithelial cells"  "Macrophages"       "Basophils"        
[13] "Tgd"               "Mast cells"        "Endothelial cells" NA                  "Stem cells"        "Stromal cells"    
[19] "B cells, pro"
```
```R
table(predictions_main$pruned.labels)
```
```
          B cells      B cells, pro         Basophils                DC Endothelial cells  Epithelial cells       Fibroblasts 
             3219                 3                33               295                67              1238               577 
              ILC       Macrophages        Mast cells         Monocytes          NK cells               NKT       Neutrophils 
              752               454                11               617               562              2241                92 
       Stem cells     Stromal cells           T cells               Tgd 
                2                18             12631               189
```
```R
table(predictions_main$labels)
```
```
          B cells      B cells, pro         Basophils                DC Endothelial cells  Epithelial cells       Fibroblasts 
             3253                 3                37               295                71              1238               589 
              ILC       Macrophages        Mast cells         Monocytes          NK cells               NKT       Neutrophils 
              763               459                11               633               565              2249                92 
       Stem cells     Stromal cells           T cells               Tgd 
                2                18             12714               193
```
```R
summary(is.na(predictions_main$pruned.labels))
```
```
   Mode   FALSE    TRUE 
logical   23001     184
```
Are there more or fewer pruned labels for the fine labels?

```R
summary(is.na(predictions_fine$pruned.labels))
```
```
   Mode   FALSE    TRUE 
logical   23005     180
```
Having gained an understanding of the structure and content of the SingleR dataframe, let's move on to visualizing the data.

```R
plotDeltaDistribution(predictions_main, ncol = 4, dots.on.top = FALSE)

plotScoreHeatmap(predictions_main)
```
Instead of solely focusing on the singleR dataframe, we have the option to incorporate the labels into our Seurat data object as a metadata field. Let's proceed with adding the cell type labels to our Seurat object.

```R
#add main labels to object
merged[['immgen_singler_main']] = rep('NA', ncol(merged))
merged$immgen_singler_main[rownames(predictions_main)] = predictions_main$labels

#add fine labels to object
merged[['immgen_singler_fine']] = rep('NA', ncol(merged))
merged$immgen_singler_fine[rownames(predictions_fine)] = predictions_fine$labels
```
Now, let's create a visualization to observe the cell typing within our data.

What are the variations in the relative cell composition among our samples?

```R
#visualizing the relative proportion of cell types across our samples
library(viridis)
library(ggplot2)
ggplot(merged[[]], aes(x = orig.ident, fill = immgen_singler_main)) + geom_bar(position = "fill") + scale_fill_viridis(discrete = TRUE)
```
We can alternatively rotate the samples and cell labels to view the data as shown below.

```R
ggplot(merged[[]], aes(x = immgen_singler_main, fill = orig.ident)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_viridis(discrete = TRUE)

ggplot(merged[[]], aes(x = immgen_singler_fine, fill = orig.ident)) + geom_bar(position = "fill") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_viridis(discrete = TRUE) 
```
How are our cluster definitions related to the annotations of cell types?

```R
#plotting cell types on our umaps
DimPlot(merged, group.by = c("immgen_singler_main"))

DimPlot(merged, group.by = c("immgen_singler_fine")) + NoLegend()
```
### What variances do we observe in our cell annotations when employing an alternate reference set?
The chosen dataset originates from the Benayoun Lab, which curated, acquired, and processed datasets from GEO matching sorted cell types as described by Benayoun et al., 2019. It comprises 358 mouse RNA-seq samples annotated across 18 primary cell types ("label.main"), further categorized into 28 subtypes ("label.fine"). These subtypes have been aligned with the Cell Ontology, similar to the ImmGen reference.

```R
ref_mouserna <- celldex::MouseRNAseqData()
```
Upon examining this reference, we can observe that it includes the `main`, `fine`, and `ont` labels similar to what we observed in the ImmGen reference.

```R
ref_mouserna
```
```
class: SummarizedExperiment 
dim: 21214 358 
metadata(0):
assays(1): logcounts
rownames(21214): Xkr4 Rp1 ... LOC100039574 LOC100039753
rowData names(0):
colnames(358): ERR525589Aligned ERR525592Aligned ... SRR1044043Aligned SRR1044044Aligned
colData names(3): label.main label.fine label.ont
```
```R
predictions_mouse_main = SingleR(test = GetAssayData(merged), 
                      ref = ref_mouserna,
                      labels = ref_mouserna$label.main)

predictions_mouse_fine = SingleR(test = GetAssayData(merged), 
                           ref = ref_mouserna,
                           labels = ref_mouserna$label.fine)
```
```
merged[['mouserna_singler_main']] = rep('NA', ncol(merged))
merged$mouserna_singler_main[rownames(predictions_mouse_main)] = predictions_mouse_main$labels

#add fine labels to object
merged[['mouserna_singler_fine']] = rep('NA', ncol(merged))
merged$mouserna_singler_fine[rownames(predictions_mouse_fine)] = predictions_mouse_fine$labels
```
Let's visualize the differences in cell labeling between the primary labels from ImmGen and MouseRNA datasets.

```R
table(predictions_main$labels)
table(predictions_mouse_main$labels)
```
```
          B cells      B cells, pro         Basophils                DC Endothelial cells  Epithelial cells 
             3253                 3                37               295                71              1238 
      Fibroblasts               ILC       Macrophages        Mast cells         Monocytes          NK cells 
              589               763               459                11               633               565 
              NKT       Neutrophils        Stem cells     Stromal cells           T cells               Tgd 
             2249                92                 2                18             12714               193 


       Adipocytes           B cells   Dendritic cells Endothelial cells       Fibroblasts      Granulocytes 
                4              3325                94               245               979               116 
      Hepatocytes       Macrophages         Monocytes          NK cells           T cells 
              652               185              1013              1341             15231 
```
```R
p1 <- DimPlot(merged, group.by = c("immgen_singler_main")) + scale_colour_viridis(option = 'turbo', discrete = TRUE)
p2 <- DimPlot(merged, group.by = c("mouserna_singler_main")) + scale_colour_viridis(option = 'turbo', discrete = TRUE)
p <- plot_grid(p1, p2, ncol = 2)
p
```

### Annotation Datasets: A Quick Insight
The choice of reference data significantly influences annotation outcomes, with a broader label spectrum in the reference dataset being crucial for accurate results compared to our test data. Trusting the original authors' labeling of reference samples can be a leap of faith, leading to variations in performance due to differences in sample preparation quality among references. Ideally, we prefer references aligned with our test dataset's technology or protocol, although this is less critical when using SingleR() for well-defined cell types.

For detailed insights and recommendations on reference selection, users are encouraged to consult the relevant vignette. Notably, while the ImmGen dataset and other references were initially part of SingleR, they are now in the celldex package for wider use across Bioconductor. Additionally, users can supply their own reference datasets, needing only log-expression values and cell/sample labels.

### Gaining insights into clustering through cell typing.
We're aware that the UMAP shape shifts when altering the number of PCs, but what drives this change in UMAP shape? Let's experiment by generating a UMAP using only 5 PCs.

```R
merged_5PC <- RunUMAP(merged, dims = 1:5)

DimPlot(merged_5PC, label = TRUE, group.by = 'immgen_singler_main')

DimHeatmap(merged, dims = 1:5, balanced = TRUE, cells = 500)
```
Upon closer observation, we notice that the cells tend to aggregate into broader clusters. For instance, immune cells amalgamate into a single large cluster. Enhancing the number of PCs allows us to leverage more information, enabling finer distinctions among specific cell types. The clear delineation of the B cell cluster in our previous analysis was achievable due to the sufficient genetic expression data included in the selected PCs.

### Saving data

```R
saveRDS(merged, file = "preprocessed_object.rds")
```























