#### IGV (Optional)
After this lab, you will be able to:

- Visualize a variety of genomic data
- Quickly navigate around the genome
- Visualize read alignments
- Validate SNP/SNV calls and structural re-arrangements by eye

We will use publicly available Illumina sequence data derived from a 52-year-old Caucasian woman with breast cancer. Specifically, we will utilize sequence read alignments filtered to the region Chromosome 21: 19,000,000-20,000,000. The data is available in the "IGV" folder at the following path: `/project/biocompworkshop/rshukla/IGV`.

File names: 
- HCC1143.normal.21.19M-20M.bam
- HCC1143.normal.21.19M-20M.bam.bai

To visualize read alignments, we will use IGV, a popular visualization tool for High-Throughput Sequencing (HTS) data. You can run IGV on your local machine.

#### Alignment Visualization - Preparation (Optional)
Before we can view our alignments in the IGV browser, we need to index our BAM files. To do this, we will use samtools index. For convenience, please index all BAM files.
```
echo $BAM_P
cd $BAM_P
samtools index HBR.bam
samtools index HBR_Rep1.bam
samtools index HBR_Rep2.bam
samtools index HBR_Rep3.bam
samtools index UHR.bam
samtools index UHR_Rep1.bam
samtools index UHR_Rep2.bam
samtools index UHR_Rep3.bam
```
#### Alignment Visualization - IGV (Optional)
Download the following files from `/project/biocompworkshop/<YOUR-FOLDER>/BAMFiles_Paired` to your local machine:
- HBR.bam
- HBR.bam.bai
- UHR.bam
- UHR.bam.bai

Once downloaded, open these files in IGV for visualization.

Explore Gene Locus on chr22

Go to an example gene locus on chr22:
- e.g. EIF3L, NDUFA6, and RBX1 have nice coverage
- e.g. SULT4A1 and GTSE1 are differentially expressed. Determine if they are up-regulated or down-regulated in the brain (HBR) compared to cancer cell lines (UHR).

Mouse over some reads and use the read group (RG) flag to determine which replicate the reads come from. Explore other details about each read and its alignment to the reference genome.

Try to find a variant position in the RNAseq data:
- HINT: DDX17 is a highly expressed gene with several variants in its 3' UTR.
- Other highly expressed genes to explore are: NUP50, CYB5R3, and EIF3L (all have at least one transcribed variant).

Are these variants previously known (e.g., present in dbSNP)? How should we interpret the allele frequency of each variant, considering that our samples are pooled RNAs from multiple individuals?

Take note of the genomic position of your variant for future reference.

