# Wyoming Computational Biology Workshop

## Table of Contents
1. [Bulk RNA-Seq](#Bulk-RNA-Seq)
2. [Module 00](#setup)
3. [Module 01](#inputs)
   - [Introduction to Inputs](#introduction-to-inputs)
   - [Reference Genomes](#reference-genomes)
   - [Annotation](#annotation)
   - [Indexing](#indexing)
   - [RNAseq Data](#rnaseq-data)
   - [Pre-alignment QC](#pre-alignment-qc)
4. [Module 02](#alignments)
   - [Introduction to Alignment](#introduction-to-alignment)
   - [Adapter Trim](#adapter-trim)
   - [Alignment](#alignment)
   - [IGV](#igv)
   - [Alignment Visualization - Preparation](#alignment-visualization-preparation)
   - [Alignment Visualization - IGV](#alignment-visualization-igv)
   - [Alignment Visualization - Allele Expression](#alignment-visualization-allele-expression)
   - [Alignment QC](#alignment-qc)

---

## Details

### Bulk RNA-Seq
This is the details section for Bulk RNA-Seq.

### Module 00
In this module, we will create the directories for storing input and output files. We will also upload all the programs required for analysis.

In the Classroom Folder: `cd /project/biocompworkshop`  
Navigate to the folder with your name; this is your home folder.  
Example command: `cd <your-folder>`

##### Created the following folders (inside your home folder)
```
mkdir BAMFiles_Paired
mkdir BAMFiles_UnPaired
mkdir FastQ
mkdir FastQC
mkdir Files
```
-   Create a nested folder `mkdir -p Grch38/{Hisat2,fasta,genes}`
-   `cd Grch38` and you will find three folders named Hisat2, fasta and genes

##### Create path shortcuts 
```
export myHOME=/project/biocompworkshop/rshukla 
export FASTQ=/project/biocompworkshop/rshukla/FastQ
export FASTQC=/project/biocompworkshop/rshukla/FastQC
export REFERENCE=/project/biocompworkshop/rshukla/Grch38/fasta
export GTF=/project/biocompworkshop/rshukla/Grch38/genes
export INDEX=/project/biocompworkshop/rshukla/Grch38/Hisat2
export BAM_P=/project/biocompworkshop/rshukla/BAMFiles_Paired
export BAM_UP=/project/biocompworkshop/rshukla/BAMFiles_UnPaired
```

##### Load all modules required for the class
```
module load samtools/1.10
module load hisat2/2.1.0
module load fastqc/0.11.7
```

### Module 01
In this module, we will delve into the essentials of reference genome files (.fa), gene annotation files (.gtf), and raw sequence files (.fastq). We'll cover the necessary commands for downloading, understanding, indexing, and performing quality control on these files.

#### Introduction to Inputs
- Reference genomes can be downloaded from Ensembl using the following links:
    - Mouse: [link]
    - Human: [link]
- The Annotation files can be downloaded using the following links:
   - Mouse: [link]
   - Human: [link]
- For the sake of time, we are going to perform the analysis using only a single chromosome (chr22) and the ERCC spike-in to make it run faster.
- The input RNAseq (.fastq file) consists of two commercially available RNA samples:
   - Universal Human Reference (UHR): Total RNA from 10 diverse cancer cell lines.
   - Human Brain Reference (HBR): Total RNA from the brains of 23 Caucasians, male and female, of varying ages but mostly 60-80 years old.

#### Reference Genomes
```
cd $REFERENCE
echo $REFERENCE
wget http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa
ls
```
View the first 10 lines of this file. Why does it look like this?
```
head chr22_with_ERCC92.fa
```
How many lines and characters are in this file? How long is this chromosome (in bases and Mbp)?
```
wc chr22_with_ERCC92.fa
```
View 10 lines from approximately the middle of this file. What is the significance of the upper and lower case characters?
```
head -n 425000 chr22_with_ERCC92.fa | tail
```
What is the count of each base in the entire reference genome file (skipping the header lines for each sequence)?

#### Annotation
Details for annotation go here.

#### Indexing
Details for indexing go here.

#### RNAseq Data
Details for RNAseq data go here.

#### Pre-alignment QC
Details for pre-alignment QC go here.

### Module 02
Details for Module 02 go here.

#### Introduction to Alignment
Details for the introduction to alignment go here.

#### Adapter Trim
Details for adapter trim go here.

#### Alignment
Details for alignment go here.

#### IGV
Details for IGV go here.

#### Alignment Visualization - Preparation
Details for alignment visualization (preparation) go here.

#### Alignment Visualization - IGV 
Details for alignment visualization







