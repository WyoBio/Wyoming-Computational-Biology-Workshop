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
module load multiqc/1.13
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

Definitions:
Reference genome - The nucleotide sequence of the chromosomes of a species. Genes are the functional units of a reference genome and gene annotations describe the structure of transcripts expressed from those gene loci.

Gene annotations - Descriptions of gene/transcript models for a genome. A transcript model consists of the coordinates of the exons of a transcript on a reference genome. Additional information such as the strand the transcript is generated from, gene name, coding portion of the transcript, alternate transcript start sites, and other information may be provided.

GTF (.gtf) file - A common file format referred to as Gene Transfer Format used to store gene and transcript annotation information. You can learn more about this format here: http://genome.ucsc.edu/FAQ/FAQformat#format4

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

The Awk approach involves utilizing Awk, an alternative scripting language included in most Linux distributions. This command is conceptually similar to the Perl approach but with a different syntax. A for loop is employed to iterate over each character until the end ("NF") is reached. The counts for each letter are then stored in a simple data structure, and once the end of the file is reached, the results are printed. 
```
time cat chr22_with_ERCC92.fa | grep -v ">" | awk '{for (i=1; i<=NF; i++){a[$i]++}}END{for (i in a){print a[i], i}}' FS= - | sort -k 2 | column -t
```
The Sed approach involves using Sed, an alternative scripting language. First, the "tr" command is used to remove newline characters. Then, Sed is utilized to split each character onto its own line, effectively creating a file with millions of lines. Afterward, Unix sort and uniq commands are applied to produce counts of each unique character, and sort is used to order the results consistently with the previous approaches.
```
time cat chr22_with_ERCC92.fa | grep -v ">" | tr -d '\n' | sed 's/\(.\)/\1\n/g'  - | sort | uniq -c | sort -k 2 | column -t
```
The grep approach utilizes the "-o" option of grep to split each match onto a line, which we then use to obtain a count. The "-i" option enables case-insensitive matching, and the "-P" option allows us to use Perl-style regular expressions with grep.
```
time cat chr22_with_ERCC92.fa | grep -v ">" | grep -i -o -P "a|c|g|t|y|n" | sort | uniq -c
```
#### Annotation
```
cd $GTF
echo $GTF
wget http://genomedata.org/rnaseq-tutorial/annotations/GRCh38/chr22_with_ERCC92.gtf
```
Take a look at the contents of the .gtf file. Press q to exit the less display.
```
less -p start_codon -S chr22_with_ERCC92.gtf
```
Note how the -S option makes it easier to veiw this file with less. Make the formatting a bit nicer still:
```
cat chr22_with_ERCC92.gtf | column -t | less -p exon -S
```
How many unique gene IDs are in the .gtf file?
We can also use grep to find this same information.
```
cat chr22_with_ERCC92.gtf | grep -w gene | wc -l
```
`grep -w gene` is telling grep to do an exact match for the string ‘gene’. This means that it will return lines that are of the feature type `gene`.

Now view the structure of a single transcript in GTF format. Press `q` to exit the `less` display when you are done.
```
grep ENST00000342247 $RNA_REF_GTF | less -p "exon\s" -S
```
#### Indexing
```
cd $GTF
echo $GTF
hisat2_extract_splice_sites.py chr22_with_ERCC92.gtf > $INDEX/splicesites.tsv
hisat2_extract_exons.py chr22_with_ERCC92.gtf > $INDEX/exons.tsv
hisat2-build -p 4 --ss $INDEX/splicesites.tsv --exon $INDEX/exons.tsv $REFERENCE/chr22_with_ERCC92.fa $INDEX/chr22
```
#### RNAseq Data
```
echo $FASTQ
cd $FASTQ
wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar
tar -xvf HBR_UHR_ERCC_ds_5pc.tar
ls
zcat UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz | head -n 8
zcat UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz | grep -P "^\@HWI" | wc -l
```
#### Pre-alignment QC
```
cd $FASTQ
fastqc -O $FASTQC/ *.fastq.gz
cd $FASTQC
multiqc ./
```

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







