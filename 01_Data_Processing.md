# Data Processing

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

### Setup
In this module, we will create the directories for storing input and output files. We will also upload all the programs required for analysis.

In the Classroom Folder: `cd /project/biocompworkshop`  
Navigate to the folder with your name; this is your home folder.  
Example command: `cd <your-folder>`

##### Created the following folders (inside your home folder)

```bash
mkdir BAMFiles_Paired
mkdir FastQ
mkdir FastQC
mkdir Counts
mkdir -p Grch38/{Hisat2,fasta,genes}
```
-   The command `mkdir -p Grch38/{Hisat2,fasta,genes}` creates a nested folder
-   `ls Grch38` and you will find three folders named Hisat2, fasta and genes

##### Create path shortcuts 

```bash
export myHOME=/project/biocompworkshop/rshukla 
export FASTQ=/project/biocompworkshop/rshukla/FastQ
export FASTQC=/project/biocompworkshop/rshukla/FastQC
export REFERENCE=/project/biocompworkshop/rshukla/Grch38/fasta
export GTF=/project/biocompworkshop/rshukla/Grch38/genes
export INDEX=/project/biocompworkshop/rshukla/Grch38/Hisat2
export BAM_P=/project/biocompworkshop/rshukla/BAMFiles_Paired
export COUNTS=/project/biocompworkshop/rshukla/Counts
```
Following the creation of path shortcuts, you can now navigate to the FastQ folder (where we will download our raw data) using `cd $FASTQ`, and then return to your home folder using `cd $myHOME`.

##### Load all modules required for the class
We will be using several Linux-based tools in this section. You can load these tools using the following commands

```bash
module load arcc/1.0
module load gcc/12.2.0
module load samtools/1.16.1
module load hisat2/2.2.1
module load fastqc/0.11.9
module load multiqc/1.13
module load openjdk/11.0.15_10
module load picard/2.26.2
module load bedops/2.4.40
module load kentutils/1.04.0
```

#### Introduction to Inputs
In the lecture, we will delve into the essentials of reference genome files (.fa), gene annotation files (.gtf), and raw sequence files (.fastq). Here we will cover the necessary commands for downloading, understanding, indexing, and performing quality control on these files.

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
**Reference genome:** The nucleotide sequence of the chromosomes of a species. Genes are the functional units of a reference genome and gene annotations describe the structure of transcripts expressed from those gene loci.

**Gene annotations:** - Descriptions of gene/transcript models for a genome. A transcript model consists of the coordinates of the exons of a transcript on a reference genome. Additional information such as the strand the transcript is generated from, gene name, coding portion of the transcript, alternate transcript start sites, and other information may be provided.

**GTF (.gtf) file:** - A common file format referred to as Gene Transfer Format used to store gene and transcript annotation information. You can learn more about this format [here](http://genome.ucsc.edu/FAQ/FAQformat#format4).


#### Reference Genomes

```bash
echo $REFERENCE
cd $REFERENCE
wget http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa
ls
```
View the first 10 lines of this file. Why does it look like this?

```bash
head chr22_with_ERCC92.fa
```
How many lines and characters are in this file? How long is this chromosome (in bases and Mbp)?

```bash
wc chr22_with_ERCC92.fa
```
View 10 lines from approximately the middle of this file. What is the significance of the upper and lower case characters?

```bash
head -n 425000 chr22_with_ERCC92.fa | tail
```
What is the count of each base in the entire reference genome file (skipping the header lines for each sequence)?

The Awk approach involves utilizing Awk, an alternative scripting language included in most Linux distributions. This command is conceptually similar to the Perl approach but with a different syntax. A for loop is employed to iterate over each character until the end ("NF") is reached. The counts for each letter are then stored in a simple data structure, and once the end of the file is reached, the results are printed. 

```bash
time cat chr22_with_ERCC92.fa | grep -v ">" | awk '{for (i=1; i<=NF; i++){a[$i]++}}END{for (i in a){print a[i], i}}' FS= - | sort -k 2 | column -t
```
The Sed approach involves using Sed, an alternative scripting language. First, the "tr" command is used to remove newline characters. Then, Sed is utilized to split each character onto its own line, effectively creating a file with millions of lines. Afterward, Unix sort and uniq commands are applied to produce counts of each unique character, and sort is used to order the results consistently with the previous approaches.

```bash
time cat chr22_with_ERCC92.fa | grep -v ">" | tr -d '\n' | sed 's/\(.\)/\1\n/g'  - | sort | uniq -c | sort -k 2 | column -t
```
The grep approach utilizes the "-o" option of grep to split each match onto a line, which we then use to obtain a count. The "-i" option enables case-insensitive matching, and the "-P" option allows us to use Perl-style regular expressions with grep.

```bash
time cat chr22_with_ERCC92.fa | grep -v ">" | grep -i -o -P "a|c|g|t|y|n" | sort | uniq -c
```
#### Annotation

```bash
echo $GTF
cd $GTF
wget http://genomedata.org/rnaseq-tutorial/annotations/GRCh38/chr22_with_ERCC92.gtf
```
Take a look at the contents of the .gtf file. Press q to exit the less display.

```bash
less -p start_codon -S chr22_with_ERCC92.gtf
```
Note how the -S option makes it easier to veiw this file with less. Make the formatting a bit nicer still:

```bash
cat chr22_with_ERCC92.gtf | column -t | less -p exon -S
```
How many unique gene IDs are in the .gtf file?
We can also use grep to find this same information.

```bash
cat chr22_with_ERCC92.gtf | grep -w gene | wc -l
```

`grep -w gene` is telling grep to do an exact match for the string ‘gene’. This means that it will return lines that are of the feature type `gene`.

Now view the structure of a single transcript in GTF format. Press `q` to exit the `less` display when you are done.

```bash
grep ENST00000342247 $GTF/chr22_with_ERCC92.gtf | less -p "exon\s" -S
```
#### Indexing

```bash
echo $GTF
cd $GTF
hisat2_extract_splice_sites.py chr22_with_ERCC92.gtf > $INDEX/splicesites.tsv
hisat2_extract_exons.py chr22_with_ERCC92.gtf > $INDEX/exons.tsv
hisat2-build -p 4 --ss $INDEX/splicesites.tsv --exon $INDEX/exons.tsv $REFERENCE/chr22_with_ERCC92.fa $INDEX/chr22
```
#### RNAseq Data

```bash
echo $FASTQ
cd $FASTQ
wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar
tar -xvf HBR_UHR_ERCC_ds_5pc.tar
ls
zcat UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz | head -n 8
zcat UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz | grep -P "^\@HWI" | wc -l
```
#### Pre-alignment QC

```bash
cd $FASTQ
fastqc -O $FASTQC/ *.fastq.gz
cd $FASTQC
multiqc ./
```

### Module 02
##### Learning objectives
- RNA-seq alignment challenges and common questions
- Alignment strategies
- HISAT2
- Introduction to the BAM and BED formats
- Basic manipulation of BAMs
- Visualization of RNA-seq alignments in IGV
- Alignment QC Assessment
- BAM read counting and determination of variant allele expression status

#### Adapter Trim (Optional)
Details for adapter trim go here.

#### Alignment
HISAT2 uses a graph-based alignment and has succeeded HISAT and TOPHAT2. The output of this step will be a SAM/BAM file for each data set.

HISAT2 options specified below:
- `-p 4` tells HISAT2 to use 4 CPUs for bowtie alignments.
- `–rna-strandness RF` specifies strandness of RNAseq library. We will specify RF since the TruSeq strand-specific library was used to make these libraries. See here for options.
- `–rg-id $ID` specifies a read group ID that is a unique identifier.
- `–rg SM:$SAMPLE_NAME` specifies a read group sample name. This together with rg-id will allow you to determine which reads came from which sample in the merged bam later on.
- `–rg LB:$LIBRARY_NAME` specifies a read group library name. This together with rg-id will allow you to determine which reads came from which library in the merged bam later on.
- `–rg PL:ILLUMINA` specifies a read group sequencing platform.
- `–rg PU:$PLATFORM_UNIT` specifies a read group sequencing platform unit. Typically this consists of FLOWCELL-BARCODE.LANE
- `–dta` Reports alignments tailored for transcript assemblers.
- `-x /path/to/hisat2/index` The HISAT2 index filename prefix (minus the trailing .X.ht2) built earlier including splice sites and exons.
- `-1 /path/to/read1.fastq.gz` The read 1 FASTQ file, optionally gzip(.gz) or bzip2(.bz2) compressed.
- `-2 /path/to/read2.fastq.gz` The read 2 FASTQ file, optionally gzip(.gz) or bzip2(.bz2) compressed.
- `-S /path/to/output.sam` The output SAM format text file of alignments.

```bash
hisat2 -p 4 --rg-id=UHR_Rep1 --rg SM:UHR --rg LB:UHR_Rep1_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-ACTGAC.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz -S $BAM_P/UHR_Rep1.sam 
hisat2 -p 4 --rg-id=UHR_Rep2 --rg SM:UHR --rg LB:UHR_Rep2_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-TGACAC.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz -S $BAM_P/UHR_Rep2.sam 
hisat2 -p 4 --rg-id=UHR_Rep3 --rg SM:UHR --rg LB:UHR_Rep3_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-CTGACA.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz -S $BAM_P/UHR_Rep3.sam 

hisat2 -p 4 --rg-id=HBR_Rep1 --rg SM:HBR --rg LB:HBR_Rep1_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-TGACAC.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz -S $BAM_P/HBR_Rep1.sam 
hisat2 -p 4 --rg-id=HBR_Rep2 --rg SM:HBR --rg LB:HBR_Rep2_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-GACACT.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz -S $BAM_P/HBR_Rep2.sam 
hisat2 -p 4 --rg-id=HBR_Rep3 --rg SM:HBR --rg LB:HBR_Rep3_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-ACACTG.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz -S $BAM_P/HBR_Rep3.sam 
```
SAM to BAM Conversion Convert HISAT2 sam files to bam files and sort by aligned position

```bash
cd $BAM_P
samtools sort -@ 4 -o UHR_Rep1.bam UHR_Rep1.sam
samtools sort -@ 4 -o UHR_Rep2.bam UHR_Rep2.sam
samtools sort -@ 4 -o UHR_Rep3.bam UHR_Rep3.sam
samtools sort -@ 4 -o HBR_Rep1.bam HBR_Rep1.sam
samtools sort -@ 4 -o HBR_Rep2.bam HBR_Rep2.sam
samtools sort -@ 4 -o HBR_Rep3.bam HBR_Rep3.sam
```
You can pipe `|` the hisat2 output to samtools to get the .bam file

```bash
hisat2 -p 4 --rg-id=UHR_Rep1 --rg SM:UHR --rg LB:UHR_Rep1_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-ACTGAC.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz | samtools view -bS - | samtools sort > $BAM_P/UHR_Rep1.bam 
hisat2 -p 4 --rg-id=UHR_Rep2 --rg SM:UHR --rg LB:UHR_Rep2_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-TGACAC.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz | samtools view -bS - | samtools sort > $BAM_P/UHR_Rep2.bam 
hisat2 -p 4 --rg-id=UHR_Rep3 --rg SM:UHR --rg LB:UHR_Rep3_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-CTGACA.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz | samtools view -bS - | samtools sort > $BAM_P/UHR_Rep3.bam 

hisat2 -p 4 --rg-id=HBR_Rep1 --rg SM:HBR --rg LB:HBR_Rep1_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-TGACAC.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz | samtools view -bS - | samtools sort > $BAM_P/HBR_Rep1.bam 
hisat2 -p 4 --rg-id=HBR_Rep2 --rg SM:HBR --rg LB:HBR_Rep2_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-GACACT.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz | samtools view -bS - | samtools sort > $BAM_P/HBR_Rep2.bam 
hisat2 -p 4 --rg-id=HBR_Rep3 --rg SM:HBR --rg LB:HBR_Rep3_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-ACACTG.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz | samtools view -bS - | samtools sort > $BAM_P/HBR_Rep3.bam 
```
Merge HISAT2 BAM files
Make a single BAM file combining all UHR data and another for all HBR data. This could be done using picard-tools.

```bash
cd $BAM_P
java -Xmx2g -jar $PICARD MergeSamFiles -OUTPUT UHR.bam -INPUT UHR_Rep1.bam -INPUT UHR_Rep2.bam -INPUT UHR_Rep3.bam
java -Xmx2g -jar $PICARD MergeSamFiles -OUTPUT HBR.bam -INPUT HBR_Rep1.bam -INPUT HBR_Rep2.bam -INPUT HBR_Rep3.bam
```
Count the alignment (BAM) files to make sure all were created successfully (you should have 8 total)

```bash
ls -l *.bam | wc -l
ls -l *.bam
```
#### Alignment QC
We will use `samtools flagstat` to get a basic summary of an alignment.

```bash
echo $BAM_P
cd $BAM_P
mkdir flagstat
samtools flagstat HBR_Rep1.bam > flagstat/HBR_Rep1.bam.flagstat
samtools flagstat HBR_Rep2.bam > flagstat/HBR_Rep2.bam.flagstat
samtools flagstat HBR_Rep3.bam > flagstat/HBR_Rep3.bam.flagstat
samtools flagstat UHR_Rep1.bam > flagstat/UHR_Rep1.bam.flagstat
samtools flagstat UHR_Rep2.bam > flagstat/UHR_Rep2.bam.flagstat
samtools flagstat UHR_Rep3.bam > flagstat/UHR_Rep3.bam.flagstat
```
View an example

```bash
cat flagstat/UHR_Rep1.bam.flagstat
```
We can summarize the `flagstat` results using `mutiqc`

```bash
cd flagstat
multiqc ./
```
We can also utilize FastQC to conduct basic quality control (QC) of your BAM file (refer to Pre-alignment QC above). This will provide output similar to when you ran FastQC on your fastq files.

```bash
echo $BAM_P
cd $BAM_P
fastqc UHR_Rep1.bam UHR_Rep2.bam UHR_Rep3.bam HBR_Rep1.bam HBR_Rep2.bam HBR_Rep3.bam
mkdir fastqc
mv *fastqc.html fastqc/
mv *fastqc.zip fastqc/
```
We can summarize the `flagstat` results using `mutiqc`

```bash
cd fastqc
multiqc ./
```
We can use Picard to generate RNA-seq specific quality metrics and figures.

```bash
# First, we can generate the necessary input files for picard.
echo $REFERENCE
cd $REFERENCE

# We will first create a .dict file for our reference
java -jar $PICARD CreateSequenceDictionary -R chr22_with_ERCC92.fa -O chr22_with_ERCC92.dict

echo $GTF
cd $GTF
grep --color=none -i -P "rrna|rrp7a" $GTF/chr22_with_ERCC92.gtf > ref_ribosome.gtf
gff2bed < ref_ribosome.gtf > ref_ribosome.bed
java -jar $PICARD BedToIntervalList -I ref_ribosome.bed -O ref_ribosome.interval_list -SD $REFERENCE/chr22_with_ERCC92.dict

gtfToGenePred -genePredExt chr22_with_ERCC92.gtf chr22_with_ERCC92.ref_flat.txt

cat chr22_with_ERCC92.ref_flat.txt | awk '{print $12"\t"$0}' | cut -d$'\t' -f1-11 > tmp.txt
mv tmp.txt chr22_with_ERCC92.ref_flat.txt

cd $BAM_P
mkdir picard
find *Rep*.bam -exec echo java -jar $PICARD CollectRnaSeqMetrics I={} O=picard/{}.RNA_Metrics REF_FLAT=$GTF/chr22_with_ERCC92.ref_flat.txt STRAND=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=$GTF/ref_ribosome.interval_list \; | sh

cd picard
multiqc ./
```
### Module 03 (Performed in R)

##### Learning objectives 
- Expression estimation for known genes and transcripts

```bash
module use /project/biocompworkshop/software/modules
module load subread/2.0.6
echo $BAM_P
cd $BAM_P

featureCounts -p -a $GTF/chr22_with_ERCC92.gtf -o $COUNTS/featurecounts.txt $BAM_P/*[Rr]ep[123].bam

featureCounts -p -a $GTF/chr22_with_ERCC92.gtf -o $COUNTS/featurecounts.txt $BAM_P/*.bam
```

we will get two files: `featurecounts.txt` & `featurecounts.txt.summary`

```bash
echo $COUNTS
cd $COUNTS

cat featurecounts.txt | less

use `q` to come out of the `cat` command.

cat featurecounts.txt | cut -f1,7- | less
```

##### Learning objectives 
- Expression estimation for known genes and transcripts
- Differential expression methods
- Downstream interpretation of expression and differential estimates

##### Preparation 
- Login to your southpass account with ARCC using this link: https://southpass.arcc.uwyo.edu/pun/sys/dashboard  
- Click on Beartooth Xfce Desktop. This will take you to resource request form.
- 
### Parallel Processing
We are using small number of samples which are prepared to study only chromosoem 22. During actuall analysis you will have several large datasets and aligning one file can take hours and aliging all the files could take days. To speed up the process we will perform parallel processing. 







