# Data Processing

## Introduction

Processing RNAseq data involves several key steps to ensure accurate and meaningful results. The primary inputs include RNAseq data files, reference genomes, and annotation files. These inputs are crucial for mapping reads to the genome and identifying gene expression levels.

In this section, we will:

1. Obtain raw RNAseq reads.
2. Check the quality of the reads.
3. Align the reads to the reference genome.
4. Check the alignment quality.
5. Obtain the expression values (counts) for each feature (gene or transcript).

Before starting, we will set up the file system required for smooth processing and load all the necessary modules for the data processing steps. Additionally, we will cover parallel processing techniques needed when dealing with several large files.

## Setup

Here, we will create directories for storing input and output files and create shortcuts for easy access. Additionally, we will upload all the necessary programs for analysis.

In the Classroom Folder: `cd /project/genomicdatasci`  
Navigate to the folder with your name; this is your home folder.  
Example command: `cd <your-folder>`

#### Created the following folders (inside your home folder)

```bash
mkdir BAMFiles_Paired
mkdir FastQ
mkdir FastQC
mkdir Counts
mkdir -p Grch38/{Hisat2,fasta,genes}
```
-   The command `mkdir -p Grch38/{Hisat2,fasta,genes}` creates a nested folder
-   `ls Grch38` and you will find three folders named Hisat2, fasta and genes

We will create some additional folders that will be required throughout day 2 of the workshop.

```bash
mkdir BatchCorrection_Results
mkdir DE_Results
mkdir Pathway_Results
mkdir PP_Results
mkdir scRNASeq_Results
```

#### Requesting an Interactive Session 

Now we will request an interactive session on the SLURM cluster under the genomicdatasci account. The session will last for up to 4 hours, use the teton partition, and request 4 CPU cores. Once the resources are allocated, you will get an interactive shell prompt on one of the compute nodes where you can run your commands interactively.

```bash
salloc --account=genomicdatasci --time=4:00:00 --mem=32G --cpus-per-task=8
```

#### Create path shortcuts 

We will now create environment variables, which are useful for storing paths, configuration options, and other values you might need to reference multiple times in your scripts or commands. Using environment variables makes your scripts more flexible and easier to maintain.

```bash
export myHOME=/project/genomicdatasci/rshukla 
export FASTQ=/project/genomicdatasci/rshukla/FastQ
export FASTQC=/project/genomicdatasci/rshukla/FastQC
export REFERENCE=/project/genomicdatasci/rshukla/Grch38/fasta
export GTF=/project/genomicdatasci/rshukla/Grch38/genes
export INDEX=/project/genomicdatasci/rshukla/Grch38/Hisat2
export BAM_P=/project/genomicdatasci/rshukla/BAMFiles_Paired
export COUNTS=/project/genomicdatasci/rshukla/Counts
```
Following the creation of path shortcuts, you can now navigate to the FastQ folder (where we will download our raw data) using `cd $FASTQ`, and then return to your home folder using `cd $myHOME`.

#### Load all modules required for the class

We will be using several Linux-based tools in this section. You can load these tools using the following commands

```bash
module load arcc/1.0
module load gcc/14.2.0
module use /project/genomicdatasci/software/modules
module load samtools/1.20
module load hisat2/2.2.1
module load fastqc/0.12.1
module load multiqc/1.24.1
module load openjdk/17.0.8.1_1
module load picard/3.2.0
module load bedops/2.4.41
module load kentutils/1.04.0
module load subread/2.0.6
module avail
```

#### Introduction to Inputs

In the lecture, we will delve into the essentials of reference genome files (.fa), gene annotation files (.gtf), and raw sequence files (.fastq). Here we will cover the necessary commands for downloading, understanding, indexing, and performing quality control on these files.

- Reference genomes can be downloaded from Ensembl using the following links:
    - Mouse: https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/
    - Human: https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/
    - Other: https://useast.ensembl.org/info/data/ftp/index.html
- The Annotation files can be downloaded using the following links:
   - Mouse: https://ftp.ensembl.org/pub/release-112/genbank/mus_musculus/
   - Human: https://ftp.ensembl.org/pub/release-112/genbank/homo_sapiens/
- For the sake of time, we are going to perform the analysis using only a single chromosome (chr22) to make it run faster.
- The input RNAseq (.fastq file) consists of two commercially available RNA samples:
   - Universal Human Reference (UHR): Total RNA from 10 diverse cancer cell lines.
   - Human Brain Reference (HBR): Total RNA from the brains of 23 Caucasians, male and female, of varying ages but mostly 60-80 years old.

Definitions:  
**Reference genome:** The nucleotide sequence of the chromosomes of a species. Genes are the functional units of a reference genome and gene annotations describe the structure of transcripts expressed from those gene loci.

**Gene annotations:** - Descriptions of gene/transcript models for a genome. A transcript model consists of the coordinates of the exons of a transcript on a reference genome. Additional information such as the strand the transcript is generated from, gene name, coding portion of the transcript, alternate transcript start sites, and other information may be provided.

**GTF (.gtf) file:** - A common file format referred to as Gene Transfer Format used to store gene and transcript annotation information. You can learn more about this format [here](http://genome.ucsc.edu/FAQ/FAQformat#format4).

## Reference Genomes

```bash
echo $REFERENCE
cd $REFERENCE
wget http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa
```

Check and verify whether the `chr22_with_ERCC92.fa` file has been successfully downloaded.

```bash
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

```bash
time cat chr22_with_ERCC92.fa | grep -v ">" | awk '{for (i=1; i<=NF; i++){a[$i]++}}END{for (i in a){print a[i], i}}' FS= - | sort -k 2 | column -t
```

## Annotation

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

## Indexing

Indexing is a crucial step before performing sequence alignment, ensuring that the alignment tool can quickly and accurately find the best matches for each read in the reference genome. It involves creating a data structure that allows for efficient mapping of sequencing reads to the reference genome by preprocessing the reference genome to generate an index, enabling quick lookup and alignment of sequences.

```bash
echo $GTF
cd $GTF
```
First, we will extract splice site information from the chr22_with_ERCC92.gtf file and writes it to a file named `splicesites.tsv` in the directory specified by the INDEX environment variable.

```bash
hisat2_extract_splice_sites.py chr22_with_ERCC92.gtf > $INDEX/splicesites.tsv
```

Next, we will extract exon information from the `chr22_with_ERCC92.gtf` file and write it to a file named `exons.tsv` in the directory specified by the INDEX environment variable.

```bash
hisat2_extract_exons.py chr22_with_ERCC92.gtf > $INDEX/exons.tsv
```

Lastly, we will build HISAT2 index for the `chr22_with_ERCC92.fa` reference genome. It uses 4 CPU threads (-p 4), the splice sites information from `splicesites.tsv`, and the exon information from `exons.tsv`. The index files are created in the directory specified by the INDEX environment variable with a base name of chr22.

```bash
hisat2-build -p 4 --ss $INDEX/splicesites.tsv --exon $INDEX/exons.tsv $REFERENCE/chr22_with_ERCC92.fa $INDEX/chr22
```

## RNAseq Data

```bash
echo $FASTQ
cd $FASTQ
wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar
```

Lets extract all files from the HBR_UHR_ERCC_ds_5pc.tar archive and display the names of the files as they are extracted

```bash
tar -xvf HBR_UHR_ERCC_ds_5pc.tar
```

Check and verify whether the files are successfully extracted.

```bash
ls
```

Let's decompress the specified FASTQ file and then display the first 8 lines of the decompressed file.

```bash
zcat UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz | head -n 8
```

Let's decompresses the specified FASTQ file, searche for lines that start with @HWI (which is typical of Illumina read identifiers), and count how many such lines exist. This effectively will count the number of reads in the file that match the specified pattern.

```bash
zcat UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz | grep -P "^\@HWI" | wc -l
```

## Pre-alignment QC

```bash
cd $FASTQ
```

We will run FastQC on all the `.fastq.gz` files in the current directory and save the resulting quality control reports to the directory specified by the `$FASTQC` variable.

```bash
fastqc -O $FASTQC/ *.fastq.gz
```
Next, we will `cd` to `$FASTQC` and run the MultiQC tool on the current directory (`./`). MultiQC is a tool used for aggregating and visualizing results from `fastqc` and other similar tools. 

```bash
cd $FASTQC
multiqc ./
```

## Alignment (to Reference Genome)

Alignment to a reference genome is the process of mapping DNA or RNA sequences to a known reference genome to identify similarities and differences. This crucial step aids in detecting mutations, structural variants, and gene expression levels. Commonly used tools for this purpose include BWA, Bowtie, and HISAT2.

In our focus on alignment, we'll delve into using HISAT2. HISAT2 employs a graph-based alignment approach and has superseded HISAT and TOPHAT2. The steps involved in HISAT2 are similar to those in many other aligners. Upon executing the commands outlined below, each dataset will yield a SAM/BAM file as output.

```bash
hisat2 -p 4 --rg-id=UHR_Rep1 --rg SM:UHR --rg LB:UHR_Rep1_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-ACTGAC.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz -S $BAM_P/UHR_Rep1.sam 
hisat2 -p 4 --rg-id=UHR_Rep2 --rg SM:UHR --rg LB:UHR_Rep2_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-TGACAC.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz -S $BAM_P/UHR_Rep2.sam 
hisat2 -p 4 --rg-id=UHR_Rep3 --rg SM:UHR --rg LB:UHR_Rep3_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-CTGACA.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz -S $BAM_P/UHR_Rep3.sam 

hisat2 -p 4 --rg-id=HBR_Rep1 --rg SM:HBR --rg LB:HBR_Rep1_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-TGACAC.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz -S $BAM_P/HBR_Rep1.sam 
hisat2 -p 4 --rg-id=HBR_Rep2 --rg SM:HBR --rg LB:HBR_Rep2_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-GACACT.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz -S $BAM_P/HBR_Rep2.sam 
hisat2 -p 4 --rg-id=HBR_Rep3 --rg SM:HBR --rg LB:HBR_Rep3_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-ACACTG.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz -S $BAM_P/HBR_Rep3.sam 
```

HISAT2 options specified above:
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
# hisat2 -p 4 --rg-id=UHR_Rep1 --rg SM:UHR --rg LB:UHR_Rep1_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-ACTGAC.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz | samtools view -bS - | samtools sort > $BAM_P/UHR_Rep1.bam 
# hisat2 -p 4 --rg-id=UHR_Rep2 --rg SM:UHR --rg LB:UHR_Rep2_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-TGACAC.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz | samtools view -bS - | samtools sort > $BAM_P/UHR_Rep2.bam 
# hisat2 -p 4 --rg-id=UHR_Rep3 --rg SM:UHR --rg LB:UHR_Rep3_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-CTGACA.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz | samtools view -bS - | samtools sort > $BAM_P/UHR_Rep3.bam 

# hisat2 -p 4 --rg-id=HBR_Rep1 --rg SM:HBR --rg LB:HBR_Rep1_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-TGACAC.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz | samtools view -bS - | samtools sort > $BAM_P/HBR_Rep1.bam 
# hisat2 -p 4 --rg-id=HBR_Rep2 --rg SM:HBR --rg LB:HBR_Rep2_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-GACACT.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz | samtools view -bS - | samtools sort > $BAM_P/HBR_Rep2.bam 
# hisat2 -p 4 --rg-id=HBR_Rep3 --rg SM:HBR --rg LB:HBR_Rep3_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-ACACTG.1 -x $INDEX/chr22 --dta --rna-strandness RF -1 $FASTQ/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $FASTQ/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz | samtools view -bS - | samtools sort > $BAM_P/HBR_Rep3.bam 
```

## Merge HISAT2 BAM files

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

## Alignment QC

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

We can summarize the `fastqc` results using `mutiqc`

```bash
cd fastqc
multiqc ./
```

We can use Picard to generate RNA-seq specific quality metrics and figures. First, we need to create the necessary input files, starting with a `.dict` file for our reference.

- A `.dict` file, or sequence dictionary file, is used in bioinformatics workflows, particularly those involving tools like Picard and GATK (Genome Analysis Toolkit). It contains metadata about the reference genome sequences, such as the names and lengths of each chromosome or contig, and their order. This information helps tools to interpret and process genomic data correctly.

```bash
echo $REFERENCE
cd $REFERENCE

java -jar $PICARD CreateSequenceDictionary -R chr22_with_ERCC92.fa -O chr22_with_ERCC92.dict
```

We will now create a BED file of the location of ribosomal sequences in our reference by first extracting them from the GTF and then converting to BED. Note that here we pull all the "rrna" transcripts from the GTF. This strategy is suitable for the whole transcriptome, but on chr22, there is very little "rrna" content, leading to 0 coverage for all samples. Therefore, we are also adding a single protein-coding ribosomal gene "RRP7A" (normally, we would not do this).

- A BED (Browser Extensible Data) file is used to store genomic coordinates of features.

```bash
echo $GTF
cd $GTF
grep --color=none -i -P "rrna|rrp7a" $GTF/chr22_with_ERCC92.gtf > ref_ribosome.gtf
gff2bed < ref_ribosome.gtf > ref_ribosome.bed
```

Create interval list file for the location of ribosomal sequences in our reference

```bash
java -jar $PICARD BedToIntervalList -I ref_ribosome.bed -O ref_ribosome.interval_list -SD $REFERENCE/chr22_with_ERCC92.dict
```

Create a genePred file for our reference transcriptome.

- A GENPHRED file stores quality scores for each base in a sequence, representing the accuracy of the base calling.

```bash
gtfToGenePred -genePredExt chr22_with_ERCC92.gtf chr22_with_ERCC92.ref_flat.txt
```

Reformat the genePred file

```bash
cat chr22_with_ERCC92.ref_flat.txt | awk '{print $12"\t"$0}' | cut -d$'\t' -f1-11 > tmp.txt
mv tmp.txt chr22_with_ERCC92.ref_flat.txt
```

Finally, we will use picard tools to collect RNA-seq metrics for each BAM file and save the results in the picard directory.

```bash
cd $BAM_P
mkdir picard
find *Rep*.bam -exec echo java -jar $PICARD CollectRnaSeqMetrics I={} O=picard/{}.RNA_Metrics REF_FLAT=$GTF/chr22_with_ERCC92.ref_flat.txt STRAND=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=$GTF/ref_ribosome.interval_list \; | sh
```

We can summarize and generate figures for alignment QC using multiqc.

```bash
cd picard
multiqc ./
```
 
## Expression estimation for known genes and transcripts

In the final step of our HPC process, we will conduct read counting on non-merged paired-end BAM files, utilizing the genomic annotation from the provided GTF file (`chr22_with_ERCC92.gtf`). The count results will be saved to `featurecounts.txt` within the `counts` directory. The `[Rr]ep[123].bam` regex pattern selectively identifies the unmerged BAM files based on their filenames.

```bash
echo $BAM_P
cd $BAM_P

featureCounts -p -a $GTF/chr22_with_ERCC92.gtf -o $COUNTS/featurecounts.txt $BAM_P/*[Rr]ep[123].bam
```

we will get two files: `featurecounts.txt` & `featurecounts.txt.summary` in the `counts` directory.

```bash
echo $COUNTS
cd $COUNTS

cat featurecounts.txt | less
```

Use `q` to come out of the `cat` command.

To enhance readability, we'll extract the first field and fields starting from the seventh onwards from each line. This modified content will then be displayed interactively using the `less` pager for better viewing.

```bash
cat featurecounts.txt | cut -f1,7- | less
```
 
## Parallel Processing

We are working with a small number of samples specifically prepared to study chromosome 22. In a real-world analysis, you will typically have much larger datasets, and aligning a single file can take hours, while aligning all files could take days. To expedite this process, we will employ parallel processing.

We'll use the same example dataset and run it using parallel processing. First, we'll create a script that generates the `hisat2` alignment commands for all the files you need to align. These commands will be saved in a file named `cmd.list`, which we will then use to run all the alignment commands in parallel.

#### Generating list of hisat2 commands

You can use the following script by pasting it into any text editor (vi or nano) and saving it as a .bash file. For this session, it has been provided as `hisat2_cmd.bash` in the `Data_Vault` folder. This script automates the generation of alignment commands for RNA-seq data using hisat2 and samtools.

```bash
#!/bin/bash

idx_dir="/project/genomicdatasci/rshukla/Grch38/Hisat2/chr22"
splice_dir="/project/genomicdatasci/rshukla/Grch38/Hisat2/splicesites.tsv"
SRC_DIR="/project/genomicdatasci/rshukla/FastQ"
out_dir="/project/genomicdatasci/rshukla/PP_Results"

cd $SRC_DIR

for file in *read1.fastq.gz; 
do
Base_Name="${file%read1.fastq.gz}"
out_file="${file%_ERCC*}"
out_name="${out_dir}/${out_file}.bam"

L1R1_File="${Base_Name}read1.fastq.gz"
L2R1_File="${Base_Name}read1.fastq.gz"
L1R2_File="${Base_Name}read2.fastq.gz"
L2R2_File="${Base_Name}read2.fastq.gz"

echo "hisat2 -p 4 --rg-id=${out_file} --rg SM:${out_file} --rg LB:${out_file} --rg PL:ILLUMINA --rg PU:CXX1234-ACTGAC.1 -x ${idx_dir} --dta --rna-strandness RF -1 ${SRC_DIR}/${L1R1_File} -2 ${SRC_DIR}/${L1R2_File} | samtools view -bS - | samtools sort > ${out_name}"

done
```
Here's a brief breakdown of what each part of the script does:

1. **Shebang and Variable Initialization:**
   ```bash
   #!/bin/bash
   idx_dir="/project/genomicdatasci/rshukla/Grch38/Hisat2/chr22"
   splice_dir="/project/genomicdatasci/rshukla/Grch38/Hisat2/splicesites.tsv"
   SRC_DIR="/project/genomicdatasci/rshukla/FastQ"
   out_dir="/project/genomicdatasci/rshukla/PP_Results"
   ```
   - The script starts with the shebang (`#!/bin/bash`), indicating that it should be run in the Bash shell.
   - It initializes several variables for directories: the index directory (`idx_dir`), splice site file (`splice_dir`), source directory containing FASTQ files (`SRC_DIR`), and output directory (`out_dir`).

2. **Navigate to the Source Directory:**
   ```bash
   cd $SRC_DIR
   ```
   - Changes the current working directory to the source directory containing the FASTQ files.

3. **Loop Through FASTQ Files:**
   ```bash
   for file in *read1.fastq.gz; 
   do
   ```
   - Loops through all files in the source directory that match the pattern `*read1.fastq.gz`.

4. **Generate Base Names and Output File Names:**
   ```bash
   Base_Name="${file%read1.fastq.gz}"
   out_file="${file%_ERCC*}"
   out_name="${out_dir}/${out_file}.bam"
   ```
   - Extracts the base name of each file by removing the `read1.fastq.gz` suffix.
   - Creates the output file name by removing the `_ERCC` suffix from the original file name.
   - Constructs the full path for the output BAM file.

5. **Set Read Pair File Names:**
   ```bash
   L1R1_File="${Base_Name}read1.fastq.gz"
   L2R1_File="${Base_Name}read1.fastq.gz"
   L1R2_File="${Base_Name}read2.fastq.gz"
   L2R2_File="${Base_Name}read2.fastq.gz"
   ```
   - Defines the file names for the read pairs (assuming a specific naming convention).

6. **Print `hisat2` Command:**
   ```bash
   echo "hisat2 -p 4 --rg-id=${out_file} --rg SM:${out_file} --rg LB:${out_file} --rg PL:ILLUMINA --rg PU:CXX1234-ACTGAC.1 -x ${idx_dir} --dta --rna-strandness RF -1 ${SRC_DIR}/${L1R1_File} -2 ${SRC_DIR}/${L1R2_File} | samtools view -bS - | samtools sort > ${out_name}"
   ```
   - Constructs a `hisat2` command for aligning the reads, specifying various options:
     - `-p 4`: Use 4 threads.
     - `--rg-id`, `--rg SM`, `--rg LB`, `--rg PL`, `--rg PU`: Include read group information.
     - `-x ${idx_dir}`: Specify the index directory.
     - `--dta`: Enable downstream transcriptome assembly.
     - `--rna-strandness RF`: Specify RNA strandness.
     - `-1`, `-2`: Specify the input read files.
   - Pipes the output to `samtools` for converting to BAM format (`view -bS -`) and sorting (`samtools sort`).
   - Prints the constructed command to the terminal (could be redirected to a file if needed).

7. **End Loop:**
   ```bash
   done
   ```
   - Ends the loop.

### Execute the script

To execute the `hisat2_cmd.bash` script, we will first copy it from `Data_Vault` to our home directory and then execute it. I have changed the permissions, but if it has not been changed, you can use `chmod +x hisat2_cmd.bash` to do so.

```bash
echo $myHOME
cd $myHOME

cp /project/genomicdatasci/Data_Vault/hisat2_cmd.bash .

ls

./hisat2_cmd.bash > cmd.list
```
### Run the command list in parallel
Below I have given the Slurm batch script designed to run commands listed in a file (cmd.list) in parallel using GNU Parallel.

```bash
#!/bin/bash
#SBATCH --job-name=hisat2-cmd
#SBATCH --account=genomicdatasci
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=6
#SBATCH --time=01:00:00
 
module load arcc/1.0  
module load gcc/12.2.0
module load hisat2/2.2.1
module load samtools/1.20
module load parallel/20220522

srun parallel --jobs 8 < /project/genomicdatasci/Data_Vault/cmd.list
```
The script does the following:
- Allocates resources on an HPC cluster (2 nodes, 4 tasks per node, 6 CPUs per task).
- Loads the necessary software modules required for the tasks.
- Uses GNU Parallel to execute commands listed in /project/genomicdatasci/Data_Vault/cmd.list concurrently, with up to 8 commands running at the same time.

### Execute the script

To execute the `ParallelProcessing.sh` script, we will first copy it from `Data_Vault` to our home directory and then execute it. Again, I have changed the permissions, but if it has not been changed, you can use `chmod +x ParallelProcessing.sh` to do so.

```bash
echo $myHOME
cd $myHOME

cp /project/genomicdatasci/Data_Vault/ParallelProcessing.sh .

ls

sbatch ParallelProcessing.sh
```

The `sbatch ParallelProcessing.sh` command will create log file named `slurm-xxxxxxxx.out` file which will have all details of the commands you are running. Where `x` represents some number. You can view this file using `cat` command. 

```bash
cat slurm-xxxxxxxx.out
```

To see the status of jobs in the queue, use the squeue command:

View All Jobs:

```bash
squeue
```

View Jobs for a Specific User:

```bash
squeue -u rshukla
```

Detailed Job Information:

```bash
scontrol show job <job_id>
```

To terminate a job that is currently running or queued, use the scancel command:

```bash
scancel job_id
```

Terminate All Jobs for a Specific User:

```bash
scancel -u rshukla
```
### Additional files required for the workshop

```bash
echo $myHOME
cd $myHOME

cp /project/genomicdatasci/Data_Vault/R_Scripts/*.Rmd .
```

### Starting R session for the workshop

```bash
module use /project/genomicdatasci/software/modules
module load r/4.4.0-biocomp
module load rstudio/2023.9.0
```
To use the R library containing the various R packages you will need to call the following from the R Console (within RStudio) and/or at the top of your R scripts:

```R
.libPaths(c('/project/genomicdatasci/software/r/libraries/4.4.0', '/apps/u/spack/gcc/12.2.0/r/4.4.0-7i7afpk/rlib/R/library'))
```

Confirm the library path:

```R
.libPaths()
```






