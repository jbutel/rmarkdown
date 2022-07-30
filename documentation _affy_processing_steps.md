# Table of contents  

- [Table of contents](#table-of-contents)
- [Software used](#software-used)
- [General processing overview with example commands](#general-processing-overview-with-example-commands)
  - [1. Import Raw Data](#1-import-raw-data)
  - [2. Build STAR Reference](#2-build-star-reference)
  - [4. Align Reads to Reference Genome then Sort and Index](#4-align-reads-to-reference-genome-then-sort-and-index)
    - [4a. Align Reads to Reference Genome with STAR](#4a-align-reads-to-reference-genome-with-star)
    - [4b. Compile Alignment Logs](#4b-compile-alignment-logs)
    - [4c. Tablulate STAR Counts in R](#4c-tablulate-star-counts-in-r)
    - [4d. Sort Aligned Reads](#4d-sort-aligned-reads)
    - [4e. Index Sorted Aligned Reads](#4e-index-sorted-aligned-reads)
  - [5. Create Reference BED File](#5-create-reference-bed-file)
    - [5a. Convert GTF to genePred File](#5a-convert-gtf-to-genepred-file)
    - [5b. Convert genePred to BED File](#5b-convert-genepred-to-bed-file)
  - [6. Assess Strandedness, GeneBody Coverage, Inner Distance, and Read Distribution with RSeQC](#6-assess-strandedness-genebody-coverage-inner-distance-and-read-distribution-with-rseqc)
    - [6a. Determine Read Strandedness](#6a-determine-read-strandedness)
    - [6b. Compile Strandedness Reports](#6b-compile-strandedness-reports)
    - [6c. Evaluate GeneBody Coverage](#6c-evaluate-genebody-coverage)
    - [6d. Compile GeneBody Coverage Reports](#6d-compile-genebody-coverage-reports)
    - [6e. Determine Inner Distance (For Paired End Datasets ONLY)](#6e-determine-inner-distance-for-paired-end-datasets-only)
    - [6f. Compile Inner Distance Reports](#6f-compile-inner-distance-reports)
    - [6g. Assess Read Distribution](#6g-assess-read-distribution)
    - [6h. Compile Read Distribution Reports](#6h-compile-read-distribution-reports)
  - [7. Build RSEM Reference](#7-build-rsem-reference)
  - [8. Quantitate Aligned Reads](#8-quantitate-aligned-reads)
    - [8a. Count Aligned Reads with RSEM](#8a-count-aligned-reads-with-rsem)
    - [8b. Compile RSEM Count Logs](#8b-compile-rsem-count-logs)
    - [8c. Calculate Total Number of Genes Expressed Per Sample in R](#8c-calculate-total-number-of-genes-expressed-per-sample-in-r)
  - [9. Normalize Read Counts, Perform Differential Gene Expression Analysis, and Add Gene Annotations in R](#9-normalize-read-counts-perform-differential-gene-expression-analysis-and-add-gene-annotations-in-r)
    - [9a. For Datasets With ERCC Spike-In](#9a-for-datasets-with-ercc-spike-in)
  - [Install R packages if not already installed](#install-r-packages-if-not-already-installed)
  - [Install annotation R packages if not already installed - only the annotation package for the organism that the data were derived from is required](#install-annotation-r-packages-if-not-already-installed---only-the-annotation-package-for-the-organism-that-the-data-were-derived-from-is-required)
        - [Set up your environment](#set-up-your-environment)
  - [Import libraries (tximport, DESeq2, tidyverse, Risa)](#import-libraries-tximport-deseq2-tidyverse-risa)

---

# Software used  

|Program|Version|Relevant Links|
|:------|:------:|:-------------|
|FastQC|0.11.9|[https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)|
|MultiQC|1.12|[https://multiqc.info/](https://multiqc.info/)|
|Cutadapt|3.7|[https://cutadapt.readthedocs.io/en/stable/](https://cutadapt.readthedocs.io/en/stable/)|
|TrimGalore!|0.6.7|[https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)|
|STAR|2.7.10a|[https://github.com/alexdobin/STAR](https://github.com/alexdobin/STAR)|
|RSEM|1.3.1|[https://github.com/deweylab/RSEM](https://github.com/deweylab/RSEM)|
|Samtools|1.15|[http://www.htslib.org/](http://www.htslib.org/)|
|gtfToGenePred|377|[http://hgdownload.cse.ucsc.edu/admin/exe/](http://hgdownload.cse.ucsc.edu/admin/exe/)|
|genePredToBed|377|[http://hgdownload.cse.ucsc.edu/admin/exe/](http://hgdownload.cse.ucsc.edu/admin/exe/)|
|infer_experiment|4.0.0|[http://rseqc.sourceforge.net/#infer-experiment-py](http://rseqc.sourceforge.net/#infer-experiment-py)|
|geneBody_coverage|4.0.0|[http://rseqc.sourceforge.net/#genebody-coverage-py](http://rseqc.sourceforge.net/#genebody-coverage-py)|
|inner_distance|4.0.0|[http://rseqc.sourceforge.net/#inner-distance-py](http://rseqc.sourceforge.net/#inner-distance-py)|
|read_distribution|4.0.0|[http://rseqc.sourceforge.net/#read-distribution-py](http://rseqc.sourceforge.net/#read-distribution-py)|
|R|4.1.3|[https://www.r-project.org/](https://www.r-project.org/)|
|Bioconductor|3.14.0|[https://bioconductor.org](https://bioconductor.org)|
|DESeq2|1.34|[https://bioconductor.org/packages/release/bioc/html/DESeq2.html](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)|
|tximport|1.22|[https://bioconductor.org/packages/release/bioc/html/tximport.html](https://bioconductor.org/packages/release/bioc/html/tximport.html)|
|tidyverse|1.3.1|[https://www.tidyverse.org](https://www.tidyverse.org)|
|Risa|1.36|[https://www.bioconductor.org/packages/release/bioc/html/Risa.html](https://www.bioconductor.org/packages/release/bioc/html/Risa.html)|
|STRINGdb|2.6.0|[https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html](https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html)|
|PANTHER.db|1.0.11|[https://bioconductor.org/packages/release/data/annotation/html/PANTHER.db.html](https://bioconductor.org/packages/release/data/annotation/html/PANTHER.db.html)|
|org.Hs.eg.db|3.14.0|[https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)|
|Oligo|3.14.0|[https://www.bioconductor.org/packages/release/bioc/html/oligo.html](https://www.bioconductor.org/packages/release/bioc/html/oligo.html)|

---

# General processing overview with example commands  

> Exact processing commands for specific datasets are provided in the [GLDS_Processing_Scripts](GLDS_Processing_Scripts) sub-directory.
> 
> All output files marked with a \# are published for each RNAseq processed dataset in the [GLDS repository](https://genelab-data.ndc.nasa.gov/genelab/projects). 

---

## 1. Import Raw Data


```bash
csv_file <- read.table(file = "<runsheet.csv>", 
                       sep=",",
                       header = TRUE,
                       check.names = FALSE
                       )
array_data_files <- csv_file["Path to Raw Data File"] 

raw_data <- oligo::read.celfiles(array_data_files)
```

**Parameter Definitions:**

- `file` – path of runsheet file to be accessed
- `sep` - setting delimiter, default is white space
- `header` - determines if first row in dataframe is column names
- `check.names` - if TRUE, variable names are checked and corrected for synaticially validity
- `oligo::read.celfiles()` - reads CEL files for Affymetrix microarrays, argument is path to CEL files

**Input Data:**

- *.csv (runsheet)

**Output Data:**

- ExpressionFeatureSet (read CEL files)

---

## 2. Build STAR Reference  

```bash
density_plot <- hist(raw_data, transfo=log2, which=c("pm", "mm", "bg", "both", "all"), nsample=10000, target = "core", main = "Density Plot")

```

**Parameter Definitions:**

- `transfo` – used to scale the data
- `which` - defines specific probe types
- `nsample` - sample size used to produce plot
- `target` - specifies group of meta-probeset
- `main` - main title of plot


**Input Data:**

- *.fasta ([genome sequence](../GeneLab_Reference_and_Annotation_Files/GL-DPPD-7101-E_ensembl_refs.csv))
- *.gtf ([genome annotation](../GeneLab_Reference_and_Annotation_Files/GL-DPPD-7101-E_ensembl_refs.csv))

**Output Data:**

STAR genome reference, which consists of the following files:

- chrLength.txt
- chrNameLength.txt
- chrName.txt
- chrStart.txt
- exonGeTrInfo.tab
- exonInfo.tab
- geneInfo.tab
- Genome
- genomeParameters.txt
- SA
- SAindex
- sjdbInfo.txt
- sjdbList.fromGTF.out.tab
- sjdbList.out.tab
- transcriptInfo.tab

<br>

---

## 4. Align Reads to Reference Genome then Sort and Index

<br>

### 4a. Align Reads to Reference Genome with STAR

```bash
STAR --twopassMode Basic \
 --limitBAMsortRAM 65000000000 \
 --genomeDir /path/to/STAR/genome/directory \
 --outSAMunmapped Within \
 --outFilterType BySJout \
 --outSAMattributes NH HI AS NM MD MC \
 --outFilterMultimapNmax 20 \
 --outFilterMismatchNmax 999 \
 --outFilterMismatchNoverReadLmax 0.04 \
 --alignIntronMin 20 \
 --alignIntronMax 1000000 \
 --alignMatesGapMax 1000000 \ # for PE only
 --alignSJoverhangMin 8 \
 --alignSJDBoverhangMin 1 \
 --sjdbScore 1 \
 --readFilesCommand zcat \
 --runThreadN NumberOfThreads \
 --outSAMtype BAM SortedByCoordinate \
 --quantMode TranscriptomeSAM GeneCounts \
 --outSAMheaderHD @HD VN:1.4 SO:coordinate \
 --outFileNamePrefix /path/to/STAR/output/directory/<sample_id> \
 --readFilesIn /path/to/trimmed_forward_reads \
 /path/to/trimmed_reverse_reads # only needed for PE studies

```

**Parameter Definitions:**

- `--twopassMode` – specifies 2-pass mapping mode; the `Basic` option instructs STAR to perform the 1st pass mapping, then automatically extract junctions, insert them into the genome index, and re-map all reads in the 2nd mapping pass
- `--limitBAMsortRAM` - maximum RAM available (in bytes) to sort the bam files, the example above indicates 65GB
- `--genomeDir` - specifies the path to the directory where the STAR reference is stored
- `--outSAMunmapped` - specifies ouput of unmapped reads in the sam format; the `Within` option instructs STAR to output the unmapped reads within the main sam file
- `--outFilterType` - specifies the type of filtering; the `BySJout` option instructs STAR to keep only those reads that contain junctions that passed filtering in the SJ.out.tab output file
- `--outSAMattributes` - list of desired sam attributes in the order desired for the output sam file; sam attribute descriptions can be found [here](https://samtools.github.io/hts-specs/SAMtags.pdf)
- `--outFilterMultimapNmax` – specifies the maximum number of loci the read is allowed to map to; all alignments will be output only if the read maps to no more loci than this value
- `--outFilterMismatchNmax` - maximum number of mismatches allowed to be included in the alignment output
- `--outFilterMismatchNoverReadLmax` - ratio of mismatches to read length allowed to be included in the alignment output; the `0.04` value indicates that up to 4 mismatches are allowed per 100 bases
- `--alignIntronMin` - minimum intron size; a genomic gap is considered an intron if its length is equal to or greater than this value, otherwise it is considered a deletion
- `--alignIntronMax` - maximum intron size
- `--alignMatesGapMax` - maximum genomic distance (in bases) between two mates of paired-end reads; this option should be removed for single-end reads
- `--alignSJoverhangMin` - minimum overhang (i.e. block size) for unannotated spliced alignments
- `--alignSJDBoverhangMin` - minimum overhang (i.e. block size) for annotated spliced alignments
- `--sjdbScore` - additional alignment score for alignments that cross database junctions
- `--readFilesCommand` - specifies command needed to interpret input files; the `zcat` option indicates input files are compressed with gzip and zcat will be used to uncompress the gzipped input files
- `--runThreadN` - indicates the number of threads to be used for STAR alignment and should be set to the number of available cores on the server node
- `--outSAMtype` - specifies desired output format; the `BAM SortedByCoordinate` options specify that the output file will be sorted by coordinate and be in the bam format
- `--quantMode` - specifies the type(s) of quantification desired; the `TranscriptomeSAM` option instructs STAR to output a separate sam/bam file containing alignments to the transcriptome and the `GeneCounts` option instructs STAR to output a tab delimited file containing the number of reads per gene
- `--outSAMheaderHD` - indicates a header line for the sam/bam file
- `--outFileNamePrefix` - specifies the path to and prefix for the output file names; for GeneLab the prefix is the sample id
- `--readFilesIn` - path to input read 1 (forward read) and read 2 (reverse read); for paired-end reads, read 1 and read 2 should be separated by a space; for single-end reads only read 1 should be indicated

**Input Data:**

- STAR genome reference (output from [Step 3](#3-build-star-reference))
- *fastq.gz (trimmed reads, output from [Step 2a](#2a-trimfilter-raw-data))

**Output Data:**

- *Aligned.sortedByCoord.out.bam (sorted mapping to genome)
- *Aligned.toTranscriptome.out.bam\# (sorted mapping to transcriptome)
- *Log.final.out\# (log file containing alignment info/stats such as reads mapped, etc)
- *ReadsPerGene.out.tab (tab deliminated file containing STAR read counts per gene with 4 columns that correspond to different strandedness options: column 1 = gene ID, column 2 = counts for unstranded RNAseq, column 3 = counts for 1st read strand aligned with RNA, column 4 = counts for 2nd read strand aligned with RNA)
- *Log.out
- *Log.progress.out
- *SJ.out.tab\# (high confidence collapsed splice junctions in tab-delimited format)
- *_STARgenome (directory containing the following:)
  - sjdbInfo.txt
  - sjdbList.out.tab
- *_STARpass1 (directory containing the following:)
  - Log.final.out
  - SJ.out.tab
- *_STARtmp (directory containing the following:)
  - BAMsort (directory containing subdirectories that are empty – this was the location for temp files that were automatically removed after successful completion)

<br>

### 4b. Compile Alignment Logs

```bash
multiqc --interactive -n align_multiqc -o /path/to/aligned_multiqc/output/directory /path/to/*Log.final.out/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*Log.final.out/files` – the directory holding the *Log.final.out output files from the [STAR alignment step](#4a-align-reads-to-reference-genome-with-star), provided as a positional argument

**Input Data:**

- *Log.final.out (log file conting alignment info/stats such as reads mapped, etc., output from [Step 4a](#4a-align-reads-to-reference-genome-with-star))

**Output Data:**

- align_multiqc.html\# (multiqc report)
- /align_multiqc_data\# (directory containing multiqc data)

<br>

### 4c. Tablulate STAR Counts in R

```R
print("Make STAR counts table")
print("")

work_dir="/path/to/working/directory/where/script/is/executed/from" ## Must contain samples.txt file
align_dir="/path/to/directory/containing/STAR/counts/files"

setwd(file.path(work_dir))

### Pull in sample names where the "samples.txt" file is a single column list of sample names ###
study <- read.csv(Sys.glob(file.path(work_dir,"samples.txt")), header = FALSE, row.names = 1, stringsAsFactors = TRUE)

##### Import Data
ff <- list.files(file.path(align_dir), pattern = "ReadsPerGene.out.tab", recursive=TRUE, full.names = TRUE)

## Reorder the *genes.results files to match the ordering of the ISA samples
ff <- ff[sapply(rownames(study), function(x)grep(paste0(x,'_ReadsPerGene.out.tab$'), ff, value=FALSE))]

# Remove the first 4 lines
counts.files <- lapply( ff, read.table, skip = 4 )

# Get counts aligned to either strand for unstranded data by selecting col 2, to the first (forward) strand by selecting col 3 or to the second (reverse) strand by selecting col 4
counts <- as.data.frame( sapply( counts.files, function(x) x[ , 3 ] ) )

# Add column and row names
colnames(counts) <- rownames(study)
row.names(counts) <- counts.files[[1]]$V1


##### Export unnormalized counts table
setwd(file.path(align_dir))
write.csv(counts,file='STAR_Unnormalized_Counts.csv')


## print session info ##
print("Session Info below: ")
print("")
sessionInfo()
```

**Input Data:**

- samples.txt (A newline delimited list of sample IDs)
- *ReadsPerGene.out.tab (STAR counts per gene, output from [Step 4a](#4a-align-reads-to-reference-genome-with-star))

**Output Data:**

- STAR_Unnormalized_Counts.csv\# (Table containing raw STAR counts for each sample)

<br>

### 4d. Sort Aligned Reads

```bash
samtools sort -m 3G \
	--threads NumberOfThreads \
	-o /path/to/*Aligned.sortedByCoord_sorted.out.bam \
  /path/to/*Aligned.sortedByCoord.out.bam
```

**Parameter Definitions:**

- `-m` - memory available per thread, `3G` indicates 3 gigabytes, this can be changed based on user resources
- `--threads` - number of threads available on server node to sort genome alignment files
- `/path/to/*Aligned.sortedByCoord.out.bam` – path to the *Aligned.sortedByCoord.out.bam output files from the [STAR alignment step](#4a-align-reads-to-reference-genome-with-star), provided as a positional argument

**Input Data:**

- *Aligned.sortedByCoord.out.bam (sorted mapping to genome file, output from [Step 4a](#4a-align-reads-to-reference-genome-with-star))

**Output Data:**

- *Aligned.sortedByCoord_sorted.out.bam\# (samtools sorted genome aligned bam file)

<br>

### 4e. Index Sorted Aligned Reads

```bash
samtools index -@ NumberOfThreads /path/to/*Aligned.sortedByCoord_sorted.out.bam
```

**Parameter Definitions:**

- `-@` - number of threads available on server node to index the sorted alignment files
- `/path/to/*Aligned.sortedByCoord_sorted.out.bam` – the path to the sorted *Aligned.sortedByCoord_sorted.out.bam output files from the [step 4d](#4d-sort-aligned-reads), provided as a positional argument

**Input Data:**

- *Aligned.sortedByCoord_sorted.out.bam (sorted mapping to genome file, ourput from [Step 4d](#4d-sort-aligned-reads))

**Output Data:**

- *Aligned.sortedByCoord_sorted.out.bam.bai\# (index of sorted mapping to genome file)

<br>

---

## 5. Create Reference BED File

<br>

### 5a. Convert GTF to genePred File  

```bash
gtfToGenePred /path/to/annotation/gtf/file \
  /path/to/output/genePred/file

```

**Parameter Definitions:**

- `/path/to/annotation/gtf/file` – specifies the file(s) containing annotated reference transcripts in the standard gtf format, provided as a positional argument
- `/path/to/output/genePred/file` – specifies the location and name of the output genePred file(s), provided as a positional argument

**Input Data:**

- *.gtf ([genome annotation](../GeneLab_Reference_and_Annotation_Files/GL-DPPD-7101-E_ensembl_refs.csv))

**Output Data:**

- *.genePred (genome annotation in genePred format)

<br>

### 5b. Convert genePred to BED File  

```bash
genePredToBed /path/to/annotation/genePred/file \
  /path/to/output/BED/file
```

**Parameter Definitions:**

- `/path/to/annotation/genePred/file` – specifies the file(s) containing annotated reference transcripts in the genePred format, provided as a positional argument
- `/path/to/output/BED/file` – specifies the location and name of the output BED file(s), provided as a positional argument

**Input Data:**

- *.genePred (genome annotation in genePred format, output from [Step 5a](#5a-convert-gtf-to-genepred-file))

**Output Data:**

- *.bed (genome annotation in BED format)

<br>

---

## 6. Assess Strandedness, GeneBody Coverage, Inner Distance, and Read Distribution with RSeQC

<br>

### 6a. Determine Read Strandedness

```bash
infer_experiment.py -r /path/to/annotation/BED/file \
 -i /path/to/*Aligned.sortedByCoord_sorted.out.bam \
 -s 15000000 > /path/to/*infer_expt.out
```

**Parameter Definitions:**

- `-r` – specifies the path to the reference annotation BED file
- `-i` - specifies the path to the input bam file(s)
- `-s` - specifies the number of reads to be sampled from the input bam file(s), 15M reads are sampled
- `>` - redirects standard output to specified file
- `/path/to/*infer_expt.out` - specifies the location and name of the file containing the infer_experiment standard output

**Input Data:**

- *.bed (genome annotation in BED format, output from [Step 5b](#5b-convert-genepred-to-bed-file))
- *Aligned.sortedByCoord_sorted.out.bam (sorted mapping to genome file, output from [Step 4d](#4d-sort-aligned-reads))
- *Aligned.sortedByCoord_sorted.out.bam.bai (index of sorted mapping to genome file, output from [Step 4e](#4e-index-sorted-aligned-reads), although not indicated in the command, this file must be present in the same directory as the respective \*Aligned.sortedByCoord_sorted.out.bam file)

**Output Data:**

- *infer_expt.out (file containing the infer_experiment standard output)

<br>

### 6b. Compile Strandedness Reports

```bash
multiqc --interactive -n infer_exp_multiqc -o /path/to/infer_exp_multiqc/output/directory /path/to/*infer_expt.out/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*infer_expt.out/files` – the directory holding the *infer_expt.out output files from the [read strandedness step](#6a-determine-read-strandedness), provided as a positional argument

**Input Data:**

- *infer_expt.out (file containing the infer_experiment standard output, output from [Step 6a](#6a-determine-read-strandedness))

**Output Data:**

- infer_exp_multiqc.html\# (multiqc report)
- /infer_exp_multiqc_data\# (directory containing multiqc data)

<br>

### 6c. Evaluate GeneBody Coverage

```bash
geneBody_coverage.py -r /path/to/annotation/BED/file \
 -i /path/to/*Aligned.sortedByCoord_sorted.out.bam \
 -o /path/to/geneBody_coverage/output/directory/<sample_id>
```

**Parameter Definitions:**

- `-r` – specifies the path to the reference annotation BED file
- `-i` - specifies the path to the input bam file(s)
- `-o` - specifies the path to the output directory
- `/path/to/geneBody_coverage/output/directory/<sample_id>` - specifies the location and name of the directory containing the geneBody_coverage output files

**Input Data:**

- *.bed (genome annotation in BED format, output from [Step 5b](#5b-convert-genepred-to-bed-file))
- *Aligned.sortedByCoord_sorted.out.bam (sorted mapping to genome file, output from [Step 4d](#4d-sort-aligned-reads))
- *Aligned.sortedByCoord_sorted.out.bam.bai (index of sorted mapping to genome file, output from [Step 4e](#4e-index-sorted-aligned-reads), although not indicated in the command, this file must be present in the same directory as the respective \*Aligned.sortedByCoord_sorted.out.bam file)

**Output Data:**

- *.geneBodyCoverage.curves.pdf (genebody coverage line plot)
- *.geneBodyCoverage.r (R script that generates the genebody coverage line plot)
- *.geneBodyCoverage.txt (tab delimited file containing genebody coverage values used to generate the line plot)

<br>

### 6d. Compile GeneBody Coverage Reports

```bash
multiqc --interactive -n genebody_cov_multiqc -o /path/to/geneBody_coverage_multiqc/output/directory /path/to/geneBody_coverage/output/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/geneBody_coverage/output/files` – the directory holding the geneBody_coverage output files from [step 6c](#6c-evaluate-genebody-coverage), provided as a positional argument

**Input Data:**

- *.geneBodyCoverage.txt (tab delimited file containing genebody coverage values, output from [Step 6c](#6c-evaluate-genebody-coverage))

**Output Data:**

- geneBody_cov_multiqc.html\# (multiqc report)
- /geneBody_cov_multiqc_data\# (directory containing multiqc data)

<br>

### 6e. Determine Inner Distance (For Paired End Datasets ONLY)

```bash
inner_distance.py -r /path/to/annotation/BED/file \
 -i /path/to/*Aligned.sortedByCoord_sorted.out.bam \
 -k 15000000 \
 -l -150 \
 -u 350 \
 -o  /path/to/inner_distance/output/directory
```

**Parameter Definitions:**

- `-r` – specifies the path to the reference annotation BED file
- `-i` - specifies the path to the input bam file(s)
- `-k` - specifies the number of reads to be sampled from the input bam file(s), 15M reads are sampled
- `-l` - specifies the lower bound of inner distance (bp).
- `-u` - specifies the upper bound of inner distance (bp)
- `/path/to/inner_distance/output/directory` - specifies the location and name of the directory containing the inner_distance output files

**Input Data:**

- *.bed (genome annotation in BED format, output from [Step 5b](#5b-convert-genepred-to-bed-file))
- *Aligned.sortedByCoord_sorted.out.bam (sorted mapping to genome file, output from [Step 4d](#4d-sort-aligned-reads))
- *Aligned.sortedByCoord_sorted.out.bam.bai (index of sorted mapping to genome file, output from [Step 4e](#4e-index-sorted-aligned-reads), although not indicated in the command, this file must be present in the same directory as the respective \*Aligned.sortedByCoord_sorted.out.bam file)

**Output Data:**

- *.inner_distance.txt (log of read-wise inner distance results)
- *.inner_distance_freq.txt (tab delimited table of inner distances mapped to number of reads with that distance)
- *.inner_distance_plot.pdf (histogram plot of inner distance distribution)
- *.inner_distance_plot.r (R script that generates the histogram plot)

<br>

### 6f. Compile Inner Distance Reports

```bash
multiqc --interactive -n inner_dist_multiqc /path/to/inner_dist_multiqc/output/directory /path/to/inner_dist/output/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/inner_dist/output/files` – the directory holding the inner_distance output files from [Step 6e](#6e-determine-inner-distance-for-paired-end-datasets-only), provided as a positional argument

**Input Data:**

- *.inner_distance_freq.txt (tab delimited table of inner distances from [step 6e](#6e-determine-inner-distance-for-paired-end-datasets-only))

**Output Data:**

- inner_distance_multiqc.html\# (multiqc report)
- /inner_distance_multiqc_data\# (directory containing multiqc data)

<br>

### 6g. Assess Read Distribution

```bash
read_distribution.py -r /path/to/annotation/BED/file \
 -i /path/to/*Aligned.sortedByCoord_sorted.out.bam > /path/to/*read_dist.out
```

**Parameter Definitions:**

- `-r` – specifies the path to the reference annotation BED file
- `-i` - specifies the path to the input bam file(s)
- `>` - redirects standard output to specified file
- `/path/to/*read_dist.out` - specifies the location and name of the file containing the read_distribution standard output

**Input Data:**

- *.bed (genome annotation in BED format, output from [Step 5b](#5b-convert-genepred-to-bed-file))
- *Aligned.sortedByCoord_sorted.out.bam (sorted mapping to genome file, output from [Step 4d](#4d-sort-aligned-reads))
- *Aligned.sortedByCoord_sorted.out.bam.bai (index of sorted mapping to genome file, output from [Step 4e](#4e-index-sorted-aligned-reads), although not indicated in the command, this file must be present in the same directory as the respective \*Aligned.sortedByCoord_sorted.out.bam file)

**Output Data:**

- *read_dist.out (file containing the read distribution standard output)

<br>

### 6h. Compile Read Distribution Reports

```bash
multiqc --interactive -n read_dist_multiqc -o /path/to/read_dist_multiqc/output/directory /path/to/*read_dist.out/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*read_dist.out/files` – the directory holding the *read_dist.out output files from [Step 6g](#6g-assess-read-distribution) provided as a positional argument

**Input Data:**

- *read_dist.out (files containing the read_distributation standard output, output from [Step 6g](#6g-assess-read-distribution))

**Output Data:**

- read_dist_multiqc.html\# (multiqc report)
- /read_dist_multiqc_data\# (directory containing multiqc data)

<br>

---

## 7. Build RSEM Reference

```bash
rsem-prepare-reference --gtf /path/to/annotation/gtf/file \
 /path/to/genome/fasta/file \
 /path/to/RSEM/genome/directory/RSEM_ref_prefix

```

**Parameter Definitions:**

- `--gtf` – specifies the file(s) containing annotated transcripts in the standard gtf format
- `/path/to/genome/fasta/file` – specifies one or more fasta file(s) containing the genome reference sequences, provided as a positional argument
- `/path/to/RSEM/genome/directory/RSEM_ref_prefix` - specifies the path to the directory where the RSEM reference will be stored and the prefix desired for the RSEM reference files, provided as a positional argument

**Input Data:**

- *.fasta ([genome sequence](../GeneLab_Reference_and_Annotation_Files/GL-DPPD-7101-E_ensembl_refs.csv))
- *.gtf ([genome annotation](../GeneLab_Reference_and_Annotation_Files/GL-DPPD-7101-E_ensembl_refs.csv))

**Output Data:**

RSEM genome reference, which consists of the following files:

- RSEM_ref_prefix.chrlist
- RSEM_ref_prefix.grp
- RSEM_ref_prefix.idx.fa
- RSEM_ref_prefix.n2g.idx.fa
- RSEM_ref_prefix.seq
- RSEM_ref_prefix.ti
- RSEM_ref_prefix.transcripts.fa

<br>

---

## 8. Quantitate Aligned Reads

<br>

### 8a. Count Aligned Reads with RSEM

```bash
rsem-calculate-expression --num-threads NumberOfThreads \
 --alignments \
 --bam \
 --paired-end \
 --seed 12345 \
 --seed-length 20 \
 --estimate-rspd \
 --no-bam-output \
 --strandedness reverse|forward|none \
 /path/to/*Aligned.toTranscriptome.out.bam \
 /path/to/RSEM/genome/directory/RSEM_ref_prefix \
 /path/to/RSEM/counts/output/directory/<sample_id>
```

**Parameter Definitions:**

- `--num-threads` – specifies the number of threads to use
- `--alignments` - indicates that the input file contains alignments in sam, bam, or cram format
- `--bam` - specifies that the input alignments are in bam format
- `--paired-end` - indicates that the input reads are paired-end reads; this option should be removed if the input reads are single-end
- `--seed` - the seed for the random number generators used in calculating posterior mean estimates and credibility intervals; must be a non-negative 32-bit integer
- `--seed-length 20` - instructs RSEM to ignore any aligned read if it or its mates' (for paired-end reads) length is less than 20bp
- `--estimate-rspd` - instructs RSEM to estimate the read start position distribution (rspd) from the data
- `--no-bam-output` - instructs RSEM not to output any bam file
- `--strandedness` - defines the strandedness of the RNAseq reads; the `reverse` option is used if read strandedness (output from [step 6](#6a-determine-read-strandedness)) is antisense, `forward` is used with sense strandedness, and `none` is used if strandedness is half sense half antisense
- `/path/to/*Aligned.toTranscriptome.out.bam` - specifies path to input bam files, provided as a positional argument
- `/path/to/RSEM/genome/directory/RSEM_ref_prefix` - specifies the path to the directory where the RSEM reference is stored and its prefix, provided as a positional argument
- `/path/to/RSEM/counts/output/directory` – specifies the path to and prefix for the output file names; for GeneLab the prefix is the sample id

**Input Data:**

- RSEM genome reference (output from [Step 7](#7-build-rsem-reference))
- *Aligned.toTranscriptome.out.bam (sorted mapping to transcriptome, output from [Step 4a](#4a-align-reads-to-reference-genome-with-star))

**Output Data:**

- *genes.results\# (counts per gene)
- *isoforms.results\# (counts per isoform)
- *stat (directory containing the following stats files)
  - *cnt
  - *model
  - *theta

<br>

### 8b. Compile RSEM Count Logs

```bash
multiqc --interactive -n RSEM_count_multiqc -o /path/to/RSEM_count_multiqc/output/directory /path/to/*stat/files
```

**Parameter Definitions:**

- `--interactive` - force reports to use interactive plots
- `-n` - prefix name for output files
- `-o` – the output directory to store results
- `/path/to/*stat/files` – the directories holding the *stat output files from the [RSEM Counts step](#8a-count-aligned-reads-with-rsem), provided as a positional argument

**Input Data:**

- *stat (directory containing the following stats files, output from [Step 8a](#8a-count-aligned-reads-with-rsem))
  - *cnt
  - *model
  - *theta

**Output Data:**

- RSEM_count_multiqc.html\# (multiqc report)
- /RSEM_count_multiqc_data\# (directory containing multiqc data)

<br>

### 8c. Calculate Total Number of Genes Expressed Per Sample in R

```R
library(tximport)
library(tidyverse)

work_dir="/path/to/working/directory/where/script/is/executed/from" ## Must contain samples.txt file
counts_dir="/path/to/directory/containing/RSEM/counts/files"

setwd(file.path(work_dir))

### Pull in sample names where the "samples.txt" file is a single column list of sample names ###
samples <- read.csv(Sys.glob(file.path(work_dir,"samples.txt")), header = FALSE, row.names = 1, stringsAsFactors = TRUE)

##### Import RSEM Gene Count Data
files <- list.files(file.path(counts_dir),pattern = ".genes.results", full.names = TRUE)

### reorder the genes.results files to match the ordering of the samples in the metadata file
files <- files[sapply(rownames(samples), function(x)grep(paste0(x,'.genes.results$'), files, value=FALSE))]

names(files) <- rownames(samples)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

##### Count the number of genes with non-zero counts for each sample 
rawCounts <- txi.rsem$counts
NumNonZeroGenes <- (as.matrix(colSums(rawCounts > 0), row.names = 1))
colnames(NumNonZeroGenes) <- c("Number of genes with non-zero counts")

##### Export the number of genes with non-zero counts for each sample
setwd(file.path(counts_dir))
write.csv(NumNonZeroGenes,file='NumNonZeroGenes.csv')

## print session info ##
print("Session Info below: ")
print("")
sessionInfo()
```

**Input Data:**

- samples.txt (A newline delimited list of sample IDs)
- *genes.results (RSEM counts per gene, output from [Step 8a](#8a-count-aligned-reads-with-rsem))

**Output Data:**

- NumNonZeroGenes.csv (A samplewise table of the number of genes expressed)

<br>

---

## 9. Normalize Read Counts, Perform Differential Gene Expression Analysis, and Add Gene Annotations in R

<br>

### 9a. For Datasets With ERCC Spike-In

```R
## Install R packages if not already installed

install.packages("tidyverse")
source("https://bioconductor.org/biocLite.R")
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("tximport")
BiocManager::install("DESeq2")
BiocManager::install("Risa")
BiocManager::install("STRINGdb")
BiocManager::install("PANTHER.db")


## Install annotation R packages if not already installed - only the annotation package for the organism that the data were derived from is required

BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("org.Rn.eg.db")
BiocManager::install("org.Dr.eg.db")
BiocManager::install("org.Dm.eg.db")
BiocManager::install("org.Ce.eg.db")
BiocManager::install("org.Sc.sgd.db")
BiocManager::install("org.At.tair.db")
BiocManager::install("org.EcK12.eg.db")
BiocManager::install("MeSH.Bsu.168.eg.db")


##### Set up your environment #####

## Import libraries (tximport, DESeq2, tidyverse, Risa)

library(tximport)
library(DESeq2)
library(tidyverse)
library(Risa)