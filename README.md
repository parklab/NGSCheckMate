# [NGSCheckMate](https://github.com/parklab/NGSCheckMate/)
NGSCheckMate is a software package for identifying next generation sequencing (NGS) data files from the same individual. It analyzes various types of NGS data files including (but not limited to) whole genome/exome sequencing, RNA-seq, ChIP-seq, targeted sequencing of various sequencing depths. It also supports checking sample pairing across data types. For example, our performance evaluation showed high accuracy of NGSCheckMate in pairing DNA sequencing files with RNA-sequencing files, or RNA-seq files and ChIP-seq files from the same individual. NGSCheckMate takes BAM (reads with alignment), VCF (variants) or FASTQ (unaligned reads) files as input. The alignment-free module of NGSCheckMate was specifically designed for FASTQ input. Since it does not perform alignment or variant calling step, it runs very fast expedited by read subsampling and parallelization (e.g., For RNA-seq data, it takes less than one minute using a single core). The idea behind how NGSCheckMate works is to evaluate similarity (correlation) of allele fractions of known single-nucleotide polymorphisms (SNPs) among input files. Currently, it works only for human data, but it can be easily extended to non-human data by providing a bed file including known SNPs for the relevant organism. 

<p align="center">
  <img src="https://parklab.github.io/NGSCheckMate/logo.svg"
       alt="NGS CheckMate" />
</p>

## Table of contents

* [Quick start](#Quick-start)
* [Usage](#Usage)
* [Supporting scripts](#Supporting-scripts)
* [Authors](#Authors)
* [Acknoledgements](#Acknoledgements)


## Quick start

### Download NGSCheckMate
* [Download the latest release](https://github.com/parklab/NGSCheckMate/)
* Clone the repo : `git clone https://github.com/parklab/NGSCheckMate.git`

### Requirements
#### Software environment
* Unix/Linux System
* Python 2.6 or above
* Samtools version 0.1.19 (required only for BAM input)
* Bcftools version 0.1.19 (required only for BAM input)
* R 3.1 or above (required only for image output of clustering dendrogram)

#### Additional files
* Required for BAM
```
Human reference genome (hg19 or GRCh37 fasta file)
A bed file(.bed) that lists the locations of selected SNPs (including in the package)
```
* Required for VCF - Output of samtools mpileup
```
A bed file(.bed) that lists the locations of selected SNPs (including in the package)
```
* Required for FASTQ input
```
A binary pattern file(.pt) that lists the k-mer sequences spanning selected SNPs (included in the package)
```

## Usage

#### 1) Input VCF
```
Usage: python ncm.py -V <–d INPUT_DIR | -I INPUT_LIST_FILE>  -bed BED_FILE –O OUTPUT_DIR [options]
```

* Required arguments
```
	-d DIR	directory that contains input files
	  or
	-I FILE	text file that lists input files (one absolute path per line) 
	
	-bed FILE	bed file that lists coordinates of known SNPs (included in the package) 
	(use SNP_hg19.bed for hg19, or SNP_GRCh37.bed for GRCh37 in feature folder)
	-O PATH		The name of the output directory
```
* Optional arguments
```
	-N NAME   The name of the output file (default: “output”)
	-f 		Use strict correlation threshold. Recommended if your data contains family members.
	-t file	A file with test sample files 
```

For examples,
- VCF files in `vcf_dir`:

   ```bash
   python ncm.py -V -d vcf_dir -O output_Dir -N outputfileName -bed bed_file
   ```

- VCF lists in `vcf_list_file`:

   ```bash
   python ncm.py -V -l vcf_list_file -O output_Dir -N outputfileName -bed bed_file
   ```

-

#### 2) Input BAM

Set paths in the configuration file (required only for BAM input)

Edit ncm.conf file in the downloaded package directory according to your environment. 

```python
REF=”absolute path”  # path for the reference fasta file
SAMTOOLS=”absolute path” # path for SAMTOOLS 
BCFTOOLS=”absolute path” # path for BCFTOOLS
```

```
Usage: python ncm.py -B <–d INPUT_DIR | -I INPUT_LIST_FILE>  -bed BED_FILE –O OUTPUT_DIR [options]
```

* Required arguments
```
	-d DIR	directory that contains input files
	  or
	-I FILE	text file that lists input files (one absolute path per line) 
	
	-bed FILE	bed file that lists coordinates of known SNPs (included in the package) 
	(use SNP_hg19.bed for hg19, or SNP_GRCh37.bed for GRCh37 in feature folder)
	-O PATH		The name of the output directory
```
* Optional arguments
```
	-N NAME   The name of the output file (default: “output”)
	-f 		Use strict correlation threshold. Recommended if your data contains family members.
	-t file	A file with test sample files 
```

For examples,
 - BAM files in `bam_dir`:

   ```bash
   python ncm.py -B -d bam_dir -O output_Dir -N outputfileName -bed bed_file
   ```

 - Bam files in `bam_list_file`:

    ```bash
    python ncm.py -B -d bam_list_file -O output_Dir -N outputfileName -bed bed_file
    ```

-

#### 3) Input FASTQ

```bash
Usage : ./ngscheckmate_fastq <options> -1 fastqfile1 [-2 fastqfile2]  patternfile(.pt) > output.vaf

	Input arguments (required)
	  patternfile : a binary file containing sequences flanking representative snv sites, along with markers indicating the snv index and whether the sequence represents reference or alternative allele.
	  fastqfile1 : see below 'Options'.

	Options
	  -1, --fastq1 <fastq_file_1> : fastq file for SE data or first fastq file for a PE data. (required)
	  -2, --fastq2 <fastq_file_2> : second fastq file for a PE data.
	  -s, --ss <subsampling_rate> : subsampling rate (default 1.0)
	  -d, --depth <desired_depth> : as an alternative to a user-defined subsampling rate, let the program compute the subsampling rate given a user-defined desired_depth and the data.
	  -R, --reference_length <reference_length> : The reference length (default : 3E9) to be used for computing subsampling rate. If the data is NOT WGS from human, and if you're using the -d option, it is highly recommended to specify the reference length. For instance, if your data is human RNA-seq, the total reference length could be about 3% of the human genome, which can be set as 1E8.
	  -L, --pattern_length <pattern_length> : The length of the flanking sequences being used to identify SNV sites. Default is 21bp. It is recommended not to change this value, unless you have created your own pattern file with a different pattern length.

	  -p, --maxthread <number_of_threads> : number of threads to use (default : 1 )
	  -j, --nodeptherror : in case estimated subsampling rate is larger than 1, do not stop but reset it to 1 and continue.
```
```bash
Usage: python vaf_ncm.py -f -I <input_directory> -O <output_directory> -N output

       -I : input directory that contains the output (vaf) files of ngscheckmate_fastq.
       -O : output directory
       -N : output_filename_tag
```


## Supporting-scripts

#### 1) Patterngenerator

This set of scripts generates the .pt file used by the fastq module, given a bed file containing a set of SNP positions. It assumes a file containing a whole genome sequence and the bowtie alignment program.


#### 2) Graph generator (Rscript)

This script with a set of xgmml templates is used for generating a graph representing matching files as connected nodes. The output format is in .xgmml, which can be read by Cytoscape.


```R
source("graph/ngscheckmate2xgmml.R")
create.xgmml.from.ngscheckmateout(label.file,ngscheckmateoutput.file,output.xgmml)
```

 - Label file : a tab-delimited text file containing a bam file name (1st column), an individual identifier (2nd column) and optionally, a file identifier (3rd column) for each line. An individual identifier must be unique to a subject (e.g. both tumor and normal samples from the same individual must have the same individual identifier). A file identifier must be unique to a file name.
 - ngscheckmateoutput.file : the output text file of NGSCheckMate. It is a tab-delimited text file containing two bam file names (1st and 2nd columns), VAF correlation (3rd column) and average depth (4th column). It may contain either all pairs or matched pairs, depending on the option used to run NGSCheckMate. Either type works.
 - Sample label file (sample.label.txt) and ngscheckmateouput.file (sample.input.txt) can be found in the subdirectory graph/.



## Authors

Software programs : Alice Lee, [Sejoon Lee][sejooning] & [Soo Lee][SooLee]


## Acknowledgments

The logos (![alt text][ncmLogo] & ![alt text][ncmIcon]) have been created by [Fritz Lekschas][flekschas]. They are composed of the following great icons:
 - [DNA][iconDna] created by Irene Hoffman ([CC BY 3.0 US][cc])
 - [King][iconKing] created by Yuri Mazursky ([CC BY 3.0 US][cc])
 - [Queen][iconQueen] created by Yuri Mazursky ([CC BY 3.0 US][cc])

[sejooning]: https://github.com/sejooning
[SooLee]: https://github.com/SooLee
[cc]: https://creativecommons.org/licenses/by/3.0/us/
[flekschas]: https://github.com/flekschas
[iconDna]: https://thenounproject.com/term/dna/57369/
[iconKing]: https://thenounproject.com/term/king/224748/
[iconQueen]: https://thenounproject.com/term/queen/224753/
[ncmLogo]: https://parklab.github.io/NGSCheckMate/logo-16px.png "NGS CheckMate Logo"
[ncmIcon]: https://parklab.github.io/NGSCheckMate/icon-16px.png "NGS CheckMate Icon"
