<p align="center">
  <img src="https://parklab.github.io/NGSCheckMate/logo.svg"
       alt="NGS CheckMate" />
</p>

# NGS CheckMate
Software program for checking sample matching for NGS data


## Usage

### Main program

#### 1) Input VCF

- VCF files in `vcf_dir`:

   ```bash
   python ncm.py -V -d vcf_dir -O output_Dir -N outputfileName -bed bed_file
   ```

- VCF lists in `bam_list_file`:

   ```bash
   python ncm.py -V -l vcf_list_file -O output_Dir -N outputfileName -bed bed_file
   ```

-

#### 2) Input BAM

Change your configruration.

Please open `ncm.py`, find `run_mpileup()` function, modifying below lists as your configurations:

```python
SAMTOOLS="/NAS/nas100-5/tools/samtools-0.1.19/samtools"

BCFTOOLS="/NAS/nas100-5/tools/samtools-0.1.19/bcftools/bcftools"

REF="/NAS/nas33-2/mpileup/hg19.fasta" or REF="/NAS/nas33-2/mpileup/GRCh37-lite.fa"
```

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


### Supporting scripts

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


## Acknoledgements

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
