# NGSCheckMate
Software program for checking sample matching for NGS data

-
1) input VCF
- VCF files in vcf_dir

python ncm.py -V -d vcf_dir -O output_Dir -N outputfileName -bed bed_file

- VCF lists in bam_list_file

python ncm.py -V -l vcf_list_file -O output_Dir -N outputfileName -bed bed_file


-
2) input BAM

change your configruration

please open ncm.py, find run_mpileup() function, modifying below lists as your configurations

SAMTOOLS="/NAS/nas100-5/tools/samtools-0.1.19/samtools"

BCFTOOLS="/NAS/nas100-5/tools/samtools-0.1.19/bcftools/bcftools"

REF="/NAS/nas33-2/mpileup/hg19.fasta" or REF="/NAS/nas33-2/mpileup/GRCh37-lite.fa"

-BAM files in bam_dir
 python ncm.py -B -d bam_dir -O output_Dir -N outputfileName -bed bed_file
-Bam files in bam_list_file
 python ncm.py -B -d bam_list_file -O output_Dir -N outputfileName -bed bed_file
 
3) input FASTQ

