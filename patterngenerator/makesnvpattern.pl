#!/usr/bin/perl
## the script assumes bowtie path is included.

if(@ARGV<5){ print "usage: $0 bedfile genomefasta genome(bowtie)index outdir outprefix\n"; exit;}
my ($bedfile,$genomefasta,$genomeindex,$outdir,$outprefix)=@ARGV;
if(!-d $outdir){ system("mkdir -p $outdir"); }
my ($scriptdir)=$0=~/(.+\/)[^\/]+$/; 
system("${scriptdir}snvbed2patternfile.3.pl $bedfile $genomefasta > $outdir/$outprefix.txt");
system("${scriptdir}patternfile2fasta.pl $outdir/$outprefix.txt > $outdir/$outprefix.fasta");
system("bowtie -f -p 4 -v 0 -a $genomeindex $outdir/$outprefix.fasta > $outdir/$outprefix.bowtieout");
system("${scriptdir}parse.pattern.bowtie.pl $outdir/$outprefix.bowtieout > $outdir/$outprefix.ntm");
system("${scriptdir}filter_patters_based_on_ntm.pl $outdir/$outprefix.txt $outdir/$outprefix.ntm > $outdir/$outprefix.uniq.txt");
system("${scriptdir}patternconverter $outdir/$outprefix.uniq.txt $outdir/$outprefix.pt");

# test
print("Checking pattern file...\n");
system("${scriptdir}patterntestreader $outdir/$outprefix.pt $outdir/$outprefix.pt-txt");
system("sort $outdir/$outprefix.pt-txt > $outdir/$outprefix.pt-txt.sorted");
system("sort $outdir/$outprefix.uniq.txt > $outdir/$outprefix.uniq.txt.sorted");
$diff=`diff $outdir/$outprefix.pt-txt.sorted $outdir/$outprefix.uniq.txt.sorted`;
if($diff!~/\S/){
  print("pattern file looks good!\n");
  exit(0); 
}
else {
  print("pattern file does not look good.\n");
  exit(1);
}

