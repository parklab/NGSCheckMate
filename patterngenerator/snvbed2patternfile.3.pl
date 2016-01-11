#!/usr/bin/perl
#Author : Soo Lee
#This version uses the correct SNV coordinate. (END instead of START). (version 3)
#This version adds reverse complemtary sequences. (version 2)

if(@ARGV<2) {print "usage: $0 bedfile fastafile > output.fa\n"; exit;}

our ($CHR,$START,$END)=(0,1,2);
our %RevCompBasePairs=qw/A T T A G C C G a t t a g c c g U A u a R Y r y Y R y r M K m k K M k m S S s s W W w w H D D H h d d h B V V B b v v b N N n n/;

my $bedfile=$ARGV[0];
my $fastafile = $ARGV[1];  ## Genome

my @offset=(10,20,0); # starting point of the extended region, relative to the snv site. (w/o minus)
my $len=21; # length of the extended region.

my %fasta=&readfastafile2hash($fastafile,1); # use the first part of the header before space

my $index=0;
$"="\t";
open IN,$bedfile or die "can't open bedfile\n";
<IN>; #discard header
while(<IN>){ 
 chomp;
 my $i=[(split/\t/)];
 my $chr=$i->[$CHR]; 
 my @seq=();
 for my $j (0..2){
   $seq[$j]=uc substr($fasta{$chr},$i->[$END]-1-$offset[$j],$len); 
   #print STDERR "$seq[$j]\n"; ##DEBUGGING
 }

 my $refbase = substr($seq[0],$offset[0],1);
 if($refbase !~ /[ACGT]/) { print STDERR "Warning: the reference base is '$refbase'. ($_)\n"; }
 if($refbase ne substr($seq[1],$offset[1],1)) { die "seq1 and seq2 not coherent: $seq[0], $seq[1]. ($_)\n"; }
 if($refbase ne substr($seq[2],$offset[2],1)) { die "seq1 and seq3 not coherent: $seq[0], $seq[2]. ($_)\n"; }
   
 my @altseqs=();
 for my $j (0..2){
  for my $b ('A','C','G','T'){
    next if($b eq $refbase);
    my $altseq=$seq[$j];
    substr($altseq,$offset[$j],1)=$b;
    push @altseqs,$altseq;
  }
 }

 for my $j (0..2){
   print "$seq[$j]\t$index\t0\n";
   printf "%s\t$index\t0\n",&RevComp($seq[$j]);
 }
 for my $j (0..$#altseqs){
   print "$altseqs[$j]\t$index\t1\n";
   printf "%s\t$index\t1\n",&RevComp($altseqs[$j]);
 }

 $index++;
}
close IN;





#################
#reads fastafile and save into hash where header is a key and sequence is the value.
#only single-line fasta file.
sub readfastafile2hash {
 my $fastafile = shift @_;
 my $option = shift @_;  ## either 0 or 1 (1 means use only the first part of the header before space)

 my %hash=();
 open $FA,$fastafile or die "can't open fasta file\n";
 while(my $read = <$FA>){
   chomp $read;
   if($option==0) { if($read=~/^>/) { my $header = $'; chomp($seq=<$FA>); $hash{$header}="$seq\n"; }}
   elsif($option==1) { if($read=~/^>/) { my ($header) = split/\s/,$'; chomp($seq=<$FA>); $hash{$header}="$seq\n"; }}
   else { die "wrong option\n"; }
 }
 close $FA;

 return %hash;
}

# RevComp(s) returns a reverse complement DNA sequence of a non-degenerate or degenerate DNA/sequences.
sub RevComp {
  my $s=shift @_;
  my $new_s='';
  for my $b (split//,$s) {
   if(!exists $RevCompBasePairs{$b}) { $new_s=$b.$new_s; }
   else { $new_s=$RevCompBasePairs{$b}.$new_s; }
  }
  return $new_s;
}



