#!/usr/bin/perl
if(@ARGV<2){ print "usage: $0 patternfile ntmfile > new.patternfile\n"; exit; }
my $patternfile = shift @ARGV;
my $ntmfile = shift @ARGV;

open NTM,$ntmfile or die "Can't open NTM file $ntmfile\n";
while(<NTM>){
  chomp;
  my ($index,$refntm,$altntm) = split/\t/;
  if($refntm==1 && $altntm==0){ $hash{$index}=1; }  ## uniquemappers with no alt mapping.
}
close NTM;


open PAT, $patternfile or die "Can't open patternfile $patternfile\n";
while(<PAT>){
  chomp;
  my ($pattern,$index,$ref_vs_alt) = split/\t/;
  if(exists $hash{$index}) { print "$_\n"; }
}
close PAT;

