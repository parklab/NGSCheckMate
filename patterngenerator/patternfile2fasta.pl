#!/usr/bin/perl
$k=0; ## k is a unique identifier.
while(<>){
  chomp;
  my ($seq,$index,$alt_vs_ref)=split/\t/;
  print ">$index.$alt_vs_ref.$k\n$seq\n";
  $k++;
}
