#!/usr/bin/perl

while(<>){
  chomp;
  my @line=split/\t/;
  $count{$line[0]}++;
}

for my $p (keys %count){
  my ($snv,$ref,$k) = split/\./,$p;
  $refmaxcount{$snv}=$count{$p} if ( $ref==0 && (!exists($refmaxcount{$snv}) || $count{$p} > $refmaxcount{$snv}) );
  $altmaxcount{$snv}=$count{$p} if ( $ref==1 && (!exists($altmaxcount{$snv}) || $count{$p} > $altmaxcount{$snv}) );
  $allsnv{$snv}=1;
}

print "snv\trefmaxcount\taltmaxcount\n";
for my $snv (sort {$a<=>$b} keys %allsnv){
  if(!exists $refmaxcount{$snv}) { $refmaxcount{$snv}=0; }
  if(!exists $altmaxcount{$snv}) { $altmaxcount{$snv}=0; }
  print "$snv\t$refmaxcount{$snv}\t$altmaxcount{$snv}\n";
}


