#!/usr/bin/env perl

use strict;
use warnings;

my $fasta = "mer_counts_dumps.fa";

my %check = ();

open(IN, '<', $fasta) || die "Could not open $fasta: $!\n";
while(<IN>){
   chomp;
   my $def = $_;
   chomp(my $seq = <IN>);

   my $num = 0;
   if ($def =~ /^>(\d+)/){
      $num = $1;
   } else {
      die "Could not extract k-mer count\n";
   }

   if (exists $check{$num}){
      next;
   } else {
      $check{$num} = 1;
      print "$def\n$seq\n";
   }

}
close(IN);

exit(0);

