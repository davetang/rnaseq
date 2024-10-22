#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Cwd 'abs_path';

my $usage = "Usage: $0 <FASTQ directory>\n";
my $fastq_dir = shift or die $usage;

my %ids = ();

opendir(DIR, $fastq_dir) or die "Could not open $fastq_dir $!\n";
while(my $in = readdir(DIR)){
   next unless $in =~ /\.fastq\.gz$/;
   my @id = split(/_/, $in);
   $ids{$id[0]} = 1;
}
closedir(DIR);

print join(",", "sample", "fastq_1", "fastq_2", "strandedness"), "\n";
foreach my $id (sort {$a cmp $b} keys %ids){
   print join(",", $id, "$fastq_dir/${id}_1.fastq.gz", "$fastq_dir/${id}_2.fastq.gz", "auto"), "\n";
}

__END__
