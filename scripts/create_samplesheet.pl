#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use Cwd 'abs_path';

my $script_dir = dirname(__FILE__);
my $sample_dir = abs_path("$script_dir/../raw/chrX_data/samples");

my %ids = ();

opendir(DIR, $sample_dir) or die "Could not open $sample_dir: $!\n";
while(my $in = readdir(DIR)){
   next unless $in =~ /\.fastq\.gz$/;
   my @id = split(/_/, $in);
   $ids{$id[0]} = 1;
}
closedir(DIR);

print join(",", "sample", "fastq_1", "fastq_2", "strandedness"), "\n";
foreach my $id (sort {$a cmp $b} keys %ids){
   print join(",", $id, "$sample_dir/${id}_chrX_1.fastq.gz", "$sample_dir/${id}_chrX_2.fastq.gz", "auto"), "\n";
}

__END__
