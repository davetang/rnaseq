## README

A workflow to test the multimapping behaviour of HISAT2 and STAR (and BWA). We will find repeated sequence on chromosome X using `jellyfish` and map them back using various aligners.

Create the Conda environment if you haven't already and activate it before running the following commands.

```bash
conda activate jellyfish

# -m Length of mer
# -s Initial hash size
# -t Number of threads
# -C Count both strand, canonical representation (false)
jellyfish count -m 75 -s 100M -t 8 -C ../raw/chrX_data/genome/chrX.fa
jellyfish dump mer_counts.jf > mer_counts_dumps.fa
```

Create a FASTA file with one representative sequence for k-mer number; the FASTA definition line indicates the number of times the k-mer was found on chromosome X.

```bash
./create_fasta.pl > eg.fa

head eg.fa
>1
AATGTATTCTTAACCATATGATCCAGCAATTGTATCTCTTAGTACTTACCAAAATTAGTTGTTGAACTTGTTCCC
>2
ACCCTCAGCTTTAGGTGGCTCAGAACAGAGACAGAGAGAGAGAGAGAGACTCTGTATGTTTGGGAGAAAGTAAGG
>4
ACAGAGTTTCCTCATGTTGGCCAGCCTGGTCTCAAACTCTTGACTTCAAGTGATCCGCCTGCCTGGGCTTCCCAG
>3
AGCTTAGTTTGGCTGGATATGAAATTCTGAGTTGAAAATTCTTTTCTTTAAGAATGTTGAATATTGGCCCCCACT
>36
AGGATTGTGTTGGAAAAGGAAATATCTTCTCCTAAAAACGACATAGAAGCATTCTCAGAAACTGCTCTGTGATGA

cat eg.fa | grep "^>" | wc -l
1224

# remove jellyfish files
rm mer_counts.jf mer_counts_dumps.fa
```

Map with STAR and HISAT2 (and BWA) using various parameters.

```bash
./map_star.sh
./map_hisat2.sh
./map_bwa.sh
```

## Summary

STAR:

* STAR will map reads to up to 10 places by default; reads mapping to more than 10 places will become unmapped
* Use the parameter `--outFilterMultimapNmax` to map reads up to *n* times
* Use the parameter `--outSAMunmapped Within` to keep unmapped reads in the SAM output file

HISAT2:

* HISAT2 will map reads to up to 10 places by default; reads mapping to more than 10 places will be mapped by only 10 locations are listed
* Use the parameter `-k` to map reads up to *k* times; reads mapping to more than *k* times will be mapped by only *k* locations are listed
* When reporting reads mapping to more than *k* times, reported locations seem to be deterministic, which suggests that the first *k* locations are reported

BWA:

* BWA MEM will map all reads and the multimapping status is reflected in the mapping quality (60 for uniquely mapping and 0 when mapped to 2 or more locations)
* BWA MEM will report the other mapping locations using the `XA` tag but will only report up to three alternate locations by default
* Use the parameter `-h` to report more mapping locations

