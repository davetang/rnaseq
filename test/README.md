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

# for BWA, you need to create the index files in ../raw/chrX_data/genome/
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
* However, the `-h` parameter does not seem to always report the additional locations. For example, for the sequence "5"

```bash
cat eg.fa | grep -A 1 "^>5$"
>5
GAATTTTTGTATAAGGTGTAAGGAAGGGATCCAGTTTCAGCTTTCTACATATGGCTAGCCAGTTTTCCCAGCACC

# no XA tag
samtools view bwa_h10.bam | grep "^5\t"
5       16      chrX    123057893       0       75M     *       0       0       GGTGCTGGGAAAACTGGCTAGCCATATGTAGAAAGCTGAAACTGGATCCCTTCCTTACACCTTATACAAAAATTC     *       NM:i:0  MD:Z:75 AS:i:75 XS:i:75

# HISAT2 result
samtools view hisat2_default.bam | grep "^5\t" | sort -k4n
5       256     chrX    112915551       1       75M     *       0       0       GAATTTTTGTATAAGGTGTAAGGAAGGGATCCAGTTTCAGCTTTCTACATATGGCTAGCCAGTTTTCCCAGCACC     IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII  AS:i:0   ZS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:75 YT:Z:UU NH:i:5
5       0       chrX    122547852       1       75M     *       0       0       GAATTTTTGTATAAGGTGTAAGGAAGGGATCCAGTTTCAGCTTTCTACATATGGCTAGCCAGTTTTCCCAGCACC     IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII  AS:i:0   ZS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:75 YT:Z:UU NH:i:5
5       272     chrX    123057893       1       75M     *       0       0       GGTGCTGGGAAAACTGGCTAGCCATATGTAGAAAGCTGAAACTGGATCCCTTCCTTACACCTTATACAAAAATTC     IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII  AS:i:0   ZS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:75 YT:Z:UU NH:i:5
5       256     chrX    148326668       1       75M     *       0       0       GAATTTTTGTATAAGGTGTAAGGAAGGGATCCAGTTTCAGCTTTCTACATATGGCTAGCCAGTTTTCCCAGCACC     IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII  AS:i:0   ZS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:75 YT:Z:UU NH:i:5
5       256     chrX    148348072       1       75M     *       0       0       GAATTTTTGTATAAGGTGTAAGGAAGGGATCCAGTTTCAGCTTTCTACATATGGCTAGCCAGTTTTCCCAGCACC     IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII  AS:i:0   ZS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:75 YT:Z:UU NH:i:5

# STAR result
cat multimap_20.Aligned.out.sam | grep "^5\t" | sort -k4n
5       272     chrX    112098806       0       74M1S   *       0       0       GGTGCTGGGAAAACTGGCTAGCCATATGTAGAAAGCTGAAACTGGATCCCTTCCTTACACCTTATACAAAAATTC     *       NH:i:9  HI:i:9  AS:i:72 nM:i:0
5       256     chrX    112915551       0       75M     *       0       0       GAATTTTTGTATAAGGTGTAAGGAAGGGATCCAGTTTCAGCTTTCTACATATGGCTAGCCAGTTTTCCCAGCACC     *       NH:i:9  HI:i:6  AS:i:73 nM:i:0
5       256     chrX    112924065       0       1S74M   *       0       0       GAATTTTTGTATAAGGTGTAAGGAAGGGATCCAGTTTCAGCTTTCTACATATGGCTAGCCAGTTTTCCCAGCACC     *       NH:i:9  HI:i:7  AS:i:72 nM:i:0
5       256     chrX    122547852       0       75M     *       0       0       GAATTTTTGTATAAGGTGTAAGGAAGGGATCCAGTTTCAGCTTTCTACATATGGCTAGCCAGTTTTCCCAGCACC     *       NH:i:9  HI:i:5  AS:i:73 nM:i:0
5       272     chrX    122919777       0       74M1S   *       0       0       GGTGCTGGGAAAACTGGCTAGCCATATGTAGAAAGCTGAAACTGGATCCCTTCCTTACACCTTATACAAAAATTC     *       NH:i:9  HI:i:4  AS:i:72 nM:i:0
5       272     chrX    123057893       0       75M     *       0       0       GGTGCTGGGAAAACTGGCTAGCCATATGTAGAAAGCTGAAACTGGATCCCTTCCTTACACCTTATACAAAAATTC     *       NH:i:9  HI:i:3  AS:i:73 nM:i:0
5       272     chrX    126951647       0       74M1S   *       0       0       GGTGCTGGGAAAACTGGCTAGCCATATGTAGAAAGCTGAAACTGGATCCCTTCCTTACACCTTATACAAAAATTC     *       NH:i:9  HI:i:8  AS:i:72 nM:i:0
5       0       chrX    148326668       0       75M     *       0       0       GAATTTTTGTATAAGGTGTAAGGAAGGGATCCAGTTTCAGCTTTCTACATATGGCTAGCCAGTTTTCCCAGCACC     *       NH:i:9  HI:i:1  AS:i:73 nM:i:0
5       256     chrX    148348072       0       75M     *       0       0       GAATTTTTGTATAAGGTGTAAGGAAGGGATCCAGTTTCAGCTTTCTACATATGGCTAGCCAGTTTTCCCAGCACC     *       NH:i:9  HI:i:2  AS:i:73 nM:i:0
```

