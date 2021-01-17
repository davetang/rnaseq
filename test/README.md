## README

Find repeated sequence using `jellyfish`.

```bash
jellyfish count -m 75 -s 100M -t 8 -C ../raw/chrX_data/genome/chrX.fa
jellyfish dump mer_counts.jf > mer_counts_dumps.fa
```

Create a FASTA file with one representative sequence for k-mer number.

```bash
./create_fasta.pl > eg.fa

# remove jellyfish files
rm mer_counts.jf mer_counts_dumps.fa
```

Map with STAR and HISAT2 using various parameters.

```bash
./map_star.sh
map_hisat2.sh
```

## Summary

STAR:

* STAR will map reads to up to 10 places by default; reads mapping to more than 10 places will become unmapped
* Use the parameter `--outFilterMultimapNmax` to map reads up to *n* times
* Use the parameter `--outSAMunmapped Within` to keep unmapped reads in the SAM output file

HISAT2:

* HISAT2 will map reads to up to 10 places by default; reads mapping to more than 10 places will be mapped by only 10 locations are listed
* Use the parameter `-k` to map reads up to *k* times; reads mapping to more than *k* times will be mapped by only *k* locations are listed

