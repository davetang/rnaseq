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

Map with STAR using various parameters.

```bash
./map_star.sh
```

* STAR will map reads to up to 10 places by default; reads mapping to more than 10 places will become unmapped
* Use the parameter `--outFilterMultimapNmax` to map reads up to *n* times
* Use the parameter `--outSAMunmapped Within` to keep unmapped reads in the SAM output file

