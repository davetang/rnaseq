# README

Install using Conda.

```console
conda env create
conda activate kallisto
kallisto version
# kallisto, version 0.50.0
```

Build index using [GENCODE](https://www.gencodegenes.org/human/) FASTA file.

```console
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz
kallisto index -t 8 -i gencode.v44.idx gencode.v44.transcripts.fa.gz
```

If you get a `Illegal instruction` error, try downloading the binary from the
[GitHub repo](https://github.com/pachterlab/kallisto/releases/tag/v0.50.0).

```console
wget https://github.com/pachterlab/kallisto/releases/download/v0.50.0/kallisto_linux-v0.50.0.tar.gz
tar -xzf kallisto_linux-v0.50.0.tar.gz
```

Quantification.

```console
kallisto quant \
   -i gencode.v44.idx \
   -t 8 \
   -o output \
   -b 100 <(zcat SRR22891573_1.fastq.gz) <(zcat SRR22891573_2.fastq.gz)
```

## Differential expression analysis

[Getting started with
sleuth](https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html).
