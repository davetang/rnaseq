## README

![STAR alignment strategy](img/star_strategy.png)

Material on STAR adapted from https://github.com/hbctraining/Intro-to-rnaseq-hpc-O2 distributed under the terms of the [Attribution 4.0 International (CC BY 4.0 license)](https://creativecommons.org/licenses/by/4.0/).

## Testing data

Data from [The transcriptional landscape and mutational profile of lung adenocarcinoma](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3483540/) deposited accession number [ERP001058](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40419) and raw data deposited at the SRA under [ERP001058](https://www.ncbi.nlm.nih.gov/sra?term=ERP001058).

Download [NCBI SRA Toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) to download data.

```bash
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/sratoolkit.2.10.8-ubuntu64.tar.gz
tar -xzf sratoolkit.2.10.8-ubuntu64.tar.gz

# add bin to your PATH
# export PATH=$PATH:/where/you/downloaded/src/sratoolkit.2.10.8-ubuntu64/bin

# then run the config tool
vdb-config --interactive
```

Download data using `fasterq-dump`.

```bash
for acc in ERR164550 ERR164559 ERR164560 ERR164563 ERR164569 ERR164585 ERR164613; do
   echo $acc
   fasterq-dump -p --outdir fastq ${acc}
done
```

## Reference

Download reference files as recommended by the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf).

```bash
mkdir index && cd index
ensembl_version=101
wget -c -N http://ftp.ensembl.org/pub/release-${ensembl_version}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget -c -N ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.primary_assembly.annotation.gtf.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip gencode.v35.primary_assembly.annotation.gtf.gz

sed -i 's/^chr//' gencode.v35.primary_assembly.annotation.gtf
sed -i 's/^M/MT/' gencode.v35.primary_assembly.annotation.gtf
```

Generate index.

```bash
mkdir index_100

../../src/STAR \
   --runMode genomeGenerate \
   --genomeDir index_100 \
   --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
   --sjdbGTFfile gencode.v35.primary_assembly.annotation.gtf \
   --sjdbOverhang 100 \
   --runThreadN 8
```

## Map

Parameters as per [processing RNA-seq data for GATK4](https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels) mostly default settings.

```bash
r1=/data/dtang/data/sra/fastq/ERR164550_1.fastq
r2=/data/dtang/data/sra/fastq/ERR164550_2.fastq

/usr/bin/time -v ../../src/STAR \
   --genomeDir ../index/index_100/ \
   --runThreadN 16 \
   --readFilesIn ${r1} ${r2} \
   --outSAMtype BAM SortedByCoordinate \
   --twopassMode Basic \
   --outFileNamePrefix ERR164550.GATK4. 2> star_gatk.log
```

Parameters adapted from [PanMutsRx](https://github.com/m081429/PanMutsRx).

```bash
r1=/data/dtang/data/sra/fastq/ERR164550_1.fastq
r2=/data/dtang/data/sra/fastq/ERR164550_2.fastq

/usr/bin/time -v ../../src/STAR \
   --genomeDir ../index/index_100/ \
   --runThreadN 16 \
   --readFilesIn ${r1} ${r2} \
   --alignSJDBoverhangMin 10 \
   --alignMatesGapMax 200000 \
   --alignIntronMax 200000 \
   --limitBAMsortRAM 31532137230 \ 
   --outSAMstrandField intronMotif \
   --outSAMtype BAM Unsorted \
   --outSAMunmapped Within KeepPairs \
   --outFileNamePrefix ERR164550.PanMutsRx. 2> star_pan_muts_rx.log
```

## Statistics

Number of reads: r1 + r2.

```bash
bc -l<<<$(cat /data/dtang/data/sra/fastq/ERR164550_1.fastq | wc -l)/2
77494486.00000000000000000000
```

Statistics using `samtools flagstat`.

```bash
samtools flagstat -@16 ERR164550.GATK4.Aligned.sortedByCoord.out.bam 
67348754 + 0 in total (QC-passed reads + QC-failed reads)
5374214 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
67348754 + 0 mapped (100.00% : N/A)
61974540 + 0 paired in sequencing
30987270 + 0 read1
30987270 + 0 read2
61974540 + 0 properly paired (100.00% : N/A)
61974540 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

samtools flagstat -@16 ERR164550.PanMutsRx.Aligned.out.bam
81902706 + 0 in total (QC-passed reads + QC-failed reads)
4408220 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
66345666 + 0 mapped (81.01% : N/A)
77494486 + 0 paired in sequencing
38747243 + 0 read1
38747243 + 0 read2
61937446 + 0 properly paired (79.92% : N/A)
61937446 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

bc -l<<<77494486+4408220
81902706
```

Almost all unmapped reads are shorter than the minimum allowed mapped length.

```bash
samtools view -f 4 ERR164550.PanMutsRx.Aligned.out.bam | perl -nle 'if (/(uT:.*)/ ){ print $1 }' | sort | uniq -c
  29568 uT:A:0
15504490 uT:A:1
  22982 uT:A:3
```

