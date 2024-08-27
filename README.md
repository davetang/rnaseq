## Table of Contents

  - [Basics](#basics)
  - [Core concepts](#core-concepts)
    - [cDNA libraries](#cdna-libraries)
  - [Count data](#count-data)
    - [The challenges of count data](#the-challenges-of-count-data)
    - [Dispersion](#dispersion)
  - [RNA-seq pipelines](#rna-seq-pipelines)
  - [Test data](#test-data)
  - [Processing](#processing)
  - [Comparing quantifications](#comparing-quantifications)
  - [nf-core/rnaseq](#nf-corernaseq)
    - [Download reference transcriptome](#download-reference-transcriptome)
    - [Running nf-core/rnaseq on the test data](#running-nf-corernaseq-on-the-test-data)
    - [nf-core/rnaseq results](#nf-corernaseq-results)
    - [Creating a reference](#creating-a-reference)
  - [Papers](#papers)

## Basics

RNA molecules in a population of cells (or homogenised tissue) are reverse complemented into complementary DNA (cDNA) and sequenced on a high-throughput sequencer. There are two main classes of data output:

1. The output of the cDNA sequence
2. The abundances of the different cDNA sequences

When most people are talking about RNA sequencing (RNA-seq), they are usually referring to the study of the abundances and how they differ under different conditions. More specifically, the main goal is to quantify systematic changes from different conditions and to assess the statistical significance of the differences. Systematic changes need to be distinguished from sampling and technical variances.

Another important point is that it is not possible to sequence and count **all** RNA molecules in a sample because the protocols are not 100% efficient; RNA and their intermediates get lost during library preparation. Instead a statistical sample is produced, in the same way a census is a sample of the population. The amount to sample, i.e., how deep should we sequence, depends on the purpose of the study and also on the complexity of the biological sample, i.e., how many different species of RNA exist. The RNA-seq protocol to use is also dependent on the purpose of the study, since some protocols can exclude the class of RNA that you may be interested in studying.

## Core concepts

* A **sequencing library** is the collection of cDNA molecules used as input for the sequencing machine
* **Fragments** are the molecules that are sequenced. Since most widely used sequencers can only deal with molecules of length 300-1000 nucleotides, longer cDNA molecules are fragmented into this size range.
* A **sequencing read** is the sequence of a fragment. Most times reads come in pairs, which is known as paired-end sequencing, and sometimes it is possible for both reads to completely cover the fragment, i.e., overlapping reads.

Reads are typically aggregated together in a transcript/gene manner; reads belonging to the same transcript/gene are grouped together. In RNA-seq, reads are usually mapped to a known reference but if a reference does not exist, reads can be clustered together based on their sequence similarity.

### cDNA libraries

A cDNA library should contain all representative sequences of an mRNA population. Furthermore, the cDNA sequences should be full-length copies of the original mRNA. Therefore the construction of a high quality cDNA library is essential for RNA-seq. A cDNA clone (one particular copy in the library) represents the fully processed RNA sequence generated from the genomic sequence.

There are many procedures for synthesising cDNA to create cDNA libraries and many focus on maximising the amount of cDNA produced from limited starting amounts of mRNA. The completeness of the cDNA synthesis is variable and unpredictable. This variation can be introduced during the following steps in the construction protocol:

* mRNA isolation
* First-strand cDNA synthesis, and
* Second-strand cDNA synthesis
* RNase contamination

The first-strand cDNA synthesis step relies on the reverse transcription of RNA and the reverse transcriptase used can create fluctuations in the quantity and quality of the first-strand cDNA.

There are many different protocols for second-strand synthesis, such as hairpin-primed synthesis and the Gubler and Hoffman procedure.

Once second-strand synthesis is complete, double-stranded cDNA can be cloned and later sequenced as a cDNA library.

## Count data

The count table tallies the number of reads mapped to genes/transcripts, where each row corresponds to a gene and each column to a particular sample. This table contains integer values and the value in the $i$ th row and $j$ th column indicates how many reads have been mapped to gene/transcript $i$ in sample $j$. These are raw counts of the sequencing reads and not some derived quantity, such as normalised counts; it is essential that the values are raw counts or else the statistical models typically used in RNA-seq analysis are not valid.

### The challenges of count data

* Count data have a large dynamic range, which starts from zero and can go up to millions. The variance and the distribution shape of the data in different parts of the dynamic range are very different, i.e., heteroscedasticity.
* The data are non-negative integers and their distribution is not symmetric, therefore normal or log-normal distribution models may be a poor fit.
* We need to understand the systematic sampling biases and adjust for them. For example adjusting for different sequencing depths.
* The estimation of dispersion parameters is difficult with the small sample sizes typically seen in RNA-seq studies.

### Dispersion

Consider a sequencing library that contains $n_1$ fragments corresponding to gene 1, $n_2$ fragments for gene 2, and so on, with a total library size of $n = n_1 + n_2 + \ldots$. This library is then sequenced and the identity of $r$ randomly sampled fragments are determined. To paint a better mental picture, the following are some typical numbers. The number of genes will be in the order of tens of thousands; the value of $n$ depends on the amount of cells that were used to prepare the library and typically this is in the order of billions or trillions; and the number of reads $r$ is usually in the tens of millions, which is much smaller than $n$. Sequencing is sampling from $n$.

From this we can conclude that the probability that a given read maps to the $i$ th gene is $p_i = n_i/n$ (ratio of a specific fragment to all fragments) and this is independent of the outcomes for all the other reads. So we can model the number of reads for gene $i$ by a Poisson distribution, where the _rate_ of the Poisson process is the product of $p_i$, the initial proportion of fragments for the $i$ th gene, times $r$, the number of reads sequenced; that is $\lambda_i = rp_i$.

In practice, we are usually not interested in modeling the read counts within a single library, but in comparing the counts between libraries. That is, we want to know whether any differences that we see between different biological conditions are larger than what we might expect even between biological replicates. Empirically, it turns out that replicate experiments vary more than the Poisson distribution predicts.

Intuitively, what happens is that $p_i$, and therefore also $\lambda_i$, varies even between biological replicates. To account for that variation, we need to add another layer of modeling on top and it turns out that the gamma-Poisson (a.k.a. negative binomial) distribution suits our modeling requirements. Instead of a single $\lambda$, which represents both mean and variance, this distribution has two parameters. In principle, these can be different for each gene and we can estimate them from the data.

## RNA-seq pipelines

* [HISAT2 + StringTie2 pipeline](https://davetang.org/muse/2017/10/25/getting-started-hisat-stringtie-ballgown/)
* [STAR, Cufflinks, RSEM](https://pubmed.ncbi.nlm.nih.gov/27662878/)
* [Kallisto](https://pachterlab.github.io/kallisto/starting)

## Test data

The data used to compare the workflows is from [Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie, and Ballgown](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5032908/). Each FASTQ file is named using the SRA RUN IDs.

```
ERR188044
ERR188104
ERR188234
ERR188245
ERR188257
ERR188273
ERR188337
ERR188383
ERR188401
ERR188428
ERR188454
ERR204916
```

You can get more information on the run using [Entrez Direct](https://www.ncbi.nlm.nih.gov/home/tools/). (The run info obtained by running the command below is provided in the `metadata` folder.)

```bash
esearch -db sra -query ERR188044 | efetch -format runinfo
Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,dbgap_study_accession,Consent,RunHash,ReadHash
ERR188044,2012-11-07 04:42:08,2012-11-07 04:41:56,36349964,5525194528,36349964,152,3596,,https://sra-downloadb.st-va.ncbi.nlm.nih.gov/sos2/sra-pub-run-3/ERR188044/ERR188044.1,ERX162864,NA18498.2.M_120131_1 extract,RNA-Seq,cDNA,TRANSCRIPTOMIC,PAIRED,280,0,ILLUMINA,Illumina HiSeq 2000,ERP001942,PRJEB3366,,204869,ERS185292,SAMEA1573216,simple,9606,Homo sapiens,SAMEA1573216,,,,,,,no,,,,,CRG,ERA169774,,public,3DDC6C2865E755D74EBB7702A5BAC58E,D5681D67D5A545BF09827BA3E3C2706D
```

From the metadata we can see the that this run ID belongs to the SRA Study [ERP001942](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?study=ERP001942), which is the "RNA-sequencing of 465 lymphoblastoid cell lines from the 1000 Genomes".

## Processing

First prepare the necessary tools; only `x86_64` and `amd64` architectures are supported.

```console
./scripts/fetch_binaries.sh

# requires various libraries for compiling
./scripts/setup_samtools.sh
./scripts/setup_rsem.sh
```

Next download the testing data and transcript references.

```console
./scripts/fetch_data.sh
```

Create the indexes for the various tools.

```console
./scripts/create_index.sh
```

Run HISAT2 and StringTie2.

```console
./scripts/hisat_stringtie.sh
```

Run STAR and RSEM.

```console
./scripts/star_rsem.sh
```

Run Kallisto.

```console
./scripts/kallisto.sh
```

Results will be in `results`.

## Comparing quantifications

The R Markdown document `compare_quant.Rmd` in `analysis` compares the quantification results.

## nf-core/rnaseq

> [nf-core/rnaseq](https://github.com/nf-core/RNAseq) is a bioinformatics pipeline that can be used to analyse RNA sequencing data obtained from organisms with a reference genome and annotation. It takes a samplesheet and FASTQ files as input, performs quality control (QC), trimming and (pseudo-)alignment, and produces a gene expression matrix and extensive QC report.

To use nf-core/rnaseq, first install `nf-core`.

```console
pip install nf-core
```

Download nf-core/rnaseq and the necessary Singularity images using `nf-core`; the images will be saved in `${HOME}/nf-core/sif`. This takes some time, so go get a coffee/drink/etc.

```console
export NXF_SINGULARITY_CACHEDIR=${HOME}/nf-core/sif

nf-core download rnaseq -r 3.14.0 --outdir ${HOME}/nf-core/rnaseq --compress none --container-system singularity -p 4
```

After downloading run a test; [install](https://www.nextflow.io/docs/latest/install.html) Nextflow if you haven't already.

```console
export NXF_SINGULARITY_CACHEDIR=${HOME}/nf-core/sif
nextflow run ${HOME}/nf-core/rnaseq/3_14_0/main.nf -profile test,singularity --outdir rnaseq_test_output
```

If everything completed successfully, you should see the following:

```
-[nf-core/rnaseq] Pipeline completed successfully -
Completed at: 19-Jul-2024 04:30:21
Duration    : 5m 7s
CPU hours   : 0.4
Succeeded   : 194
```

### Download reference transcriptome

[Download Ensembl references](https://nf-co.re/rnaseq/3.14.0/docs/usage/#reference-genome-options).

```console
RELEASE=$(curl -s 'http://rest.ensembl.org/info/software?content-type=application/json' | grep -o '"release":[0-9]*' | cut -d: -f2)
mkdir release-${RELEASE} && cd release-${RELEASE}
wget -c ftp://ftp.ensembl.org/pub/release-${RELEASE}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget -c "ftp://ftp.ensembl.org/pub/release-${RELEASE}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${RELEASE}.gtf.gz"
```

### Running nf-core/rnaseq on the test data

Generate `samplesheet.csv`.

```console
./scripts/create_samplesheet.pl > samplesheet.csv
```

Run it!

```console
./scripts/run_nfcore_rnaseq_default.sh
```

The results will be in `results/nfcore_rnaseq`.

Use `nextflow log` to check the run log.

```console
nextflow log
```
```
TIMESTAMP          	DURATION	RUN NAME       	STATUS	REVISION ID	SESSION ID                          	COMMAND
2024-08-27 11:06:24	58m 9s  	disturbed_allen	OK    	746820de9b 	e2776595-45cc-49d9-9fa7-52beeed794bb	nextflow run /home/dtang/nf-core/rnaseq/3_14_0/main.nf -resume -with-report execution_report_2024_08_27_11_06_21.html -with-trace -with-dag flowchart_2024_08_27_11_06_21.html --input /home/dtang/github/rnaseq/samplesheet.csv --outdir /home/dtang/github/rnaseq/results/nfcore_rnaseq_default --fasta /home/dtang/github/rnaseq/raw/chrX_data/genome/chrX.fa --gtf /home/dtang/github/rnaseq/raw/chrX_data/genes/gencode.v46.annotation.chrx.gtf --aligner star_rsem -profile singularity --max_cpus 6 --max_memory 16GB
```

Name of task and directory containing results of task.

```console
nextflow log -f name,workdir disturbed_allen | head -3
```
```
NFCORE_RNASEQ:RNASEQ:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:TRIMGALORE (ERR188044)    /home/dtang/github/rnaseq/work/7a/526dc6f5d710e529e2e4bb2e75f15d
NFCORE_RNASEQ:RNASEQ:PREPARE_GENOME:GTF_FILTER (chrX.fa)        /home/dtang/github/rnaseq/work/5e/8f3bd1093fd1d4d1b5947f50eb2f8b
NFCORE_RNASEQ:RNASEQ:PREPARE_GENOME:CUSTOM_GETCHROMSIZES (chrX.fa)      /home/dtang/github/rnaseq/work/3d/e24d393ecf289095608b2df1c41b50
```

### nf-core/rnaseq results

The [following steps](https://nf-co.re/rnaseq/3.14.0/) are run:

1. Merge re-sequenced FastQ files (cat)
2. Sub-sample FastQ files and auto-infer strandedness (fq, Salmon)
3. Read QC (FastQC)
4. UMI extraction (UMI-tools)
5. Adapter and quality trimming (Trim Galore!)
6. Removal of genome contaminants (BBSplit)
7. Removal of ribosomal RNA (SortMeRNA)
8. Choice of multiple alignment and quantification routes:
    * STAR -> Salmon
    * STAR -> RSEM
    * HiSAT2 -> NO QUANTIFICATION
9. Sort and index alignments (SAMtools)
10. UMI-based deduplication (UMI-tools)
11. Duplicate read marking (picard MarkDuplicates)
12. Transcript assembly and quantification (StringTie)
13. Create bigWig coverage files (BEDTools, bedGraphToBigWig)
14. Extensive quality control:
    * RSeQC
    * Qualimap
    * dupRadar
    * Preseq
    * DESeq2
15. Pseudoalignment and quantification (Salmon or ‘Kallisto’; optional)
16. Present QC for raw read, alignment, gene biotype, sample similarity, and strand-specificity checks (MultiQC, R)

To see the steps run, use `nextflow log`.

```console
nextflow log -f name,workdir disturbed_allen | grep ERR188044
```
```
NFCORE_RNASEQ:RNASEQ:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:TRIMGALORE (ERR188044)    /home/dtang/github/rnaseq/work/7a/526dc6f5d710e529e2e4bb2e75f15d
NFCORE_RNASEQ:RNASEQ:FASTQ_SUBSAMPLE_FQ_SALMON:FQ_SUBSAMPLE (ERR188044) /home/dtang/github/rnaseq/work/67/ab6ddac43779f8175031f08d82325d
NFCORE_RNASEQ:RNASEQ:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:FASTQC (ERR188044)        /home/dtang/github/rnaseq/work/f2/cab64c32f2f243f8d15f5a2276a131
NFCORE_RNASEQ:RNASEQ:FASTQ_SUBSAMPLE_FQ_SALMON:SALMON_QUANT (ERR188044) /home/dtang/github/rnaseq/work/1e/c51ad2fbf230d80e6d586959f03f97
NFCORE_RNASEQ:RNASEQ:QUANTIFY_RSEM:RSEM_CALCULATEEXPRESSION (ERR188044) /home/dtang/github/rnaseq/work/ec/a0210b2055fbda3a5538121afde2ce
NFCORE_RNASEQ:RNASEQ:QUANTIFY_RSEM:BAM_SORT_STATS_SAMTOOLS:SAMTOOLS_SORT (ERR188044)    /home/dtang/github/rnaseq/work/f4/f9e5aba2c85ac270e06a2a4bfabd7c
NFCORE_RNASEQ:RNASEQ:QUANTIFY_RSEM:BAM_SORT_STATS_SAMTOOLS:SAMTOOLS_INDEX (ERR188044)   /home/dtang/github/rnaseq/work/37/b9118d82237a5f021723107b68ce20
NFCORE_RNASEQ:RNASEQ:BAM_MARKDUPLICATES_PICARD:PICARD_MARKDUPLICATES (ERR188044)        /home/dtang/github/rnaseq/work/70/f6b8b1bb90aae8c0cd45a3293131b6
NFCORE_RNASEQ:RNASEQ:DUPRADAR (ERR188044)       /home/dtang/github/rnaseq/work/90/85df41c94ac7aefd633014492f3e0d
NFCORE_RNASEQ:RNASEQ:QUANTIFY_RSEM:BAM_SORT_STATS_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_IDXSTATS (ERR188044)     /home/dtang/github/rnaseq/work/b4/56bb47087ee916ef30b47bf6bfac81
NFCORE_RNASEQ:RNASEQ:QUANTIFY_RSEM:BAM_SORT_STATS_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_FLAGSTAT (ERR188044)     /home/dtang/github/rnaseq/work/af/5eed24a33ee7fe53e30f3fc0032bb8
NFCORE_RNASEQ:RNASEQ:QUANTIFY_RSEM:BAM_SORT_STATS_SAMTOOLS:BAM_STATS_SAMTOOLS:SAMTOOLS_STATS (ERR188044)        /home/dtang/github/rnaseq/work/01/165dded30ca01b91eeec9520902378
NFCORE_RNASEQ:RNASEQ:STRINGTIE_STRINGTIE (ERR188044)    /home/dtang/github/rnaseq/work/6d/73dc0e698153e235b2cbe83d05dec0
NFCORE_RNASEQ:RNASEQ:BEDTOOLS_GENOMECOV (ERR188044)     /home/dtang/github/rnaseq/work/f0/11fa5e004db557db47564130a3e5f0
NFCORE_RNASEQ:RNASEQ:QUALIMAP_RNASEQ (ERR188044)        /home/dtang/github/rnaseq/work/0b/791d2f9f8acf9a0f3c6372b645d438
NFCORE_RNASEQ:RNASEQ:BAM_MARKDUPLICATES_PICARD:SAMTOOLS_INDEX (ERR188044)       /home/dtang/github/rnaseq/work/ba/e3c3e01fc176262600221c733c67cb
NFCORE_RNASEQ:RNASEQ:BAM_MARKDUPLICATES_PICARD:BAM_STATS_SAMTOOLS:SAMTOOLS_IDXSTATS (ERR188044) /home/dtang/github/rnaseq/work/ae/3b941048cf11aed9c8d00e8dfbbd43
NFCORE_RNASEQ:RNASEQ:BAM_MARKDUPLICATES_PICARD:BAM_STATS_SAMTOOLS:SAMTOOLS_FLAGSTAT (ERR188044) /home/dtang/github/rnaseq/work/f2/6ea4c1d1eed6db4a9888fa96e73173
NFCORE_RNASEQ:RNASEQ:BAM_MARKDUPLICATES_PICARD:BAM_STATS_SAMTOOLS:SAMTOOLS_STATS (ERR188044)    /home/dtang/github/rnaseq/work/64/601ee5ecd8771d78d1d9c321c05a0c
NFCORE_RNASEQ:RNASEQ:BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD:UCSC_BEDCLIP (ERR188044) /home/dtang/github/rnaseq/work/58/63a013120761848763aa60316ec4b8
NFCORE_RNASEQ:RNASEQ:BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE:UCSC_BEDCLIP (ERR188044) /home/dtang/github/rnaseq/work/32/733fc1693e0c65809b4b50f7d688d6
NFCORE_RNASEQ:RNASEQ:BAM_RSEQC:RSEQC_JUNCTIONSATURATION (ERR188044)     /home/dtang/github/rnaseq/work/72/26a223bb06dcf4cc63747f3a9468dc
NFCORE_RNASEQ:RNASEQ:BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD:UCSC_BEDGRAPHTOBIGWIG (ERR188044)        /home/dtang/github/rnaseq/work/c7/132cd71ae1ddf515f63aebae9996a1
NFCORE_RNASEQ:RNASEQ:BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE:UCSC_BEDGRAPHTOBIGWIG (ERR188044)        /home/dtang/github/rnaseq/work/8c/60897426d2c86078322eab18334f58
NFCORE_RNASEQ:RNASEQ:BAM_RSEQC:RSEQC_READDUPLICATION (ERR188044)        /home/dtang/github/rnaseq/work/17/270bd314db35baa68102828396b1db
NFCORE_RNASEQ:RNASEQ:BAM_RSEQC:RSEQC_READDISTRIBUTION (ERR188044)       /home/dtang/github/rnaseq/work/ee/39e9d4b1cd296b3578c81ddaf70095
NFCORE_RNASEQ:RNASEQ:BAM_RSEQC:RSEQC_INFEREXPERIMENT (ERR188044)        /home/dtang/github/rnaseq/work/a3/8cb64c313ee38b918772a4050e882f
NFCORE_RNASEQ:RNASEQ:BAM_RSEQC:RSEQC_JUNCTIONANNOTATION (ERR188044)     /home/dtang/github/rnaseq/work/fa/533b52517472114de44f779a72ab1a
NFCORE_RNASEQ:RNASEQ:BAM_RSEQC:RSEQC_BAMSTAT (ERR188044)        /home/dtang/github/rnaseq/work/24/429d85a4510077beef6a5ea22a6d11
NFCORE_RNASEQ:RNASEQ:BAM_RSEQC:RSEQC_INNERDISTANCE (ERR188044)  /home/dtang/github/rnaseq/work/f3/e3b3b7f2fb50e338b315f305f5f71a
```

[STAR + RSEM](https://nf-co.re/rnaseq/3.14.0/docs/output/#star-via-rsem) results are :

* `rsem.merged.gene_counts.tsv`: Matrix of gene-level raw counts across all samples.
* `rsem.merged.gene_tpm.tsv`: Matrix of gene-level TPM values across all samples.
* `rsem.merged.transcript_counts.tsv`: Matrix of isoform-level raw counts across all samples.
* `rsem.merged.transcript_tpm.tsv`: Matrix of isoform-level TPM values across all samples.
* `*.genes.results`: RSEM gene-level quantification results for each sample.
* `*.isoforms.results`: RSEM isoform-level quantification results for each sample.

Raw counts.

```console
head -3 results/nfcore_rnaseq/star_rsem/rsem.merged.gene_counts.tsv
```
```
gene_id	transcript_id(s)	ERR188044	ERR188104	ERR188234	ERR188245	ERR188257	ERR188273	ERR188337	ERR188383	ERR188401	ERR188428	ERR188454	ERR204916
ENSG00000000003.16	ENST00000373020.9,ENST00000494424.1,ENST00000496771.5,ENST00000612152.4	38.00	22.00	3.00	0.00	8.00	8.00	22.00	12.00	9.00	1.00	5.00	5.00
ENSG00000000005.6	ENST00000373031.5,ENST00000485971.1	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00
```

TPM normalised.

```console
head -3 results/nfcore_rnaseq/star_rsem/rsem.merged.gene_tpm.tsv
```
```
gene_id	transcript_id(s)	ERR188044	ERR188104	ERR188234	ERR188245	ERR188257	ERR188273	ERR188337	ERR188383	ERR188401	ERR188428	ERR188454	ERR204916
ENSG00000000003.16	ENST00000373020.9,ENST00000494424.1,ENST00000496771.5,ENST00000612152.4	10.83	8.73	0.57	0.00	5.02	4.26	12.98	3.70	2.12	0.38	1.45	4.29
ENSG00000000005.6	ENST00000373031.5,ENST00000485971.1	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00
```

### Creating a reference

`nf-core/rnaseq` provides the option `--save_reference` to save the STAR index in the results directory. This is useful because you don't want to create the same index each time you run the workflow. From the documentation it isn't clear whether we can just run the workflow to create the reference. The next best thing to do is to simply skip all the steps of the workflow, except the one that creates the reference. This is what `scripts/create_nfcore_rnaseq_ref.sh` does; however we can't use `--skip_alignment` because this will skip all alignment steps.

Unfortunately, the reference does not seem to be saved in the specified results directory.

## Papers

Papers to read when deciding choice of tool, gene mdoels, and gene quantification method for RNA-seq experiments.

* [A survey of best practices for RNA-seq data analysis](https://pubmed.ncbi.nlm.nih.gov/26813401/)
* [Alignment and mapping methodology influence transcript abundance estimation](https://www.biorxiv.org/content/10.1101/657874v2)
* [A comprehensive evaluation of ensembl, RefSeq, and UCSC annotations in the context of RNA-seq read mapping and gene quantification](https://pubmed.ncbi.nlm.nih.gov/25765860/)
* [A benchmark for RNA-seq quantification pipelines](https://pubmed.ncbi.nlm.nih.gov/27107712/)

Also checkout this [list of benchmarks](https://github.com/j-andrews7/awesome-bioinformatics-benchmarks#rna-seq).
