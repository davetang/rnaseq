## Table of Contents

  - [RNA-seq pipelines](#rna-seq-pipelines)
  - [Setup](#setup)
  - [Test data](#test-data)
  - [Comparing quantifications](#comparing-quantifications)
  - [nf-core/rnaseq](#nf-corernaseq)
  - [Papers](#papers)

## RNA-seq pipelines

* [HISAT2 + StringTie2 pipeline](https://davetang.org/muse/2017/10/25/getting-started-hisat-stringtie-ballgown/)
* [STAR, Cufflinks, RSEM](https://pubmed.ncbi.nlm.nih.gov/27662878/)
* [Kallisto](https://pachterlab.github.io/kallisto/starting)

## Setup

1. In `raw` run `./fetch_data.sh` and then `./create_index.sh`.
2. In `src` run `./fetch_binaries.sh` (requires macOS or Linux), `./setup_samtools.sh`, and `./setup_rsem.sh` (requires various libraries for compiling).
3. HISAT2 and StringTie2 can be run from `hisat_stringtie` by running `./map.sh` and then `./quant.sh`.
4. STAR and RSEM can be run from `star_rsem` by running `./run.sh`.
5. Kallisto can be run from `kallisto` by running `./quant.sh`.

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

## Comparing quantifications

The R Markdown document `compare_quant.Rmd` in `analysis` compares the quantification results.

## nf-core/rnaseq

> [nf-core/rnaseq](https://github.com/nf-core/RNAseq) is a bioinformatics pipeline that can be used to analyse RNA sequencing data obtained from organisms with a reference genome and annotation. It takes a samplesheet and FASTQ files as input, performs quality control (QC), trimming and (pseudo-)alignment, and produces a gene expression matrix and extensive QC report.

To use nf-core/rnaseq, first install `nf-core`.

```console
pip install nf-core
```

Download nf-core/rnaseq and the necessary Singularity images using `nf-core`. This takes some time, so go get a coffee/drink/etc.

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

## Papers

Papers to read when deciding choice of tool, gene mdoels, and gene quantification method for RNA-seq experiments.

* [A survey of best practices for RNA-seq data analysis](https://pubmed.ncbi.nlm.nih.gov/26813401/)
* [Alignment and mapping methodology influence transcript abundance estimation](https://www.biorxiv.org/content/10.1101/657874v2)
* [A comprehensive evaluation of ensembl, RefSeq, and UCSC annotations in the context of RNA-seq read mapping and gene quantification](https://pubmed.ncbi.nlm.nih.gov/25765860/)
* [A benchmark for RNA-seq quantification pipelines](https://pubmed.ncbi.nlm.nih.gov/27107712/)

Also checkout this [list of benchmarks](https://github.com/j-andrews7/awesome-bioinformatics-benchmarks#rna-seq).
