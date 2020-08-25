# ATAC-seq
Step-by-step analysis pipeline for ATAC-seq data

The following pipeline will describe the step-by-step analysis of ATAC-seq data (the **a**ssay for **t**ransposase-**a**ccessible **c**hromatin with **seq**uencing). This has been adapted from the following resources:

- https://vallierlab.wixsite.com/pipelines/atac-seq A great tool for beginners indicating the major analysis steps 
- https://www.encodeproject.org/atac-seq/ The recommended ENCODE pipeline, for which tools are available on [github](https://github.com/ENCODE-DCC/atac-seq-pipeline) and the recommended parameters/specification are available via a [google doc](https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit)
- https://github.com/harvardinformatics/ATAC-seq An ATAC-seq pipeline from Harvard Informatics 
- https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/ A similar github page presenting an ATAC-seq pipeline 
- https://galaxyproject.github.io/training-material/topics/epigenetics/tutorials/atac-seq/tutorial.html

An recent review on the ATAC-seq analysis pipeline is reported by [Yan et al. (2020)](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-1929-3).

The following steps will be covered:

- [Pre-alignment quality control (QC)](#pre-alignment-qc) 
- [Alignment](#alignment) 
- [Post-alignment QC](#post-alignment-qc)
- [Peak Calling](#peak-calling)
- [Visualisation](#visualisation)
- Functional analysis & Motif Discovery

## Pre-alignment QC

Firstly, if the same sample has been sequenced on multiple lanes, concatenate the files:

```
cat  <sample>.lane1.R1.fastq.gz  <sample>.lane2.R1.fastq.gz  >  <sample>.R1.fastq.gz

cat  <sample>.lane1.R2.fastq.gz  <sample>.lane2.R2.fastq.gz  >  <sample>.R2.fastq.gz
```

Independent replicates should be processed separately.

### Initial QC report

The raw sequence data should first be assessed for quality. [FastQC reports](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf) can be generated for all samples to assess sequence quality, GC content, duplication rates, length distribution, K-mer content and adapter contamination. In ATAC-seq data, it is likely that Nextera sequencing adapters will be over-represented. As described by [(Yan et al. 2020)](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-1929-3), base quality should be high although may drop slightly at the 3' end, while GC content and read length should be consistent with the expected values. For paired-end reads, run fastqc on both files, with the results output to the current directory:

```
fastqc <sample>_1.fastq.gz -d . -o .

fastqc <sample>_2.fastq.gz -d . -o .
```

### Adapter trimming 

Adapters and low quality reads/bases can be trimmed using one of several programs, such as [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [cutadapt](https://cutadapt.readthedocs.io/en/stable/), or [fastp](https://github.com/OpenGene/fastp). For this pipeline, the user is advised to use fastp.

```
fastp -i <sample>_R1.fastq.gz -O <sample>_R1.trimmed.fastq.g -I <sample>_R2.fastq.gz -O <sample>_R2.trimmed.fastq.gz --detect_adapter_for_pe -l 50 -j <sample>.fastp.json -h <sample>.fastp.html
```

<img src="https://github.com/CebolaLab/ATAC-seq/blob/master/Figures/fastp-html.png" width="600">

The user should specify whether the raw data is encoded in phred+33 or phred+63 (read more [here](https://sequencing.qcfail.com/articles/incorrect-encoding-of-phred-scores/)). Most data should be encoded in the standardised phred+33. This can be confirmed using the fastQC report generated previously: the 'Encoding' field should read Sanger / Illumina 1.9 as below: 

<img src="https://github.com/CebolaLab/ATAC-seq/blob/master/Figures/fastqc1.png" width="500">

The phred encoding is specified when using trimmomatics to trim adapters. In the following example, the Nextera transposase adapter sequences are saved in the file `nextera-adapters.fa`. Other Illumina adapters sequences can be found [here](https://dnatech.genomecenter.ucdavis.edu/wp-content/uploads/2013/06/illumina-adapter-sequences_1000000002694-00.pdf).

`trimmomatic PE -phred33 -trimlog <sample>.trimLogFile <sample>_R1.fastq.gz <sample>_R2.fastq.gz <sample>_R1_trimmed_paired.fastq.gz <sample>_R1_trimmed_unpaired.fastq.gz <sample>_R2_trimmed_paired.fastq.gz <sample>_R2_trimmed_unpaired.fastq.gz ILLUMINACLIP:nextera-adapters.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36`

A QC report can be generated following trimming to compare the quality before and after trimming.

`fastqc <sample>_1_trimmed.fastq.gz -d . -o .`  

## Alignment

The processed reads should then be aligned to the reference human genome using an aligner such as [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). This pipeline will use bowtie2 to align reads to the hg19 reference genome. If the user is aligning to the more recent GRCh38 release, it is recommended to remove alternative contigs, otherwise reads may not map uniquely and will consequently be assigned a low quality score. Suggested guidelines for preparing the GRCh38 genome are discussed in [this tutorial](https://www.biostars.org/p/342482/). If the user selects an alternative alignment tool, such as bwa, they are referred to [this blog post](https://www.acgt.me/?offset=1426809676847) which discusses the resulting difference in alignment quality scores (this may require an adjustment in downstream filtering).  

The `local` parameter is used, since cutadapt removes adapter sequences >3bp; the `local` parameter 'soft clips' the end of reads to allow the best possible alignment, including any remaining adapter sequences (1 or 2bp).  By using the `--no-mixed` and `--no-discordant` parameters, reads will only be aligned if both reads align successfully as a pair (this avoids the need to later remove reads which are not properly paired, which is a common post-alignment QC step).

`bowtie2 --local --very-sensitive --no-mixed --no-discordant -x $bt2idx/hg19.masked -1 <sample>_1.paired.fastq.gz -2 <sample>_2.paired.fastq.gz) 2> <sample>.bowtie2 | samtools view -bS - > <sample>_aligned_reads.bam`

## Post-alignment QC

Remove reads aligned to the mitochondria or the ENCODE blacklisted regions, as well as reads with a low mapping quality, those which are inproperly paired and duplicate reads.

The number of uniquely mapped reads after these steps is recommended to be 25 million of 50 million paired-end reads

Specific to ATAC-seq, an additional QC step is to check the fragment size distribution, which is expected to correspond to the 

samtools view`.
Remove PCR duplicated using `samtools rmdup`
Plot fragment size using `ATACseqQC`

### ATAC-seq QC 

The QC tool [ATACseqQC](https://www.bioconductor.org/packages/release/bioc/html/ATACseqQC.html) can be used to assess the distribution of fragments. For ATAC-seq data, fragments size should correspond to the size of nucleosomes (150/200bp) and to <100bp in nculeosome-free regions. 

nucleosome-free regions (NFR), with fragments of <100bp, and mono-, di-, and tri-nucleosomes (~200, 400 and 600bp, respectively) [(Yan et al. 2020)](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-1929-3). Fragments are expected to be enriched around transcription-start sites (TSS) (corresponding to NFRs) and fragments from nucleosome-bound regions are depleted around TSS and maybe be slightly enriched in the flanking regions. 

 periodicity of 150/200 bp (nucleosome size),

To account for the 9bp duplication created by the Tn5 transposase, reads should be shifted +4bp for the positive strand and -5bp for the negative strand. 

## Peak calling  