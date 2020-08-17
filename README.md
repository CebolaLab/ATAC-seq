# ATAC-seq
Step-by-step analysis pipeline for ATAC-seq data

The following pipeline will describe the step-by-step analysis of ATAC-seq data (the **a**ssay for **t**ransposase-**a**ccessible **c**hromatin with **seq**uencing). This has been adapted from the following resources:

- https://vallierlab.wixsite.com/pipelines/atac-seq
- https://www.encodeproject.org/atac-seq/ *the recommended ENCODE pipeline*
- https://github.com/ENCODE-DCC/atac-seq-pipeline *the ENCODE pipeline available on github*
- https://github.com/harvardinformatics/ATAC-seq
- https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit

An excellent recent review on the ATAC-seq analysis pipeline is reported by [(Yan et al. 2020)](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-1929-3).

The following steps will be covered:

- Sequencing quality control (QC) 
- Alignment 
- Peak Calling
- Visualisation
- Functional analysis & Motif Discovery


## Pre-alignment QC

[FastQC reports](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf) can be generated for all samples to assess sequence quality, GC content, duplication rates, length distribution, K-mer content and adapter contamination. In ATAC-seq data, it is likely that Nextera sequencing adapters will be over-represented. As described by [(Yan et al. 2020)](https://genomebiology.biomedcentr\
al.com/track/pdf/10.1186/s13059-020-1929-3), base quality should be high although may drop slightly at the 3' end, while GC content and read length should be consistent with the expected values. For paired-end reads, run fastqc on both files:

`fastqc "$base"_1.fastq.gz -d . -o .` 
`fastqc "$base"_2.fastq.gz -d . -o .`

Adapters and low quality reads/bases can be trimmed using one of several programs, such as [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) or [cutadapt](https://cutadapt.readthedocs.io/en/stable/). 

Trim to a fixed length yes or no? Vallier lab says yes, paper suggests that one of the advantages is shorter fragments?

`trimmomatic ... `

A QC report can be generated following trimming to compare the quality before and after trimming.

`fastqc "$base"_1_trimmed.fastq.gz -d . -o .`  

## Alignment

The processed reads should then be aligned to the reference human genome using an aligner such as [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) 


`bowtie2 -x $bt2idx/hg19.masked -1 "$base"_1.paired.fastq.gz -2 "$base"_2.paired.fastq.gz) 2> "$base".bowtie2 | samtools view -bS - > "$base"_aligned_reads.bam`

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