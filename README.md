# ATAC-seq
## **UNDER CONSTRUCTION**

[Cebola Lab](https://www.imperial.ac.uk/metabolism-digestion-reproduction/research/systems-medicine/genetics--genomics/regulatory-genomics-and-metabolic-disease/)

Step-by-step analysis pipeline/tutorial for ATAC-seq data

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

Samples which have been sequenced on multiple lanes can be concatenated at this stage [optional].

```
cat  <sample>.lane1.R1.fastq.gz  <sample>.lane2.R1.fastq.gz  >  <sample>.R1.fastq.gz

cat  <sample>.lane1.R2.fastq.gz  <sample>.lane2.R2.fastq.gz  >  <sample>.R2.fastq.gz
```

Independent replicates should be processed separately.

### Initial QC report

The raw sequence data should first be assessed for quality. [FastQC reports](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf) can be generated for all samples to assess sequence quality, GC content, duplication rates, length distribution, K-mer content and adapter contamination. In ATAC-seq data, it is likely that Nextera sequencing adapters will be over-represented. As described by [Yan et al. (2020)](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-1929-3), base quality should be high although may drop slightly at the 3' end, while GC content and read length should be consistent with the expected values. For paired-end reads, run fastqc on both files, with the results output to the current directory:

```
fastqc <sample>_R1.fastq.gz -d . -o .

fastqc <sample>_R2.fastq.gz -d . -o .
```

### Adapter trimming 

Adapters and low quality reads/bases should be trimmed using one of several programs, such as [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [cutadapt](https://cutadapt.readthedocs.io/en/stable/), or [fastp](https://github.com/OpenGene/fastp). Adapter contamination can be seen in the fastqc report:

<img src="https://github.com/CebolaLab/ATAC-seq/blob/master/Figures/adapters.png" width="800">

For this pipeline, fastp is used to remove adapter sequences. The minimum fragment length is set at 35, since short ATAC-seq fragments can be observed if the transposase cuts adjacent nucleosome-free DNA. 

```bash
fastp -i <sample>_R1.fastq.gz -I <sample>_R2.fastq.gz -o <sample>_R1.trimmed.fastq.gz -O <sample>_R2.trimmed.fastq.gz --detect_adapter_for_pe -l 35 -j <sample>.fastp.json -h <sample>.fastp.html
```

The output of fastp includes a html report, part of which is shown below. This presents the total number of reads before and after filtering, including the % of high quality (Q30) bases. The report also shows the main causes of read removal. In the example below, 1.9% of reads were removed because they were shorter than the minimum read length specified above by the -l argument (35bp).

<img src="https://github.com/CebolaLab/ATAC-seq/blob/master/Figures/fastp-html.png" width="600">

## Alignment

The processed reads should then be aligned to the reference human genome using an aligner such as [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). This pipeline will use bowtie2 to align reads to the hg19 reference genome. If the user is aligning to the more recent GRCh38 release, it is recommended to remove alternative contigs, otherwise reads may not map uniquely and will be assigned a low quality score. Suggested guidelines for preparing the GRCh38 genome are discussed in [this tutorial](https://www.biostars.org/p/342482/). If the user selects an alternative alignment tool, such as bwa, they are referred to [this blog post](https://www.acgt.me/?offset=1426809676847) which discusses the resulting differences in alignment quality scores.

#### Bowtie2 alignment

The `local` parameter is used to 'soft clip' the end of reads to allow the best possible alignment, including any remaining adapter sequences (e.g. 1 or 2bp).  By using the `--no-mixed` and `--no-discordant` parameters, reads will only be aligned if both reads align successfully as a pair (this avoids the need to later remove reads which are not properly paired, which is a common post-alignment QC step). The `-I 35` and `-X 700` require fragments to be greater than 35bp and less than 700bp. The user can adjust these depending on the experiment/sequencing protocol (see the fastp html report for a plot of the estimated insert sizes). The maximum fragment length of 700bp prevents reads from aligning incorrectly outside the expected fragment range. 

***Which reference genome to use?*** See [this](http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) discussion on which reference genome to use. The recommended downloads for both hg19/b37 and GRCh38 are shown below.


```bash
#to download the build 37 (hg19) reference genome
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz

#to download the more recent GRCh38
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
```

Bowtie2 should be used to create the reference genome index files (see the bowtie2 [manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome)). After the index files have been generated, align the trimmed ATAC-seq fastq files to the genome (here using hg19/b37):

```bash
#set the bt2idx variable to the directory with the reference genome and indexes
bt2idx=/path/to/reference-genome

#Run the bowtie2 alignment and output a bam alignment file
bowtie2 --local --very-sensitive --no-mixed --no-discordant -I 35 -X 700 -x $bt2idx/human_g1k_v37.fasta -1 <sample>_R1.trimmed.fastq.gz -2 <sample>_R2.trimmed.fastq.gz) 2> <sample>.bowtie2 | samtools view -bS - > <sample>.bam
```

The output `bam` file should be sorted and indexed prior to the downstream analysis:

```bash
#Sort the output bam file by coordinate
picard SortSam I=<sample>.bam O=<sample>_sorted.bam SO=coordinate CREATE_INDEX=TRUE
```

## Post-alignment QC

The post-alignment QC steps involve several steps:

- [Remove mitochondrial reads](#remove-mitochondrial-reads)
- [Remove duplicates & low-quality alignments](#remove-duplicates-&-low-quality-alignments) (including non-uniquely mapped reads)
- [Calculate library complexity and QC](#calculate-library-complexity-and-QC)
- [Remove ENCODE black-list regions](#remove-encode-black-list-regions)

For an ATAC-seq experiment, the number of uniquely mapped reads ***after these steps*** is recommended to be 25 million of 50 million paired-end reads [(Buenrostro et al. 2015)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4374986/). Specific to ATAC-seq, an additional QC step is to check the fragment size distribution, which is expected to correspond to the length of nucleosomes:

- [Assess fragment size distribution](#assess-fragment-size-distribution)

### Remove mitochondrial reads

ATAC-seq experiments commonly include a high proportion of mitochondrial reads. These should be removed. To assess the total % of mitochondrial reads, `samtools idxstats` can be run to report the total number of reads mapping to each chromosome. `samtools flagstat` provides a short report including the total number of DNA fragments. 

```bash
#Generate the idxstats report
samtools idxstats <sample>_sorted.bam > <sample>_sorted.idxstats

#Check the number of reads mapped to the mitochondria (chrM)
grep "chrM" <sample>_sorted.idxstats
```

The second column is the length of the chromosome and the third column is the total number of reads aligned to the chromosome (chrM). To see the total number of DNA fragments, run:

```bash
#Generate the flagstat report
samtools flagstat <sample>_sorted.bam > <sample>_sorted.flagstat

head <sample>_sorted.flagstat
```

The % of DNA fragments aligned to chrM can be calculated as a % of the total DNA fragments. To remove any mitocondrial DNA, run the following:

```bash
#Remove reads aligned to the mitochondria
samtools view -h <sample>-sorted.bam | grep -v chrM | samtools sort -O bam -o <sample>.rmChrM.bam -T .
```

### Mark duplicates 

To mark duplicate reads:

```bash
picard MarkDuplicates QUIET=true INPUT=<sample>.rmChrM.bam OUTPUT=<sample>.marked.bam METRICS_FILE=<sample>.dup.metrics REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=.
```

The % of duplicates can be viewed using:

```bash
head -n 8 <sample>.dup.metrics | cut -f 7,9 | grep -v ^# | tail -n 2
```

### Remove duplicates & low-quality alignments 

The output `sam/bam` files contain several measures of quality. First, the alignment quality score. Reads which are uniquely mapped are assigned a high alignment quality score and one genomic position. If reads can map to more than one location, Bowtie2 reports one position and assigns a low quality score. The proportion of uniquely mapped reads can be assessed. In general, >70% uniquely mapped reads is expected, while <50% may be a cause for concern [(Bailey et al. 2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3828144/pdf/pcbi.1003326.pdf). Secondly, the 'flag' reports information such as whether the read is mapped as a pair or is a PCR duplicate. The individual flags are reported [here](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/04_alignment_quality.html) and are combined in a `sam/bam` file to one score, which can be deconstructed back to the original flags using [online interpretation tools](https://broadinstitute.github.io/picard/explain-flags.html). In this pipeline, the bowtie2 parameters `--no-mixed` and `--no-discordant` prevent the mapping of only one read in a pair, so these flags will not be present. All flags reported in a `sam` file can optionally be viewed using  `grep -v ^@ <sample>.sam | cut -f 2 | sort | uniq`.

**Multi-mapping:** the user should decide whether or not to retain reads which have multi-mapped, i.e. aligned to more than one position in the reference genome. When using paired-end data, it may be the case that one read aligns to a repetitive region (and therefore can map elsewhere), while the mate aligns to a unique sequence with a high quality. The bowtie2 parameters used above required reads to align within 50-700bp, so there should be no reads incorrectly aligned outside this distance. As such, the user may decide to keep multi-mapping reads on the assumption that they are likely to be mapped to the correct sequence, within the length of the DNA fragment. This may, however, cause incorrect alignments in extended repetitive regions where a read could map to multiple positions within the length of the DNA fragment. This should be minimised by the downstream removal of the [ENCODE Blacklisted regions](https://www.nature.com/articles/s41598-019-45839-z).

If a read is multi-mapped, it is assigned a low quality score by bowtie2. To view how many DNA reads which align with a quality score >30, run (divide this number by 2 to calculate the # of DNA fragments):

```bash
samtools view -q 30 -c <sample>.marked.bam
```

A low % of uniquely mapped reads map result from short reads, excessive PCR amplification or problems with the PCR (Bailey et al. 2013). The following code uses the `sam/bam` flags to retain properly mapped pairs (`-f 2`) and to remove reads which fail the platform/vendor QC checks (`-F 512`), duplicate reads (`-F 1024`) and those which are unmapped (`-F 12`). The three flags to be removed can be combined into `-F 1548`, which will remove reads which meet any of the three individual flags

To ***retain*** multi-mapped reads:

```bash
samtools view -h -b -f 2 -F 1548 <sample>.rmChrM.bam | samtools sort -n <sample>.filtered.bam 
```

To ***remove*** multi-mapped reads:

```bash
samtools view -h -b -f 2 -F 1548 -q 30 <sample>.rmChrM.bam | samtools sort -n <sample>.filtered.bam
```

The output `bam` file, which is now sorted by name, should be indexed: 

```bash
samtools index <sample>.filtered.bam
```

### Remove ENCODE black-list regions 

```bash
bedtools intersect -nonamecheck -v -abam <sample>.filtered.bam -b ${BLACKLIST} > <sample>.blacklist-filtered.bam
```

### Assess fragment size distribution and QC

The fragment size is expected to show a periodicity of 150/200 bp, reflecting the length of DNA surrounding nucleosomes, since the tagmentation typically cuts DNA in between nucleosomes. The `R` tool [ATACseqQC](https://www.bioconductor.org/packages/release/bioc/html/ATACseqQC.html) can be used to assess the distribution of fragments. Fragments sizes may also be <100bp, corresponding to nucleosome-free regions or cutting of linker DNA, as well as fragments of and mono-, di-, and tri-nucleosomes (~200, 400 and 600bp, respectively) [(Yan et al. 2020)](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-1929-3). Fragments from nucleosome-free regions are expected to be enriched around transcription-start sites (TSS) and fragments from nucleosome-bound regions are depleted around TSS and maybe be slightly enriched in the flanking regions. The ATACseqQC manual is available [here](https://www.bioconductor.org/packages/release/bioc/vignettes/ATACseqQC/inst/doc/ATACseqQC.html). A guide to R ATACseqQC is reported in a [seperate markdown] and covers:


- Estimate library complexity
- Fragment size distribution 
- GC bias
- Nucleosome positioning 
- Plot footprints
- Plot correlations between samples

## Peak calling

Peaks are identified where sequenced reads accumulate. These correspond to regions of accessible DNA.

Important considerations for ATAC-seq: there are usually no controls; the Tn5 transposase has a binding preference, resulting in a GC bias which should be corrected for during peak calling; repair of the transposase-induced nick introduces a 9bp insertion. When investigating high-resolution motif enrichment, ATAC-seq reads are shifted +4bp and -5bp for positive and negative strands to account for this. This step is not necessary for peak calling, since hundreds of bp are typically implicated.  [Yan et al. (2020)](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-1929-3) recommend using the peak caller HMMRATAC, developed by [Tarbell et al. (2019)](https://academic.oup.com/nar/article/47/16/e91/5519166), which is specifically developed for ATAC-seq data (if computational resources are sufficient). Alternatively, MACS2 is a pop\lar peak caller.


HMMRATAC is available on [github](https://github.com/LiuLabUB/HMMRATAC).

The input file should be sorted by coordinate and indexed.

```bash
samtools view -H <sample>.filtered.bam | perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > genome.info

java -jar HMMRATAC_V1.2.10_exe.jar -b <sample>.filtered.bam -i <sample>.filtered.bam.bai -g genome.info -o <sample>
```

## Visualisation

The QC-ed `bam` file can be converted to a `bedGraph` format to visualise sequencing trakcs using tools such as the UCSC browser or the integrative genomes browser. The ENCODE blacklist regions can be provided, to exclude them from the output:

```bash
bedtools bamCoverage --blackListFileName --normalizeUsing BPM -b <sample>.filtered.bam > <sample>.bedGraph
```

TmpScale=$(bc <<< "scale=6;1000000/$(samtools view -f 0 -c mybamfile.bam)")
bedtools genomecov -ibam mybamfile.bam -bg -scale $TmpScale -g mm10.chrom.sizes > myBedGraphfile.BedGraph

The `-bg` option can be replaced with `-bga` if the user wants the output file to contain regions with 0 coverage. 

## Additional analysis: motif calling

An important step with ATAC-seq data is to shift reads +4bp and -5bp for positive and negative strands, due to the 9bp duplication introducted through the repair of the Tn5 transposase nick.

A recent tool which can be used to assess motifs and transcription factor footprints is [TOBIAS](https://github.com/loosolab/TOBIAS).

```bash
picard CollectInsertSizeMetrics
```

