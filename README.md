# ATAC-seq
## **UNDER CONSTRUCTION**

A step-by-step pipeline for the analysis of ATAC-seq data (the **a**ssay for **t**ransposase-**a**ccessible **c**hromatin with **seq**uencing), written by the [Cebola Lab](https://www.imperial.ac.uk/metabolism-digestion-reproduction/research/systems-medicine/genetics--genomics/regulatory-genomics-and-metabolic-disease/). The resources used to build this pipeiline are listed at the bottom in the [Resources](#resources) section. 

Correspondence = hannah.maude12@imperial.ac.uk

The following steps will be covered:


**Pre-alignment processing**
- [Pre-alignment quality control (QC)](#pre-alignment-qc): QC reports and adapter trimming

**Alignment, QC and track visualisation**
- [Alignment](#alignment) 
- [Post-alignment QC](#post-alignment-qc): filter, check library complexity and format for peak calling
- [Bam visualisation, bam to bigWig](#bam-visualisation): generate tracks to visualise the aligned data on a genome browser

**Peak calling, QC and visualisation**
- [Peak calling](#peak-calling): call peaks using MACS2, HMMRATAC and Genrich
- [Peak calling QC](#peak-qc): FRiP, transcription start site enrichment, sample reproducibility
- [Peak visualisation](#peak-visualisation): p-value, peaks and pileup

**Differential accessibility analysis**
- [Differential accessibility (DA) analysis](#peak-QC-and-DA)
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

For this pipeline, fastp is used to remove adapter sequences. The minimum fragment length is set at 25, since short ATAC-seq fragments can be observed if the transposase cuts adjacent nucleosome-free DNA. 

```bash
fastp -i <sample>_R1.fastq.gz -I <sample>_R2.fastq.gz -o <sample>_R1.trimmed.fastq.gz -O <sample>_R2.trimmed.fastq.gz --detect_adapter_for_pe -l 25 -j <sample>.fastp.json -h <sample>.fastp.html
```

The output of fastp includes a html report, part of which is shown below. This presents the total number of reads before and after filtering, including the % of high quality (Q30) bases. The report also shows the main causes of read removal. In the example below, 1.9% of reads were removed because they were shorter than the minimum read length specified above by the -l argument (25bp).

<img src="https://github.com/CebolaLab/ATAC-seq/blob/master/Figures/fastp-html.png" width="600">

## Alignment

The processed reads should then be aligned to the reference human genome using an aligner such as [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). This pipeline will use bowtie2 to align reads to the hg19 reference genome. If the user is aligning to the more recent GRCh38 release, it is recommended to remove alternative contigs, otherwise reads may not map uniquely and will be assigned a low quality score. Suggested guidelines for preparing the GRCh38 genome are discussed in [this tutorial](https://www.biostars.org/p/342482/). If the user selects an alternative alignment tool, such as bwa, they are referred to [this blog post](https://www.acgt.me/?offset=1426809676847) which discusses the resulting differences in alignment quality scores.

#### Bowtie2 alignment

The `local` parameter is used to 'soft clip' the end of reads to allow the best possible alignment, including any remaining adapter sequences (e.g. 1 or 2bp).  By using the `--no-mixed` and `--no-discordant` parameters, reads will only be aligned if both reads align successfully as a pair (this avoids the need to later remove reads which are not properly paired, which is a common post-alignment QC step). The `-I 25` and `-X 700` require fragments to be greater than 25bp and less than 700bp. The user can adjust these depending on the experiment/sequencing protocol (see the fastp html report for a plot of the estimated insert sizes). The maximum fragment length of 700bp prevents reads from aligning incorrectly outside the expected fragment range. 

***Which reference genome to use?*** See [this](http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) discussion on which reference genome to use. The recommended downloads for both hg19/b37 and GRCh38 are shown below.


```bash
#to download the build 37 (hg19) reference genome
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz

#to download the more recent GRCh38
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

#Remember to extract the files using gunzip!
```

Bowtie2 should be used to create the reference genome index files (see the bowtie2 [manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome)). After the index files have been generated, align the trimmed ATAC-seq fastq files to the genome:

```bash
#set the bt2idx variable to the directory with the reference genome and indexes
bt2idx=/path/to/reference-genome
```

To align to hg19/b37:
```bash
#Run the bowtie2 alignment and output a bam alignment file
bowtie2 --local --very-sensitive --no-mixed --no-discordant -I 25 -X 700 -x $bt2idx/human_g1k_v37 -1 <sample>_R1.trimmed.fastq.gz -2 <sample>_R2.trimmed.fastq.gz | samtools view -bS - > <sample>.bam
```

To align to GRCh38:
```bash
bowtie2 --local --very-sensitive --no-mixed --no-discordant -I 25 -X 700 -x $bt2idx/GCA_000001405.15_GRCh38_no_alt_analysis_set -1 <sample>_R1.trimmed.fastq.gz -2 <sample>_R2.trimmed.fastq.gz | samtools view -bS - > <sample>.bam
```

The output `bam` file should be sorted and indexed prior to the downstream analysis:

```bash
#Sort the output bam file by coordinate
samtools sort <sample>.bam -o <sample>_sorted.bam 

#Generate an index file
samtools index <sample>_sorted.bam
```

## Post-alignment QC

The post-alignment QC involves several steps:

- [Remove mitochondrial reads](#remove-mitochondrial-reads)
- [Remove duplicates & low-quality alignments](#remove-duplicates-&-low-quality-alignments) (including non-uniquely mapped reads)
- [Calculate library complexity and QC](#calculate-library-complexity-and-QC)
- [Remove ENCODE blacklist regions](#remove-encode-blacklist-regions)
- [Shift read coordinates](#shift-read-coordinates)

For an ATAC-seq experiment, the number of uniquely mapped reads ***after these steps*** is recommended to be 25 million for single-end or 50 million paired-end reads [(Buenrostro et al. 2015)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4374986/). Specific to ATAC-seq, an additional QC step is to check the fragment size distribution, which is expected to correspond to the length of nucleosomes:

- [Assess fragment size distribution](#assess-fragment-size-distribution)

### Remove mitochondrial reads

ATAC-seq experiments commonly include a high proportion of mitochondrial reads. These should be removed. To assess the total % of mitochondrial reads, `samtools idxstats` can be run to report the total number of reads mapped to each chromosome. `samtools flagstat` provides a short report including the total number of DNA reads as the first line (halve this number for the total number of fragments). 

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

#check the total number of aligned fragments
head <sample>_sorted.flagstat
```

The % of DNA fragments aligned to chrM can be calculated as a % of the total DNA fragments. To remove any mitocondrial DNA, run the following:

```bash
#Remove reads aligned to the mitochondria
samtools view -h <sample>_sorted.bam | grep -v chrM | samtools sort -O bam -o <sample>.rmChrM.bam -T .
```

### Mark duplicates 

To mark duplicate reads and view the % of duplicates:

```bash
picard MarkDuplicates QUIET=true INPUT=<sample>.rmChrM.bam OUTPUT=<sample>.marked.bam METRICS_FILE=<sample>.dup.metrics REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=.

#View the % of duplicates
head -n 8 <sample>.dup.metrics | cut -f 7,9 | grep -v ^# | tail -n 2
```

### Remove duplicates & low-quality alignments 

The output `sam/bam` files contain several measures of quality. First, the alignment quality score. Reads which are uniquely mapped are assigned a high alignment quality score and one genomic position. If reads can map to more than one location, Bowtie2 reports one position and assigns a low quality score. The proportion of uniquely mapped reads can be assessed. In general, >70% uniquely mapped reads is expected, while <50% may be a cause for concern [(Bailey et al. 2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3828144/pdf/pcbi.1003326.pdf). Secondly, the 'flag' reports information such as whether the read is mapped as a pair or is a PCR duplicate. The individual flags are reported [here](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/04_alignment_quality.html) and are combined in a `sam/bam` file to one score, which can be deconstructed back to the original flags using [online interpretation tools](https://broadinstitute.github.io/picard/explain-flags.html). In this pipeline, the bowtie2 parameters `--no-mixed` and `--no-discordant` prevent the mapping of only one read in a pair, so these flags will not be present. All flags reported in a `sam` file can optionally be viewed using  `grep -v ^@ <sample>.sam | cut -f 2 | sort | uniq`.

**Multi-mapping:** the user should decide whether or not to retain reads which have multi-mapped, i.e. aligned to more than one position in the reference genome. When using paired-end data, it may be the case that one read aligns to a repetitive region (and therefore can map elsewhere), while the mate aligns to a unique sequence with a high quality. The bowtie2 parameters used above required reads to align within 50-700bp, so there should be no reads incorrectly aligned outside this distance. As such, the user may decide to keep multi-mapping reads on the assumption that they are likely to be mapped to the correct sequence, within the length of the DNA fragment. This may, however, cause incorrect alignments in extended repetitive regions where a read could map to multiple positions within the length of the DNA fragment. This should be minimised by the downstream removal of the [ENCODE Blacklisted regions](https://www.nature.com/articles/s41598-019-45839-z).

If a read is multi-mapped, it is assigned a low quality score by bowtie2. To view how many DNA reads align with a quality score >30, run the following (divide this number by 2 to calculate the # of DNA fragments):

```bash
samtools view -q 30 -c <sample>.marked.bam
```

A low % of uniquely mapped reads map result from short reads, excessive PCR amplification or problems with the PCR (Bailey et al. 2013). The following code uses the `sam/bam` flags to retain properly mapped pairs (`-f 2`) and to remove reads which fail the platform/vendor QC checks (`-F 512`), duplicate reads (`-F 1024`) and those which are unmapped (`-F 12`). The three flags to be removed can be combined into `-F 1548`, which will remove reads which meet any of the three individual flags

Here we will ***remove*** multi-mapped reads:

```bash
samtools view -h -b -f 2 -F 1548 -q 30 <sample>.rmChrM.bam | samtools sort -o <sample>.filtered.bam

samtools index <sample>.filtered.bam

#To retain multi-mapped reads:
#samtools view -h -b -f 2 -F 1548 <sample>.rmChrM.bam | samtools sort -n -o <sample>.filtered.bam 
```

### Remove ENCODE blacklist regions

The [ENCODE blacklist regions](https://github.com/Boyle-Lab/Blacklist/), most recently reported by [Amemiya et al. (2019)](https://www.nature.com/articles/s41598-019-45839-z) are defined as 'a comprehensive set of regions in the human, mouse, worm, and fly genomes that have anomalous, unstructured, or high signal in next-generation sequencing experiments independent of cell line or experiment.' These problematic regions should be removed before further analysis. Download the blacklist files for your chosen reference genome from the [Boyle Lab github repository](https://github.com/Boyle-Lab/Blacklist/tree/master/lists). Details regarding the identification of blacklist regions are reported [here](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg19-blacklist-README.pdf).

```bash
#Remove reads within the blacklist regions
bedtools intersect -nonamecheck -v -abam <sample>.filtered.bam -b hg19-blacklist.v2.bed > <sample>.blacklist-filtered.bam

#Index the bam file
samtools index <sample>.blacklist-filtered.bam
```

### Shift read coordinates

An optional step in analysing data generated using the Tn5 transposase (such as ATAC-seq, ChIPmentation etc.) is to account for a small DNA insertion, introducted as repair of the transposase-induced nick introduces a 9bp insertion. Reads aligning to the + strand should be offset by +4bp and reads aligned to the -ve strand should be offset by -5bp. For references, see the first ATAC-seq paper by [Buenrostro et al., (2013)](https://www.nature.com/articles/nmeth.2688) and the analysis by [Adey et al., (2010)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-12-r119) which showed this insertion bias. Shifting coordinates is only really important if single-base resolution is required, for example in the analysis of transcription factor motifs in ATAC-seq peak footprints. Be aware that some tools do this shifting themselves (so double check manuals!).
 

We can use the `deeptools` command  `alignmentSieve`.


```bash
#The user can set the preferred number of processors 
alignmentSieve --numberOfProcessors max --ATACshift --blackListFileName hg19-blacklist.v2.bed --bam <sample>.blacklist-filtered.bam -o <sample>.tmp.bam

#Sort and index the bam file
#Set the number of preferred threads with the -@ option
samtools sort -O bam -o <sample>.shifted.bam <sample>.tmp.bam
samtools index <sample>.shifted.bam

rm <sample>.tmp.bam
```

### Assess fragment size distribution and QC

The fragment size is expected to show a periodicity of 150/200 bp, reflecting the length of DNA surrounding nucleosomes, since the tagmentation typically cuts DNA in between nucleosomes. The `R` tool [ATACseqQC](https://www.bioconductor.org/packages/release/bioc/html/ATACseqQC.html) can be used to assess the distribution of fragments. Fragments sizes may also be <100bp, corresponding to nucleosome-free regions or cutting of linker DNA, as well as fragments of and mono-, di-, and tri-nucleosomes (~200, 400 and 600bp, respectively) [(Yan et al. 2020)](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-1929-3). Fragments from nucleosome-free regions are expected to be enriched around transcription-start sites (TSS) and fragments from nucleosome-bound regions are depleted around TSS and maybe be slightly enriched in the flanking regions. The ATACseqQC manual is available [here](https://www.bioconductor.org/packages/release/bioc/vignettes/ATACseqQC/inst/doc/ATACseqQC.html).

- Estimate library complexity
- Fragment size distribution 
- GC bias
- Nucleosome positioning 
- Plot footprints
- Plot correlations between samples

The <sample>.shifted.bam file should be analysed in the following steps.


**TO BE COMPLETED**

## Bam visualisation

Through this pipeline, two types of tracks will be generated for visualisation in genome browsers. The first, generated here, will show the aligned reads and are generated from the processed `bam` file. The second, generated after peak calling, will show the -log<sub>10</sub> p-value from the peak calling.

The `deeptools` command `bamCoverage` will be used. ATAC-seq data is commonly normalised to the effective genome size...

**TO BE COMPLETED**

```bash
# bam to bigwig
# Set your preferred number of processors
bamCoverage --numberOfProcessors 8 --binSize 10 --normalizeUsing RPGC \
  --effectiveGenomeSize $effect_genome_size --bam <sample>.shifted.bam -o <sample>.shifted.bw
```	

## Peak calling

Peaks are identified where sequenced reads accumulate. These correspond to regions of accessible DNA.

In this pipeline, peaks will be called using several alternative software: [MACS2](https://github.com/macs3-project/MACS), [HMMRATAC, developed by Tarbell et al. (2019)](https://academic.oup.com/nar/article/47/16/e91/5519166) and the popular but unpublished [Genrich](https://github.com/jsh58/Genrich). While MACS2 is the most popular peak caller and may be used to generate results consistent with other published data, [Yan et al. (2020)](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-1929-3) recommend using the peak caller HMMRATAC, developed by [Tarbell et al. (2019)](https://academic.oup.com/nar/article/47/16/e91/5519166), which is specifically developed for ATAC-seq data (if computational resources are sufficient). HMMRATAC is available on [github](https://github.com/LiuLabUB/HMMRATAC). Genrich also has a dedicated ATAC-seq mode and is particularly popular as many QC steps are built in. First, MACS2 will be used to call peaks.

Important considerations for ATAC-seq: 
- There are usually **no controls**
- The Tn5 transposase has a binding preference, resulting in a GC bias which should be corrected for during peak calling
- Repair of the transposase-induced nick introduces a 9bp insertion (completed above). 

### Peak calling - MACS2 

The following steps will be carried out:

1. Call peaks for individual replicates
2. Call replicated peaks for the pooled replicates
3. Call 'high-confidence' IDR peaks.

If using paired-end reads, MACS2 will be used with the `-f BAMPE` option, so that MACS2 calculates the pileup for full fragments. Depending on the analysis aims, there are several different options that can be used. The ENCODE3 pipeline uses the `--nomodel --shift -37 --extsize 73` options for analysing ATAC-seq data, to account for the size of nucleosomes. Nucleosomes cover \~145 bp and the ATAC-seq reads need to be shifted towards the 5' end by half this distance.

```bash 
#call peakrs
macs2 callpeak -f BAMPE --nomodel --shift -37 --extsize 73 -g hs --keep-dup all --cutoff-analysis -n <sample> -t <sample>.shifted.bam --outdir macs2/<sample> 2> macs2.log
```

The output files:

- `_broad_treat_pileup.bdg`: bedGraph format for the treatment sample
- `_broad_control_lambda.bdg`: bedGraph format for the background note
- `_broad_peaks.broadPeak`: a BED6+4 file detailing the peak locations, along with the peak summits, *p*-value and *q*-values 
- `_broad_peaks.gappedPeak`
- `_broad_peaks.xls`: a tabular file containing addition information, such as pileup and fold-enrichment.


![#1589F0](https://via.placeholder.com/15/1589F0/000000?text=+) **Output file**: the `_broad_peaks.broadPeak` is a `bed` file containing the peak information for the INDIVIDUAL replicate\*.

\*The `<sample>_peaks.narrowPeak` can be uploaded and visualised via a genome browser such as UCSC. The `bed` file of peak calls is referred to at this stage as 'relaxed' peak calls, since they are called for individual replicates. Two or more biological replicates will be combined in the next stage to generate a combined set of peaks.


The total number of peaks can be obtained using `wc -l <sample>_peaks.narrowPeak`. 

![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) **QC value**: input the total number of peaks into the QC spreadsheet.


From these output files, we will generate:

1. A `bigWig` track of the fold-enrichment (treatment over the background)
2. A `bigWig` track of the -log<sub>10</sub> *p*-value (treatment over the background)


**1. Fold-enrichment bigWig**

The following commands require an input file detailing the chromosome sizes. Use the UCSC tool `fetchChromSizes` (install via [conda](https://anaconda.org/bioconda/ucsc-fetchchromsizes)): `fetchChromSizes hg38 > hg38.chrom.sizes`

```bash
#Generate the fold-change bedGraph
macs2 bdgcmp -t <sample>.broad_treat_pileup.bdg -c <sample>.broad_control_lambda.bdg -m FE -o <sample>_FE.bdg 

#Sort the bedGraph file and convert to bigWig
sort -k1,1 -k2,2n <sample>_FE.bdg > <sample>_FE.sorted.bdg

bedGraphToBigWig <sample>_FE.sorted.bdg hg38.chrom.sizes > <sample>_macs2_FE.bw
```

![#1589F0](https://via.placeholder.com/15/1589F0/000000?text=+) **Output file**: `<sample>_macs2_FE.bw`

**2. -log<sub>10</sub> *p*-value bigWig**

```bash
#Generate the p-value bedGraph
macs2 bdgcmp -t <sample>.broad_treat_pileup.bdg -c <sample>.broad_control_lambda.bdg -m ppois --o <sample>_ppois.bdg

#Sort the bedGraph file and convert to bigWig
sort -k1,1 -k2,2n <sample>_ppois.bdg > <sample>_ppois.sorted.bdg

bedGraphToBigWig <sample>_ppois.sorted.bdg hg38.chrom.sizes > <sample>_macs2_pval.bw
```

![#1589F0](https://via.placeholder.com/15/1589F0/000000?text=+) **Output file**: `<sample>_macs2_pval.bw`

The `<sample>_macs2_pval.bw` and `<sample>_macs2_FE.bw` output files can visualised in a genome browser, such as UCSC.


### Call peaks for pooled replicates

The step assumes that the ChIP-seq expriment includes *biological replicates* for each treated condition. Best practise requires a combined set of peaks for the pooled replicates to be called. Assuming there are two biological replicates, `rep1` and `rep2`:

First, check the correlation between the replicates using the UCSC tool wigCorrelate:

```bash
wigCorrelate <sample>_rep1_macs2_FE.bw <sample>_rep2_macs2_FE.bw
```
Assuming there is a staisfactory correlation, call peaks on the combined replicates by including all the files in the `macs2 callpeak` command:


```bash 
#call peakrs
macs2 callpeak -f BAMPE --nomodel --shift -37 --extsize 73 -g hs --keep-dup all --cutoff-analysis -n <sample> -t <sample>_rep1.shifted.bam <sample>_rep2.shifted.bam --outdir macs2/<sample> 2> macs2.log
```

The number of peaks within a replicated peak file should be >150,000, though values >100,000 may be acceptable. 
The number of peaks within an IDR peak file should be >70,000, though values >50,000 may be acceptable.
A nucleosome free region (NFR) must be present.
A mononucleosome peak must be present in the fragment length distribution. These are reads that span a single nucleosome, so they are longer than 147 bp but shorter than 147*2 bp. Good ATAC-seq datasets have reads that span nucleosomes (which allows for calling nucleosome positions in addition to open regions of chromatin).


1. **Call peaks**: MACS2
2. **Visualise**: generate -log<sub>10</sub> p-value bigwig tracks for MACS2 peaks
3. **Quality control**
4. Repeat above, calling peaks with HMMRATAC and Genrich


The number of peaks within a replicated peak file should be >150,000, though values >100,000 may be acceptable. The number of peaks within an IDR peak file should be >70,000, though values >50,000 may be acceptable.


### HMMRATAC

```bash
samtools view -H <sample>.filtered.bam | perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > genome.info

java -jar HMMRATAC_V1.2.10_exe.jar -b <sample>.filtered.bam -i <sample>.filtered.bam.bai -g genome.info -o <sample>
```

### Genrich 

See this [post](https://informatics.fas.harvard.edu/atac-seq-guidelines.html#peak) from Harvard FAS Informatics.
Genrich does everything in one go!


### Other peak callers

[NucleoATAC](https://nucleoatac.readthedocs.io/en/latest/) ? 

... check reproducibility of peaks between replicates... then re-run MACS2 with the merged bam file.


#### Peak QC

Quality control steps should be carried out to assess the called peaks, as well as reproducibility between samples.

Quality control of the peaks, along with differential accessiblity analysis (if this is an aim of your project) will be carried out. 

***The following analysis will be completed in R***

See [ENCODE ATAC-seq data standards and prototype processing pipeline](https://www.encodeproject.org/atac-seq/)

- The number of peaks
- Fraction of reads in peaks (FRiP) score 
- Transcription start site (TSS) enrichment 
- Biological replicate peak overlap
- Compare normalisation methods and systematic biases

Number of peaks should be >150,000 and not less than 100,000 (>70,000 in an IDR file)


## Differential accessibility 

[Yan et al. (2020)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1929-3) recommendations for ATAC-seq data analysis:
csaw for peak differential analysis
ChIPseeker for annotation and visualisation
MEME suite for motic detection and enrichment
HMMRATAC for nucleosome detection
HINT-ATAC for footprint analysis
PCEA for regulatory network reconstruction with RNA-seq
HOMER for motif discovery.


TSS enrichment should be ideally >10 for hg19 and not <6 and ideally >7 and not <5 for GRCh38.

HOMER can be used for genomic annotation.

FRiP score should be >0.3 and not less than 0.2. 

Biological replicate peak overlap: ENCODE-define naive overlap which calls peaks on pooled replicates and then identifies peaks with >50% overlap with all single replicate peaks.


## Visualisation

The QC-ed `bam` file can be converted to a `bedGraph` format to visualise sequencing trakcs using tools such as the UCSC browser or the integrative genomes browser. The ENCODE blacklist regions can be provided, to exclude them from the output:

```bash
bedtools bamCoverage --normalizeUsing BPM -b <sample>.filtered.bam > <sample>.bedGraph
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

## Resources

This pipeline has been developed using guidance from the following resources:

- https://vallierlab.wixsite.com/pipelines/atac-seq A great tool for beginners indicating the major analysis steps 
- https://www.encodeproject.org/atac-seq/ The recommended ENCODE pipeline, for which tools are available on [github](https://github.com/ENCODE-DCC/atac-seq-pipeline) and the recommended parameters/specification are available via a [google doc](https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit)
- https://github.com/harvardinformatics/ATAC-seq An ATAC-seq pipeline from Harvard Informatics 
- https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/ A similar github page presenting an ATAC-seq pipeline 
- https://galaxyproject.github.io/training-material/topics/epigenetics/tutorials/atac-seq/tutorial.html
- Standardized ENCODE pipeline devides by Kundaje et al. https://libraries.io/github/kundajelab/atac_dnase_pipelines
- [Hbctaining tutorial for peak calling](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html)

An recent review on the ATAC-seq analysis pipeline is reported by [Yan et al. (2020)](https://genomebiology.biomedcentral.com/track/pdf/10.1186/s13059-020-1929-3).