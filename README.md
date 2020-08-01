# ezRAD-Pool-seq

## An implementation of the ezRAD pipeline described in Phair _et. al._ (2019).

#### Getting the publication raw sequence data
The sequence data for this study has deposited into the NCBI Sequence read archive (SRA) and may be accessed from the [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/) using the [PRJNA503110](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA503110) Bioproject accession number. Sequence read metadata and sequence read accession is available for download. We use the list of accessions to download the raw sequence read data using the [SRA toolkit](https://www.ncbi.nlm.nih.gov/books/NBK158900/#SRA_download.how_do_i_download_and_insta).

To downwload the sequences reads, this may take a while as at all sequences are ~7.5 GB.
``` 
xargs -a  SraAccList.txt -I{} prefetch {}

```
To generate paired-ends of the downloaded sequence data.

```
xargs -a ../SraAccList.txt -I{} fastq-dump –split-e  {}
 
```


#### Downloading the _Z. marina_ genome from NCBI Genomes. 

The [_Zostera marina_](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=29655&lvl=3&lin=f&keep=1&srchmode=1&unlock) assembly is available by [ftp](https://ftp.ncbi.nih.gov/genomes/genbank/plant/Zostera_marina/all_assembly_versions/GCA_001185155.1_Zosma_marina.v.2.1/). The [genomic.fna.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/185/155/GCA_001185155.1_Zosma_marina.v.2.1/GCA_001185155.1_Zosma_marina.v.2.1_genomic.fna.gz). Interesting piece on contigs and scaffolds [here](https://www.pacb.com/blog/genomes-vs-gennnnes-difference-contigs-scaffolds-genome-assemblies/).


#### Indexing the reference genome.

Create an index of the reference genome using the the BWA algorithm (https://github.com/lh3/bwa)

```
  bwa index GCA_001185155.1_Zosma_marina.v.2.1_genomic.fna

```
This index is used in mapping of reads to reference step.


Create an index of the reference for the ```samtools mpileup``` step.



#### Mapping reads to the reference 

The gzipped FastQ reads are mapped to the reference genome generated by the BWA algorithm to generate SAM files with minimum score of 20 parameter. The location of the BWA index is set in the _params.bwa_ref_ parameter of the nextflow file. SAM files are converted to BAM files using [SAMtools](http://www.htslib.org/doc/samtools.html) ([samtools cheatsheet](https://coding4medicine.com/cheatsheets/samtools.html)).
Mapped reads and unmapped reads were calculated using the samtools _idxstats_ option which reports the summary statistics. There is a nice explanation with visualisations on sequencing depth and coverage from [ecSeq](https://www.ecseq.com/support/ngs/how-to-calculate-the-coverage-for-a-sequencing-experiment). Sorted bam files are then used generate pileup file using the ```mpileup```  command to samtools.


### SNP calling with Popoolation2
Pileup files are then analysed using [Popoolation2](https://sourceforge.net/projects/popoolation2/files/latest/download) which is well suited for pool data such as ezRAD. The pileup file is converted into a sync file and SNPs called. SNPs  are called from multiple regions, however, the SNP to Genepop tool can only be used on a single regions. As SNPs are callled from multiple regions, a custom Python script ```get_region.py``` is used to extract regions and coordinates based on the *_rc file. This script is such that a single contig and coordinates coveraging all SNPs are covered and all regions meeting the target and maximum coverage parameters.

Multiple Genepop files are generated with each contig represented by a Genepop file. These file need to be merged to make a master Genepop file and the  custom Python```merge_genepop.py``` is used.

## Citations ##

Phair NL, Toonen RJ, Knapp I, von der Heyden S. 2019. Shared genomic outliers across two divergent population clusters of a highly threatened seagrass. PeerJ 7:e6806 http://doi.org/10.7717/peerj.6806
