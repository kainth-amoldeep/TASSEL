# TASSEL (Transcript Assembly using Short and Strand Emended Long reads)
<img align="left" width="350" src="https://user-images.githubusercontent.com/66103719/196807655-e1bc74a1-cf67-47eb-ad26-90985af0fbae.png">
TASSEL is a hybrid transcript assembly pipeline that merges transcriptome from short-read RNA-seq and long-read RNA-seq. The output is a merged transcritome file (gtf) which combines high depth of short-read sequencing with long-range information from long-read RNA-seq. The unique feature about TASSEL is that it strands the otherwise unstranded long reads using inbuilt *SLURP* methodology and then use them for transcript assembly. 

### Usage
Create two directories - one that contains short-read fastq files and another that contains long-read fastq files that you want to merge.<br/>
<br/>
**User input is required for the following variables within the script:**<br/>
**hisat_indices:** provide path to hisat indices<br/>
**rnastranded:** provided type of strandedness for short-reads (options: F, R; default: R)<br/>
**referenceGTF:** provide path to reference_annotation_file (required for guided assembly)<br/>
**referenceFASTA:** provide path to reference_fasta_file<br/>
**shortread_fastq_dir:**  provide path to directory that contains only short-read fastq files<br/>
**longread_fastq_dir:** provide path to directory that contains only short-read fastq files<br/>
**library_type:** for library strandedness (options: --rf (for first strand), --fr (for second strand); default --rf)<br/>
**primer1:** sequence of primer used for first strand synthesis during long-read cDNA library prep; limit to 15 nt (default: GCTCTATCTTCTTT)<br/>
**primer2:** sequence of primer used for second strand synthesis (strand switching) during long-read cDNA library prep; limit to 15 nt (default: CTGATATTGCTGGG)<br/>
**rc_primer2:** reverse complement of primer2 (default: CCCAGCAATATCAG)<br/>
**processor:** number of processors (default: 4)<br/>

Run ```bash TASSEL.sh``` in the directory that contains directories for short-read and long-read fastq files.<br/>

#### Dependencies
**hisat2:** hisat2 can be obtained from (http://daehwankimlab.github.io/hisat2/download/) or ```conda install -c bioconda hisat2```<br/>
**samtools:** samtools can be obtained from (http://www.htslib.org/download/) or ```conda install -c bioconda samtools```<br/>
**StringTie2:** StringTie2 can be obtained from (https://ccb.jhu.edu/software/stringtie/index.shtml) or ```conda install -c bioconda stringtie```<br/> 
**seqkit:** seqkit can be obtained from (https://bioinf.shenwei.me/seqkit/) or ```conda install -c bioconda seqkit```<br/>
**minimap2:** seqkit can be obtained from (https://github.com/lh3/minimap2#install) or ```conda install -c bioconda minimap2```<br/>
