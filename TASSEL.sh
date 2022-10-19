#!/bin/bash
# user input required
hisat_indices="/path/to/hisatindices"
rnastranded=R
referenceGTF="/path/to/reference_annotation_file"
referenceFASTA="/path/to/reference_fasta_file"
shortread_fastq_dir="/path/to/directory_with_only_short-read_fastqfiles"
longread_fastq_dir="/path/to/directory_with_only_long-read_fastqfiles"
library_type=--rf
primer1=GCTCTATCTTCTTT
primer2=CTGATATTGCTGGG
rc_primer2=CCCAGCAATATCAG
processor=4

# Aligning short-reads
mkdir short-read_bam
for shortfastq in ${shortread_fastq_dir}/*.fastq
do
  base=`basename $shortfastq '.fastq'`
  echo "Aligning reads from" $base'.fastq'
  hisat2 -x $hisat_indices --rna-strandness $rnastranded --dta -p $processor -U ${shortread_fastq_dir}/${base}.fastq | samtools view -bS - > ./short-read_bam/${base}.bam
  samtools sort ./short-read_bam/${base}.bam -o ./short-read_bam/${base}.sorted.bam -@$processor
done
# assembling short-read transcripts from sorted bam files
mkdir gtf_dir
for shortbam in short-read_bam/*.sorted.bam
do
  base=`basename $shortbam '.sorted.bam'`
  echo "Assembling transcripts from" $base'.sorted.bam'
  stringtie -p $processor -G $referenceGTF $library_type -o ./gtf_dir/${base}.gtf ./short-read_bam/${base}.sorted.bam
done

# Stranding long reads
mkdir long-read_strandedfastq
for longfastq in ${longread_fastq_dir}/*.fastq
do
base=`basename $longfastq '.fastq'`
temp="$longfastq"
basename=${temp%.*}
# calculate number of reads in the parent file
parent_reads=$(wc -l ${temp} | awk '{print $1 * 0.25}')
echo $parent_reads "reads found in" $temp
# extract reads that are minimum 100 nt in length
seqkit seq -g -m 100 ${temp} > "$basename"_minimum100bp.fq 
# extract reads with q1 sequence in the first 100 bp of positive strand and allowing 2 mismatches 
cat "$basename"_minimum100bp.fq | seqkit grep -s -i -P -R 1:100 -m 2 -p $primer1 -o "$basename"_q1m2.fq
# extract reads with q3 sequence in the first 100 bp of positive strand and allowing 2 mismatches 
cat "$basename"_minimum100bp.fq | seqkit grep -s -i -P -R 1:100 -m 2 -p $primer2 -o "$basename"_q3m2.fq
# extract reads with q3rc sequence in the last 100 bp of positive strand and allowing 2 mismatches 
cat "$basename"_minimum100bp.fq | seqkit grep -s -i -P -R -100:-1 -m 2 -p $rc_primer2 -o "$basename"_q3m2rc.fq
# making reverse complement of q3m2
seqkit seq  -t dna -r -p "$basename"_q3m2.fq > "$basename"_rcq3m2.fq
# combining q1m2 and q3m2rc and removing duplicates
cat "$basename"_q1m2.fq "$basename"_q3m2rc.fq | seqkit rmdup -o "$basename"_q1m2_q3m2rc.fq
# combining q1m2_q3m2rc and rcq3m2 and removing duplicates
cat "$basename"_q1m2_q3m2rc.fq "$basename"_rcq3m2.fq > "$basename"_q1m2_q3m2rc_rcq3m2.fq 
seqkit rmdup "$basename"_q1m2_q3m2rc_rcq3m2.fq -o long-read_strandedfastq/"$base"_SLURP_stranded.fastq
wc -l long-read_strandedfastq/"$base"_SLURP_stranded.fastq | awk '{print $1 * 0.25, "reads found in '$base'_SLURP_stranded.fastq", $1*0.25/'$parent_reads'*100, "% of parent"}'
rm ${longread_fastq_dir}/*.fq
done

# alignment of stranded long reads
mkdir long-read_bam
# Aligning long-reads from fastq files
for longfastq in long-read_strandedfastq/*.fastq
do
  base=`basename $longfastq '.fastq'`
  echo "Aligning reads from" $base'.fastq'
  minimap2 -ax splice -un -t $processor --MD $referenceFASTA long-read_strandedfastq/${base}.fastq | samtools view -bS - > ./long-read_bam/${base}.bam
  samtools sort ./long-read_bam/${base}.bam -o ./long-read_bam/${base}.sorted.bam -@$processor
done

# assembling long-read transcripts from sorted bam files
for longbam in long-read_bam/*.sorted.bam
do
  base=`basename $longbam '.sorted.bam'`
  echo "Assembling transcripts from" $base'.sorted.bam'
  stringtie -p $processor -G $referenceGTF -L -o ./gtf_dir/${base}.gtf ./long-read_bam/${base}.sorted.bam
done
find gtf_dir -type f > mergelist.txt
echo "Creating TASSEL-merged gtf..."
stringtie --merge --rf -c 0 -F 0 -T 0 -f 0 mergelist.txt -o TASSEL-merged.gtf
echo "TASSEL completed!"
