#!/bin/bash

# Usage : ./launch.sh GM12878_RNAPII GM12878_RNAPII_R1.fastq GM12878_RNAPII_R2.fastq

module load R/R-3.1.0 python/2.7.8 parallel/20140922

prefix=$1
fastq1=$2
fastq2=$3

id1=$(basename $fastq1 .fastq)
id2=$(basename $fastq2 .fastq)

#split -l 400000000 $fastq1 ${id1}.
#split -l 400000000 $fastq2 ${id2}.

#for file in $(ls ${id1}.[a-z][a-z]; ls ${id2}.[a-z][a-z])
#do
	#grep "CGCGATATCTTATCTGACT\|GTCAGATAAGATATCGCGT" -B 1 -A 2 --no-group-separator $file > ${file}.linkerONLYreads.fastq
#done

for fileR1 in $(ls *_R1.[a-z][a-z].linkerONLYreads.fastq)
do
	fileR2=$(echo $fileR1 | sed 's/_R1\./_R2\./g')
	tmpid1=$(basename $fileR1 .fastq)
	tmpid2=$(basename $fileR2 .fastq)

	#python fastqCombinePairedEnd.py $fileR1 $fileR2

	#cutadapt -a ACGCGATATCTTATCTGACT -A AGTCAGATAAGATATCGCGT --minimum-length 17 --overlap 10 -o $tmpid1.cutadapt.fastq -p $tmpid2.cutadapt.fastq ${fileR1}_pairs_R1.fastq ${fileR2}_pairs_R2.fastq
done

#cat ${prefix}_R1.*.cutadapt.fastq > ${prefix}_1.same.fastq
#cat ${prefix}_R2.*.cutadapt.fastq > ${prefix}_2.same.fastq

#Rscript mango/mango/mango.R --fastq1 ${prefix}_1.same.fastq --fastq2 ${prefix}_2.same.fastq --prefix $prefix --bowtieref /is1/commonDatasets/mugqic_space/genomes/species/Homo_sapiens.hg19/Sequence/BowtieIndex/genome --bedtoolsgenome /is1/commonDatasets/mugqic_space/genomes/species/Homo_sapiens.hg19/Sequence/WholeGenomeFasta/genome.fa --chromexclude chrM,chrY --stages 2:5 --shortreads FALSE --peakslop 1500 --reportallpairs TRUE
Rscript mango/mango/mango.R --fastq1 ${prefix}_1.same.fastq --fastq2 ${prefix}_2.same.fastq --prefix $prefix --bowtieref /is1/commonDatasets/mugqic_space/genomes/species/Homo_sapiens.hg19/Sequence/BowtieIndex/genome --bedtoolsgenome bedtools/hg19.genome --chrominclude chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 --stages 4:5 --shortreads FALSE --peakslop 1500 --reportallpairs TRUE --keepempty TRUE --maxlength 1000


