#!/bin/bash

# Usage : ./launch.sh GM12878_RNAPII GM12878_RNAPII_R1.fastq GM12878_RNAPII_R2.fastq

# Parse arguments #############################################################

# Use -gt 0 to consume one or more arguments per pass in the loop (e.g.
# some arguments don't have a corresponding value to go with it such
# as in the --default example).
# note: if this is set to -gt 0 the /etc/hosts part is not recognized ( may be a bug )
while [[ $# -gt 1 ]]
do
    key="$1"
    
    case $key in
        -f|--foward)
        fastq1="$2"
        shift # past argument
        ;;
        -r|--reverse)
        fastq2="$2"
        shift # past argument
        ;;
        -b|--bam)
        bam="$2"
        shift # past argument
        ;;
        -o|--output)
        prefix="$2"
        shift # past argument
        ;;
        -s|--species)
        species="$2"
        shift # past argument
        ;;    
        *)
                # unknown option
        ;;
    esac
    shift # past argument or value
done

case $species in
    hg19)
    bowtieindex=/is1/commonDatasets/mugqic_space/genomes/species/Homo_sapiens.hg19/Sequence/BowtieIndex/genome
    genomefasta=/is1/commonDatasets/mugqic_space/genomes/species/Homo_sapiens.hg19/Sequence/WholeGenomeFasta/genome.fa
    chrMName=chrM
    ;;
    hg38)
    ;;
    mm9)
    ;;
    mm10)
    ;;
    *)
    ;;
esac


module load R/R-3.3.0 python/2.7.8 parallel/20140922

mkdir -p $prefix/fastq
if [[ "$bam" != "" ]]
then
    # Find out the sample identifier.
    id=$(basename $bam .bam)

    # Split the bam file into two FASTQ files.
#    java -jar $PICARD_PATH SamToFastq I=$bam F=$prefix/fastq/$id.R1.fastq.gz F2=$prefix/fastq/$id.R2.fastq.gz
    
    # Set the input fastq file names.
    fastq1=$prefix/fastq/$id.R1.fastq.gz
    fastq2=$prefix/fastq/$id.R2.fastq.gz
fi

gzregex="(.*\/)*(.*)(\.fq|\.fastq)(\.gz)*$"
if [[ $fastq1 =~ $gzregex ]]
then
    id1=${BASH_REMATCH[2]}
else
    echo "Could not detect fastq1 file type. Aborting."
fi

if [[ $fastq2 =~ $gzregex ]]
then
    id2=${BASH_REMATCH[2]}
else
    echo "Could not detect fastq2 file type. Aborting."
fi

mkdir -p $prefix/cut
#~/.local/bin/cutadapt -a ACGCGATATCTTATCTGACT -A AGTCAGATAAGATATCGCGT \
#                      --minimum-length 17 --overlap 10 -o $prefix/cut/$id1.fastq -p $prefix/cut/$id2.fastq \
#                      $fastq1 $fastq2

rm ${prefix}/mango/mango_1.same.fastq
rm ${prefix}/mango/mango_2.same.fastq
ln -s `pwd`/$prefix/cut/$id1.fastq ${prefix}/mango/mango_1.same.fastq
ln -s `pwd`/$prefix/cut/$id2.fastq ${prefix}/mango/mango_2.same.fastq
                      
                      
mkdir -p $prefix/mango
Rscript mango/mango/mango.R --bedtoolspath bedtools/bedtools2/bin/bedtools \
                            --bowtiepath bowtie/bowtie-1.1.2/bowtie \
                            --fastq1 $prefix/cut/$id1.fastq \
                            --fastq2 $prefix/cut/$id2.fastq \
                            --outdir $prefix/mango \
                            --bowtieref $bowtieindex \
                            --bedtoolsgenome $genomefasta \
                            --chromexclude $chrMName,chrY \
                            --stages 2:5 \
                            --shortreads FALSE \
                            --peakslop 1500 \
                            --reportallpairs TRUE
                      
#for fileR1 in $(ls *_R1.[a-z][a-z].linkerONLYreads.fastq)
#do
#	fileR2=$(echo $fileR1 | sed 's/_R1\./_R2\./g')
#	tmpid1=$(basename $fileR1 .fastq)
#	tmpid2=$(basename $fileR2 .fastq)
#
#	#python fastqCombinePairedEnd.py $fileR1 $fileR2
#
#	#cutadapt -a ACGCGATATCTTATCTGACT -A AGTCAGATAAGATATCGCGT --minimum-length 17 --overlap 10 -o $tmpid1.cutadapt.fastq -p $tmpid2.cutadapt.fastq ${fileR1}_pairs_R1.fastq ${fileR2}_pairs_R2.fastq
#done
#
##cat ${prefix}_R1.*.cutadapt.fastq > ${prefix}_1.same.fastq
##cat ${prefix}_R2.*.cutadapt.fastq > ${prefix}_2.same.fastq
#
##Rscript mango/mango/mango.R --fastq1 ${prefix}_1.same.fastq --fastq2 ${prefix}_2.same.fastq --prefix $prefix --bowtieref /is1/commonDatasets/mugqic_space/genomes/species/Homo_sapiens.hg19/Sequence/BowtieIndex/genome --bedtoolsgenome /is1/commonDatasets/mugqic_space/genomes/species/Homo_sapiens.hg19/Sequence/WholeGenomeFasta/genome.fa --chromexclude chrM,chrY --stages 2:5 --shortreads FALSE --peakslop 1500 --reportallpairs TRUE
#Rscript mango/mango/mango.R --fastq1 ${prefix}_1.same.fastq --fastq2 ${prefix}_2.same.fastq --prefix $prefix --bowtieref /is1/commonDatasets/mugqic_space/genomes/species/Homo_sapiens.hg19/Sequence/BowtieIndex/genome --bedtoolsgenome bedtools/hg19.genome --chrominclude chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22 --stages 4:5 --shortreads FALSE --peakslop 1500 --reportallpairs TRUE --keepempty TRUE --maxlength 1000


