#!/bin/bash

# Usage : ./launch.sh GM12878_RNAPII GM12878_RNAPII_R1.fastq GM12878_RNAPII_R2.fastq

# Parse arguments #############################################################

# Use -gt 0 to consume one or more arguments per pass in the loop (e.g.
# some arguments don't have a corresponding value to go with it such
# as in the --default example).
# note: if this is set to -gt 0 the /etc/hosts part is not recognized ( may be a bug )
while [[ $# -gt 0 ]]
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
        -l|--forcelinker)
        forcelinker="Y"
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
    genomefasta=/is1/commonDatasets/mugqic_space/genomes/species/Homo_sapiens.hg19/Sequence/WholeGenomeFasta/genome.fa.fai
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


module load R/R-3.3.0_BioC-3.3 python/2.7.8 parallel/20140922

# Convert the BAM file into fastq files #################################################

# If the input is a BAM file, convert it to two fastq files.
mkdir -p $prefix/fastq
if [[ "$bam" != "" ]]
then
    # Find out the sample identifier.
    id=$(basename $bam .bam)

    # Determine the output file names.
    fastq1=$prefix/fastq/$id.R1.fastq.gz
    fastq2=$prefix/fastq/$id.R2.fastq.gz

    # If the output file already exist, do not regenerate them.
    if [ ! -e $fastq1 ] || [ ! -e $fastq2 ]
    then
        # Split the bam file into two FASTQ files.    
        echo "Extracting fastq from bam file."
        java -jar $PICARD_PATH SamToFastq I=$bam \
           F=$prefix/fastq/$id.R1.fastq.gz \
           F2=$prefix/fastq/$id.R2.fastq.gz &
    fi
fi
wait

# Figure out the fastq file basenames. ########################################
gzregex="(.*\/)*(.*)(\.fq|\.fastq)(\.gz)*$"
if [[ $fastq1 =~ $gzregex ]]
then
    id1=${BASH_REMATCH[2]}
else
    echo "Could not detect fastq1 file type. Aborting."
    exit
fi

if [[ $fastq2 =~ $gzregex ]]
then
    id2=${BASH_REMATCH[2]}
else
    echo "Could not detect fastq2 file type. Aborting."
    exit
fi


# Remove linker sequences #####################################################
mkdir -p $prefix/cut

cut1=$prefix/cut/$id1.fastq.gz
cut2=$prefix/cut/$id2.fastq.gz
if [ ! -e $cut1 ] || [ ! -e $cut2 ]
then
    echo "Removing linker sequences.."
    
    # Are we keeping only trimmed reads?
    extraoptions=""
    if [[ "$forcelinker" == "Y" ]]
    then
        extraoptions="--discard-untrimmed"
    fi
    
    ~/.local/bin/cutadapt -a ACGCGATATCTTATCTGACT -A AGTCAGATAAGATATCGCGT \
                          --minimum-length 17 --overlap 10 $extraoptions  --info-file $prefix/cut/cutadapt.infofile \
                          -o $cut1 -p $cut2 $fastq1 $fastq2 1> $prefix/cut/cutadapt.stdout 2> $prefix/cut/cutadapt.stderr
fi


# Perform alignment #####################################################

# Mango can generate alignments itself, but it is limited to 1 thread, which
# makes the process painfully slow. We'll do it outside of mango, using the same
# parameters.
# Mango uses the mapq 40 bowtie parameter to set the MAPQ of uniquely mapped reads to
# a specific value. This option is not present in alternate aligners, some
# we're stuck with bowtie1 unless we modify mango itself.
mkdir -p $prefix/alignment
align1=$prefix/alignment/$id1.sam
align2=$prefix/alignment/$id2.sam
for alignout in $align1 $align2
do
    if [ ! -e $alignout.unsorted ]
    then
        fastqid=`basename $alignout .sam`

        gzip -dc $prefix/cut/$fastqid.fastq.gz |
            bowtie/bowtie-1.1.2/bowtie -S -n 2 -l 50 -k 1 --chunkmbs 500 --sam-nohead \
                                   --mapq 40 -m 1 --best --threads 8 --phred33-quals \
                                   $bowtieindex \
                                   - \
                                   $alignout.unsorted &> $alignout.log &
    fi
done
wait

for alignout in $align1 $align2
do
    if [ ! -e $alignout ]
    then
        # Mango requires sam files to have reads in the exact same order, and to have the 
        # exact same name. Bowtie somehow mixes up the order (probably because of multithreading).
        # So we need to resort by query name. To do that and get the right results,
        # we'll first need to remove the read markers (/1, /2) at the end of the read names.
        sed -i -e 's/\/[12]//' $alignout.unsorted &
    fi
done
wait

for alignout in $align1 $align2
do
    if [ ! -e $alignout ]
    then
        # We'll also need to sort numerically on the very last part of the read name,
        # because entries ending with, for example, 27892 and 2789 have no guaranteed sort
        # order under the default parameters.
        sort -t: -k 1,6d -k 7,7n $alignout.unsorted > $alignout &
    fi
done
wait




# Run mango #########################################################
mkdir -p $prefix/mango

# Generate symbolic links for the sam files.
ln -s `pwd`/$align1 ${prefix}/mango/mango_1.same.sam
ln -s `pwd`/$align2 ${prefix}/mango/mango_2.same.sam

Rscript mango/mango/mango.R --bedtoolspath bedtools/bedtools2/bin/bedtools \
                            --bowtiepath bowtie/bowtie-1.1.2/bowtie \
                            --fastq1 $fastq1 \
                            --fastq2 $fastq2 \
                            --outdir $prefix/mango \
                            --bowtieref $bowtieindex \
                            --bedtoolsgenome $genomefasta \
                            --chromexclude $chrMName,chrY \
                            --stages 3:5 \
                            --shortreads FALSE \
                            --peakslop 1500 \
                            --reportallpairs TRUE > $prefix/mango/mango.stdout 2> $prefix/mango/mango.stderr
