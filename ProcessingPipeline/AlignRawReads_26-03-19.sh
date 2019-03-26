#!/bin/bash

###############
# SCRIPT INFO #
###############

# Align paired end reads against reference genome
# Author: Joseph Crispell

# Command Line Structure:
# bash AlignRawReads_DATE.sh forward.fastq.gz reverse.fastq.gz reference

# Arguments:
# forward.fastq.gz	Compressed FASTQ file containing the forward reads (file name should contain _R1_)
# reverse.fastq.gz	Compressed FASTQ file containing the reverse reads (file name should contain _R2_)
# reference         Reference genome file - must be indexed (see end of TOOL INSTALLATION section below)

#####################
# TOOL INSTALLATION #
#####################

## Installing cutadapt
# $ sudo apt install python3-pip
# $ sudo pip3 install cutadapt

## Installing bwa
# Download tar file from: https://sourceforge.net/projects/bio-bwa/files/
# $ tar -xvf bwa-x.x.x.tar.bz2
# $ cd bwa-x.x.x
# $ make
# $ export PATH=\$PATH:/path/to/bwa-x.x.x # Add this to ~/.bashrc file

## Installing samtools
# Download tar file from: http://www.htslib.org/download/
# $ tar -xvf samtools-x.x.tar.bz2
# $ cd samtools-x.x
# $ make
# $ export PATH=\$PATH:/path/to/samtools-x.x # Add this to ~/.bashrc file

## Installing bcftools
# Download tar file from: http://www.htslib.org/download/
# $ tar -xvf bcftools-x.x.tar.bz2
# $ cd bcftools-x.x
# $ make
# $ export PATH=\$PATH:/path/to/bcftools-x.x # Add this to ~/.bashrc file

## Adding to ~/.bashrc file
# $ nano ~/.bashrc
# Paste in line interested in and save and exit file
# $ source ~/.bashrc

## Reference sequence
# Download the FASTA file for your desired reference genome
# $ bwa index reference.fasta

#############
# ARGUMENTS #
#############

# Get the forward and reverse reads
FORWARD=$1
REVERSE=$2

# Get the reference sequence
REFERENCE=$3

# Set the number of threads to use throughout - how many tasks can be done simultaneously
NTHREADS=3

#########
# UNZIP #
#########

# Print progress information
echo -e "\e[0;34m Unzipping raw reads... \e[0m"

# Unzip the forward and reverse reads
gunzip $FORWARD & gunzip $REVERSE

# Note names of unzipped files
FORWARD=`echo ${FORWARD:0:-3}` # Remove the .gz from the file name
REVERSE=`echo ${REVERSE:0:-3}`

############
# TRIMMING #
############

# Print progress information
echo -e "\e[0;34m Starting to trim raw reads... \e[0m"

# Set your trimming parameters here
QUALITY=25 # Quality threshold used to filter poor quality reads
MAXPROPN=0.25 # Threshold proportion of Ns in read (sites with no information) used to filter poor quality reads
MINLENGTH=50 # Threshold to filter reads out that are too short
TRIMLEFTFORWARD=10 # Number of sites to remove from left of forward reads
TRIMRIGHTFORWARD=10 # Number of sites to remove from right of forward reads
TRIMLEFTREVERSE=10 # Number of sites to remove from left of forward reads
TRIMRIGHTREVERSE=10 # Number of sites to remove from right of forward reads

# Name the output files
TRIMMEDFORWARD="trimmed_forward.fastq"
TRIMMEDREVERSE="trimmed_reverse.fastq"

# Run cutadapt
cutadapt -b AGATCGGAAGAG -B AGATCGGAAGAG -b CTGTCTCTTATA -B CTGTCTCTTATA -o $TRIMMEDFORWARD -p $TRIMMEDREVERSE $FORWARD $REVERSE --minimum-length=$MINLENGTH --trim-n --max-n=$MAXPROPN --quality-cutoff=$QUALITY,$QUALITY -u $TRIMLEFTFORWARD -u -$TRIMRIGHTFORWARD -U $TRIMLEFTREVERSE -U -$TRIMRIGHTREVERSE --cores=$NTHREADS

############
# ALIGNING #
############

# Print progress information
echo -e "\e[0;34m Starting to align trimmed reads... \e[0m"

# Name the output file
SAMFILE="aligned.sam"

# Run the alignment
bwa mem -t $NTHREADS $REFERENCE $TRIMMEDFORWARD $TRIMMEDREVERSE > $SAMFILE

# Convert SAM to BAM
BAMFILE="aligned.bam"
samtools view --threads $NTHREADS -b $SAMFILE > $BAMFILE

# Sort the aligned reads
SORTEDBAMFILE="aligned_sorted.bam"
samtools sort $BAMFILE -o $SORTEDBAMFILE --threads $NTHREADS

# Index the sorted reads
samtools index $SRTDBAMFILE -@ $NTHREADS

# Remove duplicates
NODUPLICATES="aligned_sorted_rmDuplicates.bam"
samtools rmdup $SRTDBAMFILE $NODUPLICATES

###################
# VARIANT CALLING #
###################

# Print progress information
echo -e "\e[0;34m Starting variant calling... \e[0m"

# Set some quality filters
ADJUST=50 # Parameter to adjust quality scores
MAPQAULITY=30 # Skip aligned reads with quality score less than this
BASEQUALITY=20 # Skip bases with quality score less than this

# Create a BCF file (compressed form of VCF)
BCFFILE="variants.bcf"
samtools mpileup --adjust-MQ $ADJUST --min-MQ $MAPQUALITY --min-BQ $BASEQUALITY --uncompressed --fasta-ref $REFERENCE $NODUPLICATES > $BCFFILE

# Convert BCF to VCF
VCFFILE="variants.vcf"
bcftools call $BCFFILE --ploidy 1 --multiallelic-caller --output-type v --threads $NTHREADS > $VCFFILE

echo -e "\e[0;32m Finished! :-) \e[0m"