#!/bin/bash

###############
# SCRIPT INFO #
###############

# Align paired end reads against reference genome
# Author: Joseph Crispell

# Command Line Structure:
# bash AlignRawReads_DATE.sh forward.fastq.gz reverse.fastq.gz reference.fasta annotation.gff3

# Arguments:
# forward.fastq.gz	Compressed FASTQ file containing the forward reads (file name should contain _R1_)
# reverse.fastq.gz	Compressed FASTQ file containing the reverse reads (file name should contain _R2_)
# reference         Reference genome file - must be indexed (see end of TOOL INSTALLATION section below)
# annotations.gff3	Full path to annotation file for the reference genome (NOTE M. bovis NCBI gff3 file doesn't work properly)

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

##################
# HELP STATEMENT #
##################

# Print help statement when more/less than three input arguments provided
if test "$#" -ne 4; then

    echo -e "\e[0;34m bash script designed to process a pair of FASTQ files \e[0m"
	echo " Command line structure:"
	echo "	bash AlignRawReads_DATE.sh forward.fastq.gz reverse.fastq.gz reference.fasta annotations.gff3"
	echo "		forward.fastq.gz      Full path to FASTQ files containing the forward reads (usually has \"_R1_\" in name)"
	echo "		reverse.fastq.gz      Full path to FASTQ files containing the reverse reads (usually has \"_R2_\" in name)"
	echo "		reference.fasta       Full path to reference genome FASTA file. Note this should be indexed with bwa index."
	echo "		annotations.gff3      Full path to annotation file for the reference genome"
	exit 0
fi

#############
# ARGUMENTS #
#############

# Get the forward and reverse reads
FORWARD=$1
REVERSE=$2

# Get the reference sequence
REFERENCE=$3

# Get the annotation file for the reference
ANNOTATION=$4

# Set the number of threads to use throughout - how many tasks can be done simultaneously
NTHREADS=10

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
TRIMRIGHTFORWARD=0 # Number of sites to remove from right of forward reads
TRIMLEFTREVERSE=10 # Number of sites to remove from left of forward reads
TRIMRIGHTREVERSE=0 # Number of sites to remove from right of forward reads

# Name the output files
TRIMMEDFORWARD="trimmed_forward.fastq"
TRIMMEDREVERSE="trimmed_reverse.fastq"

# Run cutadapt
cutadapt -b AGATCGGAAGAG -B AGATCGGAAGAG -b CTGTCTCTTATA -B CTGTCTCTTATA -o $TRIMMEDFORWARD -p $TRIMMEDREVERSE $FORWARD $REVERSE --minimum-length=$MINLENGTH --trim-n --max-n=$MAXPROPN --quality-cutoff=$QUALITY,$QUALITY -u $TRIMLEFTFORWARD -u -$TRIMRIGHTFORWARD -U $TRIMLEFTREVERSE -U -$TRIMRIGHTREVERSE --cores=$NTHREADS

echo -e "\e[0;34m Finished trimming. Zipping up original FASTQ files... \e[0m"

# Zip up the original untrimmed FASTQ files
gzip $FORWARD & gzip $REVERSE

############
# ALIGNING #
############

# Print progress information
echo -e "\e[0;34m Finished zipping FASTQ files. Starting to align trimmed reads - creating SAM file... \e[0m"

# Name the output file
SAMFILE="aligned.sam"

# Run the alignment
bwa mem -t $NTHREADS $REFERENCE $TRIMMEDFORWARD $TRIMMEDREVERSE > $SAMFILE

echo -e "\e[0;34m Finished aligning and created SAM file. Converting SAM file to BAM file... \e[0m"

# Convert SAM to BAM
BAMFILE="aligned.bam"
samtools view --threads $NTHREADS -b $SAMFILE > $BAMFILE

echo -e "\e[0;34m Finished creating BAM file. Sorting BAM file... \e[0m"

# Sort the aligned reads
SORTEDBAMFILE="aligned_sorted.bam"
samtools sort $BAMFILE -o $SORTEDBAMFILE --threads $NTHREADS

echo -e "\e[0;34m Finished sorting BAM file. Indexing BAM file... \e[0m"

# Index the sorted reads
samtools index $SORTEDBAMFILE -@ $NTHREADS

echo -e "\e[0;34m Finished indexing BAM file. Removing duplicates from BAM file... \e[0m"

# Remove duplicates
NODUPLICATES="aligned_sorted_rmDuplicates.bam"
samtools rmdup $SORTEDBAMFILE $NODUPLICATES

###################
# VARIANT CALLING #
###################

# Print progress information
echo -e "\e[0;34m Finished removing duplicates. Starting variant calling - creating BCF file... \e[0m"

# Set some quality filters
ADJUST=50 # Parameter to adjust quality scores
MAPQUALITY=30 # Skip aligned reads with quality score less than this
BASEQUALITY=20 # Skip bases with quality score less than this

# Create a BCF file (compressed form of VCF)
BCFFILE="variants.bcf"
bcftools mpileup $NODUPLICATES --threads $NTHREADS --adjust-MQ $ADJUST --min-MQ $MAPQUALITY --min-BQ $BASEQUALITY --fasta-ref $REFERENCE --output-type b --output $BCFFILE

echo -e "\e[0;34m Created BCF file. Converting BCF file to VCF file... \e[0m"

# Convert BCF to VCF
VCFFILE="variants.vcf"
bcftools call $BCFFILE --ploidy 1 --multiallelic-caller --output-type v --threads $NTHREADS --variants-only --output-type v --output $VCFFILE

echo -e "\e[0;34m Finished converting BCF file to VCF file. Running consequence calling on variants... \e[0m"

# Consequence calling
# Note that the organism name in the the FASTA file (the #CHROM name in the VCF file) must match the first column of the gff file
CONSEQUENCEFILE="consequence.tsv"
bcftools csq --fasta-ref $REFERENCE --gff-annot $ANNOTATION --output-type t --output $CONSEQUENCEFILE $VCFFILE

# Remove the intermediate files
#echo -e "\e[0;34m Finished consequence calling. Removing unecessary intermediate files... \e[0m"
#rm aligned.*
#rm trimmed_*
#rm aligned_*.bam*
#rm variants.bcf

echo -e "\e[0;32m Finished! :-) \e[0m"
