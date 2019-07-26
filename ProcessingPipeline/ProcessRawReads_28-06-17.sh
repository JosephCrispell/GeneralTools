#!/bin/bash

# Processing Paired Raw Reads
# Author: Joseph Crispell

# Command Line Structure:
# bash ProcessRawReads.sh fileEnding trimGalore prinseq pathToRef pickRandomReads examineBlastOutput

# RUN SCRIPT USING bash script.sh (in Directory containing Raw Reads)
# Note create permissions first ---> chmod +x script.sh

####################
# DEFINE FUNCTIONS #
####################

function checkIfHelpStatementRequested {

	# Get the number of arguments from the inpur parameter
	NUMARGS=$1

	# Check if correct number of input arguments were provided
	if [ $NUMARGS -ne 5 ]
	then
		echo -e "\e[0;34m bash script designed to process paired raw read FASTQ files \e[0m"
		echo " Command line structure:"
		echo "	bash ProcessRawReads.sh fileEnding cutadapt pathToRef pickRandomReads examineBlastOutput"
		echo "		fileEnding: 		The ending of each FASTQ file that is common to all"
		echo "		cutadapt:		Would you like to use cutadapt to remove adaptors and trim the reads? [pathToCutadapt/false]"
		echo "		pathToRef:		Provide full path to indexed (bwa index) reference fasta sequence"
		echo "		pickRandomReads:	Provide full path to perl script that randomly picks unmapped reads to blast"
		echo "		examineBlastOutput:	Provide full path to perl script that examines blast output"
		echo
		echo " Requires the following tools:"
		echo "	cutadapt"
		echo "	bwa (capable of aln and mem)"
		echo "	samtools"
		echo "	BLASTn"
		echo "	bcftools"
		echo
		echo " Notes:"
		echo "	Needs access and sudo rights to the internet to run BLAST searches on de novo assembled unmapped reads of poorly mapped"
		echo "	Threshold for BLASTING unmapped reads is hard coded = 0.9. Edit script to change."
		echo "	Uses the 1st column of file names to create unique file name"
		echo "	If wanting to specify Cutadapt settings outside of this bash script. Put settings into file with start: CutadaptSettings. Use same format as used in this script"

		# Exit without an error
		exit 0
	fi
}

function checkToolsAreInstalled {

	# cutadapt - used by trim galore
	if ! type "cutadapt" > /dev/null 2>&1;
	then
		echo -e "\e[0;31m ERROR! A necessary program is not accessible: \e[0m""cutadapt"
		echo "	Follow the instructions available here: http://cutadapt.readthedocs.io/en/stable/installation.html"
		echo "	Requires python-dev:"
		echo "	sudo apt install python3-pip"
		echo "	"
		echo "	sudo pip3 install cutadapt"
		exit 0
	fi

	# bwa
	if ! type "bwa" > /dev/null 2>&1;
	then
		echo -e "\e[0;31m ERROR! A necessary program is not accessible: \e[0m""bwa"
		echo "	Download tar file from: https://sourceforge.net/projects/bio-bwa/files/"
		echo "	tar -xvf bwa-x.x.x.tar.bz2"
		echo "	cd bwa-x.x.x"
		echo "	make"
		echo "	nano ~/.bashrc"
		echo "		Add: export PATH=\$PATH:/path/to/bwa-x.x.x"
		echo "	source ~/.bashrc"
		exit 0
	fi

	# samtools
	if ! type "samtools" > /dev/null 2>&1;
	then
		echo -e "\e[0;31m ERROR! A necessary program is not accessible: \e[0m""samtools"
		echo "	Download tar file from: http://www.htslib.org/download/"
		echo "	tar -xvf samtools-x.x.tar.bz2"
		echo "	cd samtools-x.x"
		echo "	make"
		echo "	nano ~/.bashrc"
		echo "		Add: export PATH=\$PATH:/path/to/samtools-x.x"
		echo "	source ~/.bashrc"
		exit 0
	fi

	# BLASTn
	if ! type "blastn" > /dev/null 2>&1;
	then
		echo -e "\e[0;31m ERROR! A necessary program is not accessible: \e[0m""blastn"
		echo "	(Note if previously installed blastn with suao apt-get install ncbi-blast+. Remove using sudo apt-get remove ncbi-blast+)"
		echo "	Download tar file from: http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/blast/executables/LATEST/"
		echo "	Follow the instructions available here: https://www.ncbi.nlm.nih.gov/books/NBK52640/"
		echo "	tar zxvpf ncbi-blast+2.2.29-x64-linux.tar.gz"
		echo "	nano ~/.bashrc"
		echo "		Add: export PATH=\$PATH:/path/to/ncbi-blast-2.2.29+/bin"
		echo "	source ~/.bashrc"
		exit 0
	fi

	# bcftools
	if ! type "bcftools" > /dev/null 2>&1;
	then
		echo -e "\e[0;31m ERROR! A necessary program is not accessible: \e[0m""bcftools"
		echo "	Download tar file from: http://www.htslib.org/download/"
		echo "	tar -xvf bcftools-x.x.tar.bz2"
		echo "	cd bcftools-x.x"
		echo "	make"
		echo "	nano ~/.bashrc"
		echo "		Add: export PATH=\$PATH:/path/to/bcftools-x.x"
		echo "	source ~/.bashrc"
		exit 0
	fi

	# java - not used in pipeline but necessary
	if ! type "java" > /dev/null 2>&1;
	then
		echo -e "\e[0;31m ERROR! A necessary program is not accessible: \e[0m""java"
		echo "Note that java isn't necessary for this pipeline but is for next step."
		echo "	Follow the instructions available here: http://tipsonubuntu.com/2016/07/31/install-oracle-java-8-9-ubuntu-16-04-linux-mint-18/"
		echo "	sudo add-apt-repository ppa:webupd8team/java"
		echo "	sudo apt update"
		echo "	sudo apt install oracle-java8-installer"
		exit 0
	fi
}

function checkForCutadaptSettingsFile {

	# Check for Cutadapt settings file
	FOUNDCUTADAPTSETTINGS=`ls | grep "CutadaptSettings" | wc -l`
	if [ $FOUNDCUTADAPTSETTINGS == "1" ]
	then

		# Get the Cutadapt settings file
		CUTADAPTSETTINGS=`ls | grep "CutadaptSettings"`

		# Notify that found Cutadapt Settings
		echo
		echo "Found Cutadapt settings in "$CUTADAPTSETTINGS

		# Read settings from file
		while read LINE
		do

			# First adapter sequence - Illumina universal adapter
			if [ `echo $LINE | grep "UNIVERSALADAPTER1" | wc -l` == "1" ]
			then
				UNIVERSALADAPTER1=`echo $LINE | awk '{ split($0, array, "="); split(array[2], output, " "); print output[1] }'`

				# Remove quotes from string
				UNIVERSALADAPTER1="${UNIVERSALADAPTER1//\"/}"

			# First adapter sequence - Nextera Transposase sequence
			elif [ `echo $LINE | grep "UNIVERSALADAPTER2" | wc -l` == "1" ]
			then
				UNIVERSALADAPTER2=`echo $LINE | awk '{ split($0, array, "="); split(array[2], output, " "); print output[1] }'`

				# Remove quotes from string
				UNIVERSALADAPTER2="${UNIVERSALADAPTER2//\"/}"

			# Read length filter
			elif [ `echo $LINE | grep "LENGTH" | wc -l` == "1" ]
			then
				LENGTH=`echo $LINE | awk '{ split($0, array, "="); split(array[2], output, " "); print output[1] }'`

			# Read mean quality filter
			elif [ `echo $LINE | grep "QUAL" | wc -l` == "1" ]
			then
				QUAL=`echo $LINE | awk '{ split($0, array, "="); split(array[2], output, " "); print output[1] }'`

			# Number of sites to trim off left
			elif [ `echo $LINE | grep "TRIML=" | wc -l` == "1" ]
			then
				TRIML=`echo $LINE | awk '{ split($0, array, "="); split(array[2], output, " "); print output[1] }'`

			# Number of sites to trim off right
			elif [ `echo $LINE | grep "TRIMR=" | wc -l` == "1" ]
			then
				TRIMR=`echo $LINE | awk '{ split($0, array, "="); split(array[2], output, " "); print output[1] }'`

			# Threshold for proportion of Ns in reads
			elif [ `echo $LINE | grep "PROPN=" | wc -l` == "1" ]
			then
				PROPN=`echo $LINE | awk '{ split($0, array, "="); split(array[2], output, " "); print output[1] }'`
			fi

		done <$CUTADAPTSETTINGS

	fi

	echo "Using the following Cutadapt settings:"
	echo "	UNIVERSALADAPTER1 = "$UNIVERSALADAPTER1"	Illumina universal adapter sequence"
	echo "	UNIVERSALADAPTER2 = "$UNIVERSALADAPTER2"	Nextera Transposase sequence"
	echo "	READLENGTH = "$LENGTH"				The minimum length of read to be accepted"
	echo "	QUAL = "$QUAL"				Trim low-quality bases from 5' and 3' ends"
	echo "	TRIML = "$TRIML"				Trim sequence at the 5' end by x positions"
	echo "	TRIMR = "$TRIMR"				Trim sequence at the 3' end by x positions"
	echo "	PROPN = "$PROPN"				Discard reads with more than this proportion of Ns"
}

###########################
# Print Help if requested #
###########################

NUMARGS=$#
checkIfHelpStatementRequested $NUMARGS

#############################
# Get the input information #
#############################

# Get the command line information
FILEENDING=$1 # Get the common FASTQ file ending
CUTADAPT=$2 # Note whether to trim reads and remove adapters
REFERENCE=$3 # Get the path to the Reference fasta file
PICKREADS=$4 # Get path to perl script that extracts reads from SAM file
EXAMINEBLASTOUTPUT=$5 # Get path to perl script that examines blast output

# Get number of threads of computer
NTHREADS=`nproc --all`

# Create the output isolate mapping file
DATE=`date +"%d-%m-%y"`
SAMSUMMARY="isolateMappingSummary_"$DATE".txt" # Get the name of the file to input mapping summary into

# Report input settings
echo -e "\e[0;34m The following input information was provided: \e[0m"
echo -e "\e[0;32m 	Searching for FASTQ files with ending: \e[0m"$FILEENDING
echo -e "\e[0;32m 	Path to cutadapt tool: \e[0m"$CUTADAPT
echo -e "\e[0;32m 	Path to M. bovis reference sequence: \e[0m"$REFERENCE
echo -e "\e[0;32m 	Path to perl script that extracts sequences from FASTA: \e[0m"$PICKREADS
echo -e "\e[0;32m 	Path to perl script that examines BLAST output: \e[0m"$EXAMINEBLASTOUTPUT
echo
echo -e "\e[0;32m 	Produces summary File for Isolate Mapping Info: \e[0m"$SAMSUMMARY

# Add a header into the sam mapping summary file
echo "Isolate	NumberMappedReads	NumberUnmappedReads	NumberMultimappedReads" > $SAMSUMMARY # Wipe and add header to the file

##########################
# Read trimming settings #
##########################


if [ ! $CUTADAPT == "false" ]
then

	# Define Prinseq settings - USE THIS FORMAT FOR INPUT FILE IF USING ONE
	UNIVERSALADAPTER1="AGATCGGAAGAG" # Illumina universal adapter sequence (source: https://github.com/s-andrews/FastQC/blob/master/Configuration/adapter_list.txt)
	UNIVERSALADAPTER2="CTGTCTCTTATA" # Nextera Transposase sequence (source:https://github.com/s-andrews/FastQC/blob/master/Configuration/adapter_list.txt)
	LENGTH=50                        # The minimum length of read to be accepted
	QUAL=25                          # Trim low-quality bases from 5' and 3' ends
	PROPN=0.5                        # Discard reads with more than this proportion of Ns
	TRIML=20                         # Trim sequence at the 5' end by x positions
	TRIMR=5                          # Trim sequence at the 3' end by x positions

	# Check if Prinseq settings file is available with settings
	checkForCutadaptSettingsFile $CUTADAPT
fi


###########################################
# Check all necessary tools are installed #
###########################################

checkToolsAreInstalled

########################
# Find the FASTQ files #
########################

echo -e "\e[0;34m Finding raw read files... \e[0m" # Note the extra \e[0;34m ... \e[0m changes writing colour:

RAWREADS=(`ls | grep "$FILEENDING$"`) # Note the $ denotes ends with
NFILES=${#RAWREADS[@]}
NPAIRS=`expr $NFILES / 2`

#########################
# BEGIN DATA PROCESSING #
#########################

echo -e "\e[0;34m Beginning Sequencing Data Processing of $NPAIRS Samples... \e[0m"

# Create Directory to Store VCF Files - if doesn't exist
if [ ! -d "vcfFiles" ]
then
	mkdir vcfFiles
fi

# Start examine each of the FASTQ files in turn
RUN=0
for (( i=0; i<${NFILES}; i+=2))
do
	###### Prepare files #######

	# Keep track of progress and time
	RUN=`expr $RUN + 1`
	TIME=`date +"%T"`

	# Get the names of the forward and reverse reads
	FILE1=${RAWREADS[$i]}
	FILE2=${RAWREADS[$i+1]}

	#################################################################################################################
	# !!!!!!!!!!!!!!!!!!!! ADD IN A CHECK THAT COMPARES IDs AND CHECKS FOR R1 and R2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
	#################################################################################################################

	# Create Unique Prefix using first two columns
	PAIRID=`echo $FILE1 | awk '{ split($0, array, "_"); print array[1] }'`

	# Note progress and check file names
	echo -e "\e[0;34m Beginning Read Processing for Read Pair: $PAIRID ---> $RUN of $NPAIRS... \e[0m""	"$TIME
	echo "	$FILE1"
	echo "	$FILE2"

	# Unzip Files
	echo -e "\e[0;34m Unzipping Read Sequence Files... \e[0m"
	if type "pigz" > /dev/null 2>&1;
	then
		pigz --decompress $FILE1
		pigz --decompress $FILE2
	else
		gunzip $FILE1
		gunzip $FILE2
	fi

	FILE1=`echo ${FILE1:0:-3}` # Remove the .gz from the file name
	FILE2=`echo ${FILE2:0:-3}`

	####### Removing adapter sequences and trim reads #######

	if [ ! $CUTADAPT == "false" ]
	then
		echo -e "\e[0;34m Trimming reads with cutadapt... \e[0m"
		
		# Build the trimmed output files
		TRIMMED1=$PAIRID"_"$RUN"_forward_trimmed.fastq"
		TRIMMED2=$PAIRID"_"$RUN"_reverse_trimmed.fastq"

		# Trim the reads and remove adapters with cutadapt
		# -b[-B]             Adapter sequence to remove - can be ligated to 3' or 5' end on FORWARD[REVERSE] read
		# -o[-p]             Output file name for trimmed FORWARD[REVERSE] reads
		# -u[-U]             Remove bases from FORWARD[REVERSE] read. If positive - remove from beginning. If negative - remove from end.
		# --minimum-length   Discard reads shorter than this
		# --trim-n           Trim N's on ends of reads
		# --max-n            Discard reads with more than this number of N bases. If ranges between 0 and 1 - interpreted as a fraction of read length
		# --quality-cutoff   Trim low-quality bases from 5' and 3' ends of each read before adapter removal. Two values=5',3'
		# --cores            Number of CPU cores to use
		cutadapt -b $UNIVERSALADAPTER1 -B $UNIVERSALADAPTER1 -b $UNIVERSALADAPTER2 -B $UNIVERSALADAPTER2 -o $TRIMMED1 -p $TRIMMED2 $FILE1 $FILE2 --minimum-length=$LENGTH --trim-n --max-n=$PROPN --quality-cutoff=$QUAL","$QUAL -u $TRIML -u "-"$TRIMR -U $TRIML -U "-"$TRIMR --cores=$NTHREADS

		# Zip original fastq files
		if type "pigz" > /dev/null 2>&1;
		then
			pigz $FILE1
			pigz $FILE2
		else
			gzip $FILE1
			gzip $FILE2
		fi

		# Get the output files from trim galore
		echo -e "\e[0;34m Finished read trimming and adapter removal. \e[0m"
		FILE1=$TRIMMED1
		FILE2=$TRIMMED2

		# Check found the correct files
		echo "	$FILE1"
		echo "	$FILE2"
	fi

	# Get the read length - of first forward read
	READLENGTH=`sed '2q;d' $FILE1 | wc | awk '{ split($0, array, " "); print array[3] }'`

	####### Alignment #######
	
	# Generate SAM file
	# Create a name for the SAM file
	SAMFILE=$PAIRID"_"$RUN"_aln-pe.sam"

	# Run bwa mem for reads > 70bp and bwa aln for <= 70bp
	if [ $(echo "$READLENGTH > 70" | bc) -eq 1 ]
	then
		echo -e "\e[0;34m Beginning Read Alignment with bwa mem... \e[0m"

		# Generate SAM file
		# -t Flag: specify the number of cores to use
		bwa mem -t $NTHREADS $REFERENCE $FILE1 $FILE2 > $SAMFILE

	else
		echo -e "\e[0;34m Beginning Read Alignment with bwa aln... \e[0m"

		# Index the FASTQ Files
		# -I Illumina 1.3+ encoding (NOTE: This encoding covers 1.3+ to 1.8)
		# REMOVE -I flag if Illumina 1.9
		FILE1SAI=$FILE1".sai"
		FILE2SAI=$FILE2".sai"
		bwa aln -t $NTHREADS -I $REFERENCE $FILE1 > $FILE1SAI
		bwa aln -t $NTHREADS -I $REFERENCE $FILE2 > $FILE2SAI

		# Generate SAM file
		# -f Flag: SAM file to output results
		bwa sampe -f $SAMFILE $REFERENCE $FILE1SAI $FILE2SAI $FILE1 $FILE2

		# Remove unwanted files
		rm $FILE1SAI
		rm $FILE2SAI
	fi

	echo -e "\e[0;34m SAM file created. \e[0m"

	# Remove un-needed trimmed files if present
	if [ $CUTADAPT == "false" ]
	then
		if type "pigz" > /dev/null 2>&1;
		then
			pigz $FILE1
			pigz $FILE2
		else
			gzip $FILE1
			gzip $FILE2
		fi
	else
		rm $FILE1
		rm $FILE2
	fi

	####### Check Unmapped Reads #######
	# Get the information for the reads in the current SAM file
	# -F Flag: only include reads with none of the FLAGS
	# -f Flag: only include reads with all of the FLAGS
	# -c Flag: print only a count of matching records
	MAPPED=`samtools view $SAMFILE -c -F 4 --threads $NTHREADS`
	UNMAPPED=`samtools view $SAMFILE -c -f 4 --threads $NTHREADS`
	MULTIMAPPED=`samtools view $SAMFILE -c -f 256 --threads $NTHREADS`

	# HEX code table for samtools view: Base 16
	# 4096 256 16 1
	# 0    0   0  0
	# 0x0001      p       the read is paired in sequencing
	# 0x0002      P       the read is mapped in a proper pair
	# 0x0004      u       the query sequence itself is unmapped
	# 0x0008      U       the mate is unmapped
	# 0x0010      r       strand of the query (1 for reverse)
	# 0x0020      R       strand of the mate
	# 0x0040      1       the read is the first read in a pair
	# 0x0080      2       the read is the second read in a pair
	# 0x0100      s       the alignment is not primary
	# 0x0200      f       the read fails platform/vendor quality checks
	# 0x0400      d       the read is either a PCR or an optical duplicate
	# Help here: http://broadinstitute.github.io/picard/explain-flags.html

	# Print out the information for the reads in the current SAM file
	echo $PAIRID"	"$MAPPED"	"$UNMAPPED"	"$MULTIMAPPED >> $SAMSUMMARY

	# Check the proportion of mapped reads
	PROPMAPPED=`perl -E "say $MAPPED / ($MAPPED + $UNMAPPED + $MULTIMAPPED)"` # Changed to include multi-mapped reads!

	if [ $(echo " $PROPMAPPED < 0.9" | bc) -eq 1 ]
	then

		echo -e "\e[0;31m Examining unmapped reads of poorly mapped isolate: \e[0m"$PAIRID

		UNMAPPEDFILE=$PAIRID"_"$RUN"_unmapped.sam"
		RANDOMREADS=$PAIRID"_"$RUN"_randomReads.txt"
		BLASTOUTPUT=$PAIRID"_"$RUN"_BLAST.txt"
		BLASTHITS=$PAIRID"_"$RUN"_unmappedReadHits.txt"

		# BLAST some random reads from the unmapped reads file
		# -f Flag: onlu include reads with all of the FLAGS in INT present
		samtools view $SAMFILE -f 4 > $UNMAPPEDFILE # Store unmapped reads in file
		perl $PICKREADS 10 $UNMAPPEDFILE $RANDOMREADS # Pick 10 random unmapped reads
		blastn -query $RANDOMREADS -out $BLASTOUTPUT -db nr -remote # CHANGE THE NUMBER OF ALIGNMENTS!?! -num_alignments OR -max_target_seqs
		perl $EXAMINEBLASTOUTPUT 1 $BLASTOUTPUT > $BLASTHITS

		# Remove the unneccesary files
		rm $UNMAPPEDFILE
		rm $RANDOMREADS
		rm $BLASTOUTPUT
	fi

	# Convert the SAM file to a BAM file
	BAMFILE=$PAIRID"_"$RUN".bam"
	SRTDBAMFILE=$PAIRID"_"$RUN"_srtd.bam"
	BAMFILEINDX=$SRTDBAMFILE".bai"
	NDUPBAMFILE=$PAIRID"_"$RUN"_srtd_ndup.bam"

	# -b Flag: output BAM
	samtools view --threads $NTHREADS -b $SAMFILE > $BAMFILE
	samtools sort $BAMFILE -o $SRTDBAMFILE --threads $NTHREADS
	samtools index $SRTDBAMFILE -@ $NTHREADS
	samtools rmdup $SRTDBAMFILE $NDUPBAMFILE			# Should I be removing the duplicates????!!!!!
	echo -e "\e[0;34m and Converted to BAM File -> Alignment Complete. \e[0m"

	# Remove unneccessary files
	rm $SAMFILE
	rm $BAMFILE
	rm $BAMFILEINDX
	rm $SRTDBAMFILE

	####### Identify Variants #######
	echo -e "\e[0;34m Identifying Variants... \e[0m"

	# Create BCF File
	# --adjust-MQ: parameter for adjusting mapQ
	# --min-BQ: skip alignments with mapQ smaller than INT
	# --min-BQ: skip bases with baseQ/BAQ smaller than INT
	# --uncompressed: generate uncompress BCF output
	# --fasta-ref: faidx indexed reference sequence file
	BCFFILE=$PAIRID"_"$RUN".bcf"
	samtools mpileup --adjust-MQ 50 --min-MQ 30 --min-BQ 20 --uncompressed --fasta-ref $REFERENCE $NDUPBAMFILE > $BCFFILE
	echo -e "\e[0;34m BCF File Created. \e[0m"

	# Remove unnecessary files
	rm $NDUPBAMFILE

	# Convert BCF File to VCF File
	# --output-type v: output format is uncompressed vcf
	# --multiallelic-caller:  Updated and recommended calling method
	# --ploidy 1: Treat sites as haploid
	# Notes:
	#	- Calls varying and non-varying sites by default
	#	- Ignores dubious reference (N) sites
	VCFFILE=$PAIRID"_"$RUN".vcf"
	bcftools call $BCFFILE --ploidy 1 --multiallelic-caller --output-type v --threads $NTHREADS > $VCFFILE
	echo -e "\e[0;34m VCF File Created. \e[0m"
	echo -e "\e[0;34m Finished Identifying Variants. \e[0m"

	####### Moving VCF File to the VCF Directory #######
	echo -e "\e[0;34m Zipping up and moving VCF File... \e[0m"
	if type "pigz" > /dev/null 2>&1;
	then
		pigz $VCFFILE
	else
		gzip $VCFFILE
	fi
	mv $VCFFILE".gz" vcfFiles

	# Remove unnecessary files
	rm $BCFFILE

	echo -e "\e[0;34m Completed Read Processing for Read Pair: $PAIRID ---> $RUN of $NPAIRS. \e[0m"
done

echo -e "\e[0;32m Sequencing Data Processing Complete. \e[0m"
