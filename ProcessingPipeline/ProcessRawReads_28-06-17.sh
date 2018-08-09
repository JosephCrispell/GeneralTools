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
	if [ $1 == "-h" -o $1 == "-help"  -o $1 == "help" -o $1 == "" ]
	then
		echo -e "\e[0;34m bash script designed to process paired raw read FASTQ files \e[0m"
		echo " Command line structure:"
		echo "	bash ProcessRawReads.sh fileEnding trimGalore prinseq pathToRef pickRandomReads examineBlastOutput"
		echo "	fileEnding: 		The ending of each FASTQ file that is common to all"
		echo "	trimGalore:		Would you like to use trim_galore to remove adaptors? [pathToTrimGalore/false]"
		echo "	prinseq:		Would you like to use prinseq to trim the reads? [pathToPrinseq/false]"
		echo "	pathToRef:		Provide full path to indexed (bwa index) reference fasta sequence"
		echo "	pickRandomReads:	Provide full path to perl script that randomly picks unmapped reads to blast"
		echo "	examineBlastOutput:	Provide full path to perl script that examines blast output"
		echo 
		echo " Requires the following tools:"
		echo "	cutadapt"
		echo "	bwa (capable of aln and mem)"
		echo "	samtools"
		echo "	BLASTn"
		echo "	bcftools"
		echo 
		echo " Notes:"
		echo "	Needs access to the internet to run BLAST searches on unmapped reads of poorly mapped"
		echo "	Threshold for BLASTING unmapped reads is hard coded. Edit script to change."
		echo "	Uses the 1st and 2nd columns of file names to create unique file name"
		echo "	prinseq settings are hard-coded. Edit script to change."
		echo "	Only examines unmapped reads if proportion mapped drops below 0.9"
		echo "	Requires sudo rights and internet to examine unmapped reads"
		echo "	If wanting to specify Prinseq settings outside of this bash script. Put settings into file with start: PrinseqSettings. Use same format as used in this script"
	
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
		echo "	sudo apt-get install python-dev"
		echo "	sudo apt install python-pip"
		echo "	pip install --user --upgrade cutadapt"
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

function checkForPrinseqSettingsFile {
	
	# Check for Prinseq settings file
	FOUNDPRINSEQSETTINGS=`ls | grep "PrinseqSettings" | wc -l`
	if [ $FOUNDPRINSEQSETTINGS == "1" ]
	then
   
		# Get the Prinseq settings file
		PRINSEQSETTINGS=`ls | grep "PrinseqSettings"`
	
		# Notify that found Prinseq Settings
		echo
		echo "Found Prinseq settings in "$PRINSEQSETTINGS
   
		# Read settings from file
		while read LINE
		do
		
			# Read length filter
			if [ `echo $LINE | grep "LENGTH" | wc -l` == "1" ]
			then
				LENGTH=`echo $LINE | awk '{ split($0, array, "="); split(array[2], output, " "); print output[1] }'`
		
			# Read mean quality filter
			elif [ `echo $LINE | grep "MEANQUAL" | wc -l` == "1" ]
			then
				MEANQUAL=`echo $LINE | awk '{ split($0, array, "="); split(array[2], output, " "); print output[1] }'`
			
			# Number of sites to trim off left
			elif [ `echo $LINE | grep "TRIML=" | wc -l` == "1" ]
			then
				TRIML=`echo $LINE | awk '{ split($0, array, "="); split(array[2], output, " "); print output[1] }'`
			
			# Number of sites to trim off right
			elif [ `echo $LINE | grep "TRIMR=" | wc -l` == "1" ]
			then
				TRIMR=`echo $LINE | awk '{ split($0, array, "="); split(array[2], output, " "); print output[1] }'`
		
			# Quality filter for sliding window working from left
			elif [ `echo $LINE | grep "TRIMQUALL" | wc -l` == "1" ]
			then
				TRIMQUALL=`echo $LINE | awk '{ split($0, array, "="); split(array[2], output, " "); print output[1] }'`
				
			# Quality filter for sliding window working from right
			elif [ `echo $LINE | grep "TRIMQUALR" | wc -l` == "1" ]
			then
				TRIMQUALR=`echo $LINE | awk '{ split($0, array, "="); split(array[2], output, " "); print output[1] }'`
				
			# Summary statistic for sliding window
			elif [ `echo $LINE | grep "TRIMTYPE" | wc -l` == "1" ]
			then
				TRIMTYPE=`echo $LINE | awk '{ split($0, array, "="); split(array[2], output, " "); print output[1] }'`
				TRIMTYPE="${TRIMTYPE//\"/}"
				
			# Sliding window size
			elif [ `echo $LINE | grep "WINDSIZE" | wc -l` == "1" ]
			then
				WINDSIZE=`echo $LINE | awk '{ split($0, array, "="); split(array[2], output, " "); print output[1] }'`
				
			# Size of nucleotide repeat to filter from left
			elif [ `echo $LINE | grep "TRIMLTAIL" | wc -l` == "1" ]
			then
				TRIMLTAIL=`echo $LINE | awk '{ split($0, array, "="); split(array[2], output, " "); print output[1] }'`
			
			# Size of nucleotide repeat to filter from right
			elif [ `echo $LINE | grep "TRIMRTAIL" | wc -l` == "1" ]
			then
				TRIMRTAIL=`echo $LINE | awk '{ split($0, array, "="); split(array[2], output, " "); print output[1] }'`
			fi
		
		done <$PRINSEQSETTINGS
		
		echo "Using the following Prinseq settings:"
		echo "	READLENGTH = "$LENGTH"			The minimum length of read to be accepted"
		echo "	MEANQUAL = "$MEANQUAL"	Filter sequence if mean quality score below x"
		echo "	TRIML = "$TRIML"	Trim sequence at the 5' end by x positions"
		echo "	TRIMR = "$TRIMR"	Trim sequence at the 3' end by x positions"
		echo "	TRIMQUALL = "$TRIMQUALL"	Trim sequence by quality score from the 5' end with this threshold score"
		echo "	TRIMQUALR = "$TRIMQUALR"	Trim sequence by quality score from the 3' end with this threshold score"
		echo "	TRIMTYPE = "$TRIMTYPE"		Type of quality score calculation to use [min, mean, max, sum]"
		echo "	WINDSIZE = "$WINDSIZE"		The sliding window size used to calculate quality score by type"
		echo "	TRIMLTAIL = "$TRIMLTAIL"	Trim poly A/T > X length at 5' end"
		echo "	TRIMRTAIL = "$TRIMRTAIL"	Trim poly A/T > X length at 3' end"
   
	else
	
		echo
		echo "No Prinseq settings file found, using the following default settings:"
		echo "	READLENGTH = "$LENGTH"	The minimum length of read to be accepted"
		echo "	MEANQUAL = "$MEANQUAL"		Filter sequence if mean quality score below x"
		echo "	TRIML = "$TRIML"		Trim sequence at the 5' end by x positions"
		echo "	TRIMR = "$TRIMR"		Trim sequence at the 3' end by x positions"
		echo "	TRIMQUALL = "$TRIMQUALL"		Trim sequence by quality score from the 5' end with this threshold score"
		echo "	TRIMQUALR = "$TRIMQUALR"		Trim sequence by quality score from the 3' end with this threshold score"
		echo "	TRIMTYPE = "$TRIMTYPE"	Type of quality score calculation to use [min, mean, max, sum]"
		echo "	WINDSIZE = "$WINDSIZE"		The sliding window size used to calculate quality score by type"
		echo "	TRIMLTAIL = "$TRIMLTAIL"		Trim poly A/T > X length at 5' end"
		echo "	TRIMRTAIL = "$TRIMRTAIL"		Trim poly A/T > X length at 3' end"
	fi
}

###########################
# Print Help if requested #
###########################

checkIfHelpStatementRequested $1

#############################
# Get the input information #
#############################

FILEENDING=$1 # Get the common FASTQ file ending
TRIMGALORE=$2 # Note whether to remove adaptors with trim_galore
PRINSEQ=$3 # Get the path to the prinseq tool
REFERENCE=$4 # Get the path to the Reference fasta file
PICKREADS=$5 # Get path to perl script that randomly picks unmapped reads to blast
EXAMINEBLASTOUTPUT=$6 # Get path to perl script that examines blast output

# Create the output isolate mapping file
DATE=`date +"%d-%m-%y"`
SAMSUMMARY="isolateMappingSummary_"$DATE".txt" # Get the name of the file to input mapping summary into

# Report input settings
echo -e "\e[0;34m The following input information was provided: \e[0m"
echo -e "\e[0;32m 	Searching for FASTQ files with ending: \e[0m"$FILEENDING
echo -e "\e[0;32m 	Path to trim_galore tool: \e[0m"$TRIMGALORE
echo -e "\e[0;32m 	Path to prinseq tool: \e[0m"$PRINSEQ
echo -e "\e[0;32m 	Path to M. bovis reference sequence: \e[0m"$REFERENCE
echo -e "\e[0;32m 	Path to perl script that randomly picks reads for BLAST: \e[0m"$PICKREADS
echo -e "\e[0;32m 	Path to perl script that examines BLAST output: \e[0m"$EXAMINEBLASTOUTPUT
echo
echo -e "\e[0;32m 	Produces summary File for Isolate Mapping Info: \e[0m"$SAMSUMMARY

# Add a header into the sam mapping summary file
echo "Isolate	NumberMappedReads	NumberUnmappedReads	NumberMultimappedReads" > $SAMSUMMARY # Wipe and add header to the file

##########################
# Read trimming settings #
##########################


if [ ! $PRINSEQ == "false" ]
then
	
	# Define Prinseq settings - USE THIS FORMAT FOR INPUT FILE IF USING ONE
	LENGTH=50			# The minimum length of read to be accepted
	MEANQUAL=20			# Filter sequence if mean quality score below x
	TRIML=20			# Trim sequence at the 5' end by x positions
	TRIMR=5				# Trim sequence at the 3' end by x positions
	TRIMQUALL=20		# Trim sequence by quality score from the 5' end with this threshold
	TRIMQUALR=20		# Trim sequence by quality score from the 3' end with this threshold score
	TRIMTYPE="mean"		# Type of quality score calculation to use [min, mean, max, sum]
	WINDSIZE=10			# The sliding window size used to calculate quality score by type
	TRIMLTAIL=5			# Trim poly A/T > X length at 5' end
	TRIMRTAIL=5			# Trim poly A/T > X length at 3' end

	# Check if Prinseq settings file is available with settings
	checkForPrinseqSettingsFile $PRINSEQ
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
	
	# Create Unique Prefix using first two columns
	PAIRID=`echo $FILE1 | awk '{ split($0, array, "_"); print array[1]"_"array[2] }'`
	
	# Note progress and check file names
	echo -e "\e[0;34m Beginning Read Processing for Read Pair: $PAIRID ---> $RUN of $NPAIRS... \e[0m""	"$TIME
	echo "	$FILE1"
	echo "	$FILE2"
		
	# Unzip Files
	echo -e "\e[0;34m Unzipping Read Sequence Files... \e[0m"
	gunzip $FILE1
	gunzip $FILE2
	FILE1=`echo ${FILE1:0:-3}` # Remove the .gz from the file name
	FILE2=`echo ${FILE2:0:-3}`
	
	####### Removing the adapter sequence #######
	# --paired: Used paired-end reads
	if [ ! $TRIMGALORE == "false" ]
	then
		echo -e "\e[0;34m Removing adapter sequences if present... \e[0m"
		$TRIMGALORE --paired $FILE1 $FILE2
	
		# Zip original fastq files
		gzip $FILE1
		gzip $FILE2
	
		# Get the output files from trim galore
		echo -e "\e[0;34m Adapter sequences removed. \e[0m"
		FILE1=`ls | grep $PAIRID".*_val_1"`
		FILE2=`ls | grep $PAIRID".*_val_2"`

		# Remove the unecessary output files
		rm *trim*
		
		# Check found the correct files
		echo "	$FILE1"
		echo "	$FILE2"
	fi
	
	####### Trim the Reads #######
	# Using program Prinseq -> Trimming reads to remove low quality bases according to the Phred Scoring system
	# min_len: The minimum length of read to be accepted
	# min_qual_mean: Filter sequence if mean quality score below x
	# trim_left: Trim sequence at the 5' end by x positions
	# trim_right: Trim sequence at the 3' end by x positions	# trim_qual_left: Trim sequence by quality score from the 5' end with this threshold score
	# trim_qual_left: Trim sequence by quality score from the 3' end with this threshold score
	# trim_qual_right: Trim sequence by quality score from the 3' end with this threshold score
	# trim_qual_type: Type of quality score calculation to use [min, mean, max, sum]
	# trim_qual_window: The sliding window size used to calculate quality score by type. To stop at the first base that fails the rule defined,
	#					use a window size of 1
	if [ ! $PRINSEQ == "false" ]
	then
		echo -e "\e[0;34m Beginning Read Trimming... \e[0m"
		perl $PRINSEQ -fastq $FILE1 -fastq2 $FILE2 -min_len $LENGTH -min_qual_mean $MEANQUAL -trim_left $TRIML -trim_right $TRIMR -trim_qual_left $TRIMQUALL -trim_qual_right $TRIMQUALR -trim_qual_type $TRIMTYPE -trim_qual_window $WINDSIZE  

		# Remove files previously stored under FILE1 and FILE2
		if [ ! $TRIMGALORE == "false" ]
		then
			rm $FILE1
			rm $FILE2
		else
			gzip $FILE1
			gzip $FILE2
		fi
		
		# Find the output files
		echo -e "\e[0;34m Read Trimming Complete. \e[0m"
		TRIMFILES=(`ls | grep $PAIRID".*prinseq_good" | grep -v "singletons"`)
		FILE1=${TRIMFILES[0]}
		FILE2=${TRIMFILES[1]}
	
		echo "	$FILE1"
		echo "	$FILE2"
		
		# Remove unecessary files
		rm *prinseq_bad*
		rm *prinseq_good_singletons*
	fi
	
	# Get the read length - of first forward read
	READLENGTH=`sed '2q;d' $FILE1 | wc | awk '{ split($0, array, " "); print array[3] }'`
	
	####### Alignment #######
	# Get number of cores of computer
	NCORES=`nproc --all`
	
	# Generate SAM file
	# Create a name for the SAM file
	SAMFILE=$PAIRID"_"$RUN"_aln-pe.sam"
	
	# Run bwa mem for reads > 70bp and bwa aln for <= 70bp
	if [ $(echo "$READLENGTH > 70" | bc) -eq 1 ] 
	then
		echo -e "\e[0;34m Beginning Read Alignment with bwa mem... \e[0m"
		
		# Generate SAM file
		# -t Flag: specify the number of cores to use
		bwa mem -t $NCORES $REFERENCE $FILE1 $FILE2 > $SAMFILE
	
	else
		echo -e "\e[0;34m Beginning Read Alignment with bwa aln... \e[0m"
	
		# Index the FASTQ Files
		# -I Illumina 1.3+ encoding (NOTE: This encoding covers 1.3+ to 1.8)
		# REMOVE -I flag if Illumina 1.9
		FILE1SAI=$FILE1".sai"
		FILE2SAI=$FILE2".sai"
		bwa aln -I $REFERENCE $FILE1 > $FILE1SAI
		bwa aln -I $REFERENCE $FILE2 > $FILE2SAI
	
		# Generate SAM file
		# -f Flag: SAM file to output results
		bwa sampe -f $SAMFILE $REFERENCE $FILE1SAI $FILE2SAI $FILE1 $FILE2
		
		# Remove unwanted files
		rm $FILE1SAI
		rm $FILE2SAI
	fi
	
	echo -e "\e[0;34m SAM file created. \e[0m"
	
	# Remove un-needed files
	if [ $TRIMGALORE == "false" ] && [ $PRINSEQ == "false" ]
	then
		gzip $FILE1
		gzip $FILE2
	else
		rm $FILE1
		rm $FILE2
	fi
	
	####### Check Unmapped Reads #######
	# Get the information for the reads in the current SAM file
	# -F Flag: only include reads with none of the FLAGS
	# -f Flag: only include reads with all of the FLAGS
	# -c Flag: print only a count of matching records
	MAPPED=`samtools view $SAMFILE -c -F 4`
	UNMAPPED=`samtools view $SAMFILE -c -f 4`
	MULTIMAPPED=`samtools view $SAMFILE -c -f 256`
	
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
		blastn -query $RANDOMREADS -out $BLASTOUTPUT -db nr -remote
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
	samtools view -b $SAMFILE > $BAMFILE
	samtools sort $BAMFILE -o $SRTDBAMFILE
	samtools index $SRTDBAMFILE
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
	bcftools call $BCFFILE --ploidy 1 --multiallelic-caller --output-type v > $VCFFILE 
	echo -e "\e[0;34m VCF File Created. \e[0m"
	echo -e "\e[0;34m Finished Identifying Variants. \e[0m"
	
	####### Moving VCF File to the VCF Directory #######
	echo -e "\e[0;34m Zipping up and moving VCF File... \e[0m"
	gzip $VCFFILE
	mv $VCFFILE".gz" vcfFiles
	
	# Remove unnecessary files
	rm $BCFFILE
	
	echo -e "\e[0;34m Completed Read Processing for Read Pair: $PAIRID ---> $RUN of $NPAIRS. \e[0m"
done

echo -e "\e[0;32m Sequencing Data Processing Complete. \e[0m"
