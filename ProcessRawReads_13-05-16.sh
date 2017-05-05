#!/bin/sh
#!/bin/bash

# Processing the Paired Raw Sequencing Reads
# Author: Joseph Crisp

# Command Line Structure:
# sudo bash ProcessRawReads.sh nthreads IsolateReadSummary.txt propMapped 
# NOTE need to run this script as sudo - for the BLAST to work

# RUN SCRIPT USING bash script.sh (in Directory containing Raw Reads)
# Note create permissions first ---> chmod +x script.sh

###################
# Necessary Input #
###################

# Get the number of threads
NTHREADS=$1

# Prepare to examine the Read Mapping
SAMSUMMARY=$2
echo "Isolate	NumberMappedReads	NumberUnmappedReads	NumberMultimappedReads" > $SAMSUMMARY # Wipe and add header to the file
PICKREADS="/media/sf_UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/Tools/PickRandomReadsFromSAM.pl"
EXAMINEBLASTOUTPUT="/media/sf_UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/Tools/ExamineBLASTOutput.pl"
PROPMAPPED=$3

# Check the input
echo "Number of threads to use: "$NTHREADS
echo "Summary File for Isolate Mapping Info: "$SAMSUMMARY
echo "Min Proportion Mapped Reads to Avoid Blasting: "$PROPMAPPED

# M. Bovis Reference Fasta File - AF2122/97 (Garnier et al. 2003)
# Reference File Taken From: ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Mycobacterium_bovis_AF2122_97_uid57695/
# Create Associated Index Files: bwa index Reference.fasta
MbREF="/media/sf_UbuntuSharedFolder/Reference/NC_002945.3_AF2122-97.fasta"

# Trimming
LENGTH=50			# The minimum length of read to be accepted
MEANQUAL=20			# Filter sequence if mean quality score below x
TRIML=20			# Trim sequence at the 5' end by x positions
TRIMR=5			# Trim sequence at the 3' end by x positions
TRIMQUALL=20		# Trim sequence by quality score from the 5' end with this threshold
TRIMQUALR=20		# Trim sequence by quality score from the 3' end with this threshold score
TRIMTYPE="mean"		# Type of quality score calculation to use [min, mean, max, sum]
WINDSIZE=10			# The sliding window size used to calculate quality score by type
TRIMLTAIL=5			# Trim poly A/T > X length at 5' end
TRIMRTAIL=5			# Trim poly A/T > X length at 3' end
PRINSEQ="~/prinseq-lite-0.20.4/prinseq-lite.pl"

####### Investigate the Input Raw Read Files #######
echo -e "\e[0;34m Examining the Raw Read Files... \e[0m" # Note the extra \e[0;34m ... \e[0m changes writing colour:

# Find Files
RAWREADS=(`ls | grep "001.fastq.gz$"`) # Note the $ denotes ends with

# Create a Array containing only the FIRST PAIR of each read pair
PAIREDREADS=()
POS=-1
for VAR in "${RAWREADS[@]}"
do
	PART=`echo $VAR | awk '{ split($0, array, "_"); print array[4] }'` # Split file name by _ and select Pair number col - 3rd
	if [ $PART == "R1" ] # Only Interested in the FIRST of the read Pair
	then
		POS=`expr $POS + 1`
		PAIREDREADS[$POS]=$VAR # Append the File name in an Array
	fi
done
NOPAIREDREADS=`expr $POS + 1`

#########################
# BEGIN DATA PROCESSING #
#########################
echo -e "\e[0;34m Beginning Sequencing Data Processing of $NOPAIREDREADS Samples... \e[0m" # Note the extra \e[0;34m ... \e[0m changes writing colour:

# Create Directory to Store VCF Files
mkdir vcfFiles
RUN=0
for FILE1 in "${PAIREDREADS[@]}" # Loop through each Read Pair
do
	RUN=`expr $RUN + 1` # Keep track of which pair working on (Note expr needed to do arithmetic)

	TIME=`date +"%T"`
	
	# Identify Reciprocal Second Paired Read File
	PAIRID=`echo $FILE1 | awk '{ split($0, array, "_"); print array[1]"_"array[2]}'` # Split by _ and select 1st and 2nd cols
	FILE2=`ls | grep $PAIRID".*R2"` # Find the file name of the 2nd Pair

	# Create Unique Prefix
	PREFIX=`echo $FILE1 | awk '{ split($0, array, "_"); print array[1]"_"array[2]}'`"_$RUN" # First 2 columns
	echo -e "\e[0;34m Beginning Read Processing for Read Pair: $PREFIX ---> $RUN of $NOPAIREDREADS... \e[0m""	"$TIME
	echo "	$FILE1"
	echo "	$FILE2"

	####### File Preparation #######
	# Unzip Files
	echo -e "\e[0;34m Unzipping Read Sequence Files... \e[0m"
	gunzip $FILE1
	gunzip $FILE2
	FILE1=`echo ${FILE1:0:-3}` # Remove the .gz from the file name
	FILE2=`echo ${FILE2:0:-3}`

	####### Removing the adapter sequence #######
	# --paired: Used paired-end reads
	# --nextera: Remove Nextera adapter sequence
    # Note tim_galore is a wrapper for cutadapt - need to be installed as well
	
	echo -e "\e[0;34m Removing adapter sequences if present... \e[0m"
	~/../../media/sf_UbuntuSharedFolder/trim_galore_zip/trim_galore --paired $FILE1 $FILE2
	
	echo -e "\e[0;34m Adapter sequences removed. \e[0m"
	FILE1NOADAPT=`ls | grep $PAIRID".*_val_1"`
	FILE2NOADAPT=`ls | grep $PAIRID".*_val_2"`

	echo "	$FILE1NOADAPT"
	echo "	$FILE2NOADAPT"
	
	####### Trimming the Reads #######
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
	
	echo -e "\e[0;34m Beginning Read Trimming... \e[0m"
	perl ~/prinseq-lite-0.20.4/prinseq-lite.pl -fastq $FILE1NOADAPT -fastq2 $FILE2NOADAPT -min_len $LENGTH -min_qual_mean $MEANQUAL -trim_left $TRIML -trim_right $TRIMR -trim_qual_left $TRIMQUALL -trim_qual_right $TRIMQUALR -trim_qual_type $TRIMTYPE -trim_qual_window $WINDSIZE  

	# Find the output files
	echo -e "\e[0;34m Read Trimming Complete. \e[0m"
	TRIM1=`ls | grep $PAIRID".*R1.*prinseq_good" | grep -v "singletons"` 	# Remove .fastq and add tag from prinseq
	TRIM2=`ls | grep $PAIRID".*R2.*prinseq_good" | grep -v "singletons"`
	
	echo "	$TRIM1"
	echo "	$TRIM2"
	
	####### Alignment #######
	echo -e "\e[0;34m Beginning Read Alignment... \e[0m"
	
	# Generate SAM file
	# Note: Using the latest version of bwa algorithm so I can use the mem command
	# Have to run the bwa from where it was installed. mem is used for aligning paired end long reads > 70bp
	# Install instructions:
	# Download tar file from - http://sourceforge.net/projects/bio-bwa/files/
	# Find tar file and type into command line: 
	#		tar -xvf bwa-0.7.10.tar.bz2
	# 		cd bwa-0.7.10
	# 		make
	# -t Flag: number of threads (taken in from the command line)
	
	SAMFILE=$PREFIX"_aln-pe.sam"
	~/bwa-0.7.10/bwa mem -t $1 $MbREF $TRIM1 $TRIM2 > $SAMFILE
	echo -e "\e[0;34m SAM file created. \e[0m"

	# Get the information for the reads in the current SAM file
	MAPPED=`samtools view -S $SAMFILE -c -F 4`
	UNMAPPED=`samtools view -S $SAMFILE -c -f 4`
	MULTIMAPPED=`samtools view -S $SAMFILE -c -f 256`
	
	# Print out the information for the reads in the current SAM file
	echo $PREFIX"	"$MAPPED"	"$UNMAPPED"	"$MULTIMAPPED >> $SAMSUMMARY
	
	# Check the proportion of mapped reads
	VALUE=`perl -E "say $MAPPED / ($MAPPED + $UNMAPPED)"`

	if [ $(echo " $VALUE < $PROPMAPPED" | bc) -eq 1 ] 
	then
		UNMAPPEDFILE=$PREFIX"_unmapped.sam"
		RANDOMREADS=$PREFIX"_randomReads.txt"
		BLASTOUTPUT=$PREFIX"_BLAST.txt"
		BLASTHITS=$PREFIX"_unmappedReadHits.txt"
		
		# BLAST some random reads from the unmapped reads file
		samtools view -S $SAMFILE -f 4 > $UNMAPPEDFILE # Store unmapped reads in file
		perl $PICKREADS 10 $UNMAPPEDFILE $RANDOMREADS # Pick 10 random unmapped reads
		blastn -query $RANDOMREADS -out $BLASTOUTPUT -db nr -remote
		perl $EXAMINEBLASTOUTPUT 1 $BLASTOUTPUT > $BLASTHITS
		
		# Remove the unneccesary files
		rm $UNMAPPEDFILE
		rm $RANDOMREADS
		rm $BLASTOUTPUT
	fi
	
	# Convert the SAM file to a BAM file
	BAMFILE=$PREFIX".bam"
	SRTDBAMPREFIX=$PREFIX"_srtd"
	SRTDBAMFILE=$SRTDBAMPREFIX".bam"
	BAMFILEINDX=$SRTDBAMFILE".bai"
	NDUPBAMFILE=$SRTDBAMPREFIX"_ndup.bam"

	# -S Flag: input is SAM
	# -u Flag: uncompressed BAM output (force -b)
	# -o Flag: output file name [stdout]
	samtools view -S -u -o $BAMFILE $SAMFILE
	samtools sort $BAMFILE $SRTDBAMPREFIX
	samtools index $SRTDBAMFILE
	samtools rmdup $SRTDBAMFILE $NDUPBAMFILE			# Should I be removing the duplicates????!!!!!
	echo -e "\e[0;34m and Converted to BAM File -> Alignment Complete. \e[0m"

	####### Identify Variants #######h
	echo -e "\e[0;34m Identifying Variants... \e[0m"
	# Create BCF File
	# -P Flag: comma separated list of platforms for indels [all]
	# -C Flag: parameter for adjusting mapQ
	# -q Flag: skip alignments with mapQ smaller than INT
	# -Q Flag: skip bases with baseQ/BAQ smaller than INT
	# -u Flag: generate uncompress BCF output
	# -f Flag: faidx indexed reference sequence file
	BCFFILE=$PREFIX".bcf"
	samtools mpileup -P ILLUMINA -C50 -q30 -Q20 -u -f $MbREF $NDUPBAMFILE > $BCFFILE
	echo -e "\e[0;34m BCF File Created. \e[0m"

	# Convert BCF File to VCF File
	# -c Flag: SNP calling
	# -g Flag: call genotypes at variant sites
	# -N Flag: skip sites where REF is not A/C/G/T
	# -I Flag: skip indels
	
	VCFFILE=$PREFIX".vcf"
	bcftools view -c -g -N $BCFFILE > $VCFFILE
	echo -e "\e[0;34m VCF File Created. \e[0m"
	echo -e "\e[0;34m Finished Identifying Variants. \e[0m"

	####### Copy VCF File to the VCF Directory #######
	echo -e "\e[0;34m Moving VCF File... \e[0m"
	mv $VCFFILE vcfFiles

	####### Zip up Raw Read Files #######
	gzip $FILE1
	gzip $FILE2
	
	####### Remove all Unnecessary Files Created #######
	echo -e "\e[0;34m Removing Unnecessary Output Files... \e[0m"
	rm $SAMFILE
	rm $BAMFILE
	rm $BAMFILEINDX
	rm $SRTDBAMFILE
	rm $NDUPBAMFILE
	rm $BCFFILE
	rm *fastq # Removes all files produced from the trimming Process
	rm $FILE1NOADAPT
	rm $FILE2NOADAPT
	rm *trimming_report*
	rm *trimmed*
	
	echo -e "\e[0;34m Completed Read Processing for Read Pair: $PREFIX ---> $RUN of $NOPAIREDREADS. \e[0m"
	echo ""
done

echo -e "\e[0;32m Sequencing Data Processing Complete. \e[0m"
