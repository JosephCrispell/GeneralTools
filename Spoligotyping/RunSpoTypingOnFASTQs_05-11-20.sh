#!/bin/bash

# Run SpoTyping on raw paired end FASTQ files
# Author: Joseph Crispell

# Spotyping can be downloaded from here: https://github.com/xiaeryu/SpoTyping/releases/

# Command Line Structure (run in a directory of FASTQ files):
# bash RunSpoTypingOnFASTQs_DATE.sh pathToSpotyper

# RUN SCRIPT USING bash script.sh (in Directory containing Raw Reads)
# Note create permissions first ---> chmod +x script.sh

###########################################
# Get path to spotyping from command line #
###########################################

SPOTYPING=$1

########################
# Find the FASTQ files #
########################

echo -e "\e[0;34m Finding raw read files... \e[0m" # Note the extra \e[0;34m ... \e[0m changes writing colour:
FILEENDING=".fastq.gz"
RAWREADS=(`ls | grep "$FILEENDING$"`) # Note the $ denotes ends with
NFILES=${#RAWREADS[@]}
NPAIRS=`expr $NFILES / 2`

#########################
# Create an output file #
#########################

RESULTS="SpotypingResults.txt"
echo "Files	Code	Something" > $RESULTS

##############################################
# EXAMINE EACH SET OF PAIRED END FASTQ FILES #
##############################################

echo -e "\e[0;34m Beginning to examine $NPAIRS paired end FASTQ files... \e[0m"

# Start examine each of the FASTQ files in turn
RUN=0
for (( i=0; i<${NFILES}; i+=2))
do
	###### Getting file names #######

	# Keep track of progress and time
	RUN=`expr $RUN + 1`
	TIME=`date +"%T"`

	# Get the names of the forward and reverse reads
	FILE1=${RAWREADS[$i]}
	FILE2=${RAWREADS[$i+1]}

    # Create Unique Prefix using first two columns
	PAIRID=`echo $FILE1 | awk '{ split($0, array, "_"); print array[1] }'`

	# Note progress and check file names
	echo -e "\e[0;34m Beginning to run Spotying on read Pair: $PAIRID ---> $RUN of $NPAIRS... \e[0m""	"$TIME
	echo "	$FILE1"
	echo "	$FILE2"

    # Run Spotyping
    python $SPOTYPING $FILE1 $FILE2 --output=Spotyping_temp.txt

    # Store the results
    cat Spotyping_temp.txt >> $RESULTS
    
    # Remove temporary files
    rm Spotyping_temp.txt*
	rm SITVIT_ONLINE*


	echo -e "\e[0;34m Completed running Spotyping on current pair of FASTQ files: $PAIRID ---> $RUN of $NPAIRS. \e[0m"
done

echo -e "\e[0;32m Complete. \e[0m"
