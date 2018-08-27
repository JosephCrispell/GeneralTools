#!/bin/bash

# Download FASTQ files from NCBI using SRA toolkit
# Author: Joseph Crispell

# Command Line Structure:
# bash DownloadFASTQsUsingSRAToolkit.sh linkTable.csv pathToFasterq-dump

# Downloaded linkTable.csv from: https://doi.org/10.1371/journal.pone.0189838.s009 (paste content into excel and coverted to csv)

# Requires SRA Toolkit
# Downloaded from here: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software

# Get the command line arguments
ACCESSIONSFILE=$1
FASTERQDUMP=$2

# Get the isolate accession numbers from the file - using fourth column - secondary accession
ACCESSIONS=(`cat $ACCESSIONSFILE | egrep -v "ENAProjectAccesion" | awk '{ split($0, array, ","); print array[4]}'`)
N=${#ACCESSIONS[@]}

# Report how many Accessions found
echo -e "\e[0;34m Found $N Accessions in the file provided ... \e[0m"

# Examine each accession number
COUNT=0
for ACCESSION in "${ACCESSIONS[@]}"
do

	# Keep track of progress and time
	COUNT=`expr $COUNT + 1`
	TIME=`date +"%T"`

	# Check that file hasn't already been downloaded
	FOUND=`ls | grep $ACCESSION | wc -l`
	if [ $FOUND == "2" ]
	then
		echo -e "\e[0;34m Files for $ACCESSION already present ($COUNT of $N)\e[0m""	"$TIME
		continue
	fi

	# Note progress
	echo -e "\e[0;34m Beginning file download for $ACCESSION ($COUNT of $N)\e[0m""	"$TIME

	# Download the FASTQ files for the current ACCESSION
	$FASTERQDUMP $ACCESSION --outdir . --progress --split-files --threads 12

	# Zip up the files downloaded
	echo "Zipping up downloaded files..."
	pigz *fastq

done

# Run fastqc on the output
#TIME=`date +"%T"`
#echo -e "\e[0;34m Finished file downloads. Running FASTQC...\e[0m""	"$TIME
#fastqc *fastq.gz
#rm *fastqc.zip
#mkdir FASTQC
#mv *fastqc.html FASTQC

# FINISHED!!!
TIME=`date +"%T"`
echo -e "\e[0;32m Finished! \e[0m""	"$TIME
