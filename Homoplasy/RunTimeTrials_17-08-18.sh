#!/bin/bash

# Running time trials to compare HomoplasyFinder to Treetime's homoplasy_scanner
# Author: Joseph Crispell

# Command Line Structure:
# bash RunJavaAndTreetimeTrials_11-05-18.sh

# RUN SCRIPT USING bash script.sh (in Directory containing Raw Reads)
# Note create permissions first ---> chmod +x script.sh

####################
# DEFINE FUNCTIONS #
####################

####################################################################
# Run time trials to compare HomoplasyFinder and homoplasy_scanner #
####################################################################

# Initialise an output file
OUTPUT="timeTaken_17-08-18.csv"
echo "NSequences	Replicate	NSites	Java	R	Phangorn	TreeTime" > $OUTPUT

# Set the date the datasets were created
DATE="17-08-18"

# Loop through the number of sequences in each dataset
for (( NSEQUENCES=100; NSEQUENCES<=1000; NSEQUENCES+=50 ))
do
	TIME=`date +"%T"`
	echo -e "\e[0;34m Processing example data set with $NSEQUENCES sequences... \e[0m""	"$TIME
	
	# Loop through the replicates for each dataset
	for (( REPLICATE=1; REPLICATE<=10; REPLICATE+=1 ))
	do
		
		echo "    Replicate: $REPLICATE"
		
		# Build fasta file
		FASTA="Example_"$((2 * $NSEQUENCES))"-"$NSEQUENCES"_"$REPLICATE"_"$DATE".fasta"
		
		# Get the number of sites in the fasta file - calculate rescale for homoplasy_scanner
		NSITES=`head -n 1 $FASTA | awk '{split($0, array, " "); print array[2]}'`
		RESCALE=`perl -E "say 1 / $NSITES"`
		
		# Build tree file
		TREE="Example_"$((2 * $NSEQUENCES))"-"$NSEQUENCES"_"$REPLICATE"_"$DATE".tree"
		
		# Run homoplasy_scanner
		START="$(date -u +%s.%N)"
		homoplasy_scanner.py --aln $FASTA --tree $TREE --rescale $RESCALE
		END="$(date -u +%s.%N)"
		TREETIME="$(bc <<<"$END-$START")"

		# Run HomoplasyFinder Java tool
		START="$(date -u +%s.%N)"
		java -jar ../HomoplasyFinder.jar 0 $FASTA $TREE
		END="$(date -u +%s.%N)"
		JAVATIME="$(bc <<<"$END-$START")"
		
		# Run HomoplasyFinder in R
		START="$(date -u +%s.%N)"
		Rscript ../../GeneralTools/Homoplasy/SpeedTestConsistencyIndexInHomoplasyFinder_17-08-18.R ./ $TREE $FASTA
		END="$(date -u +%s.%N)"
		RTIME="$(bc <<<"$END-$START")"
		
		# Run phangorn in R
		START="$(date -u +%s.%N)"
		Rscript ../../GeneralTools/Homoplasy/SpeedTestConsistencyIndexInPhangorn_10-08-18.R . $TREE $FASTA
		END="$(date -u +%s.%N)"
		PHANGORN="$(bc <<<"$END-$START")"
		
		# Store the times
		echo $NSEQUENCES"	"$REPLICATE"	"$NSITES"	"$JAVATIME"	"$RTIME"	"$PHANGORN"	"$TREETIME >> $OUTPUT	
		
	done
done