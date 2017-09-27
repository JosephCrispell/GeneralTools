#!/bin/bash

# Create the nexus and isolate trait files necessary for the tip randomisation test
# Author: Joseph Crispell

# Command Line Structure:
# bash CreateTipRandomisationFiles.sh RandomlyShuffleYearsInXML.pl nTimes date input.xml outPrefix

# RUN SCRIPT USING bash script.sh (in Directory containing Raw Reads)
# Note create permissions first ---> chmod +x script.sh

# Get the input information from the command line
SHUFFLE=$1
NTIMES=$2
DATE=$3
XML=$4
PREFIX=$5

# Do NTIMES random shuffles of isolate years
COUNT=0
while [ $COUNT -lt $NTIMES ]
do
	COUNT=`expr $COUNT + 1`
	echo -e "\e[0;34m Randomly shuffling isolate Years: \e[0m"$COUNT
	
	# Create a directory to store the shuffled data
	DIRECTORY="TipRandomisation_"$DATE"_"$COUNT
	mkdir $DIRECTORY
	
	# Move into new directory
	cd $DIRECTORY
	
	# Randomly shuffle the isolate sampling years - Note "1" specifies that shuffling is to be done between years (0 for normal)
	perl "../"$SHUFFLE 1 "../"$XML $PREFIX"_"$DATE"_"$COUNT".xml"
	
	# Move back into the original directory
	cd ..
	
done

echo -e "\e[0;32m Finished. \e[0m"
