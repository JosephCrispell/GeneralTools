#!/bin/bash

# Combining FASTQ files for isolates sequenced across multiple lanes
# Author: Joseph Crispell

# Command Line Structure:
# bash CombineFASTQsAcrossLanes.sh

# Get the prefix numbers from the directory names
SAMPLENUMBERS=(`ls | grep "L001" | egrep -v "fastq" | awk '{ split($0, array, "_"); print array[1] }'`)

# Examine each sample
for NUMBER in "${SAMPLENUMBERS[@]}"
do
   
   echo "Combining files from run: $NUMBER..."
   
   # Note the directory for the lane 1 fastq files
   DIR=$NUMBER"_L001*"
   
   # Build the output file to store the FORWARD reads
   FILEFORWARD=`ls $DIR | grep "R1" | awk '{ split($0, array, "_"); print array[1]"_"array[2]"_L001-4_"array[4]"_"array[5]}'`
   FILEFORWARD=`echo ${FILEFORWARD:0:-3}` # trims away ".gz"
   
   # Build the output file to store the REVERSE reads
   FILEREVERSE=`ls $DIR | grep "R2" | awk '{ split($0, array, "_"); print array[1]"_"array[2]"_L001-4_"array[4]"_"array[5]}'`
   FILEREVERSE=`echo ${FILEREVERSE:0:-3}`
   
   # Examine each lane
   for LANE in {1..4}
   do
      
	  # Progress
	  echo "  Lane: $LANE"
	  
	  # Note the directory name for the current lane and move into it
	  DIR=$NUMBER"_L00"$LANE"-*"
	  cd $DIR
	  
	  # Unzip all the files for the current lane
	  pigz --decompress *
	  #gunzip *
	  	  
	  # Send the FORWARD reads into the FORWARD output file
	  FILE="*_R1_*"
	  cat $FILE >> "../"$FILEFORWARD
	  
	  # Send the REVERSE reads into the REVERSE output file
	  FILE="*_R2_*"
	  cat $FILE >> "../"$FILEREVERSE
	  
	  # Zip up the original files and exit the directory
	  pigz *
	  #gzip *
	  cd ..
   done
   
   # Progress
   echo "  Zipping up output files..."
   
   # Zip up the output files
   pigz *fastq
   #gzip *fastq
   
done
