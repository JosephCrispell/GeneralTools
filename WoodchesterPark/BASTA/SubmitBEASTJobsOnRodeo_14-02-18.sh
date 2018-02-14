#!/bin/sh
#!/bin/bash

# Submit BEAST jobs on the Rodeo server
# Author: Joseph Crispell

# RUN SCRIPT USING bash script.sh (in Directory containing Raw Reads)
# Note create permissions first ---> chmod +x script.sh

# Command Line Structure:
# bash SubmitBEASTJobsOnRodeo.sh

#############
# FUNCTIONS #
#############

####################
# Find Directories #
####################

# Get a list of the directories in teh current directory
DIRECTORIES=(`ls -d */`)

#######################################################
# Run the XML file present in each directory in BEAST #
#######################################################

# Create a file to store the PIDs of the jobs submitted
DATE=`date +"%d-%m-%y"`
PIDFILE="PIDs_"$DATE".txt"
echo "FileName	PID" > $PIDFILE

for DIRECTORY in "${DIRECTORIES[@]}"
do
	# Get the name of the XML file present in the current directory without suffix
	XML=`ls $DIRECTORY | grep ".xml" | awk '{ split($0, array, "."); print(array[1]) }'`
	
	# Move into current directory
	cd $DIRECTORY
	
	# Run that XML in BEAST and note the PID
	nohup beast -threads 1 $XML".xml" > $XML"_ScreenOutput.txt" &
	
	# Wait a second to let java initialise
	sleep 1
	
	# Get the PID of tyhe submitted job
	PID=`ps -u | grep "java" | grep $XML | awk '{ split($0, array, " +"); print(array[2]) }'`
	
	# Move out of current directory
	cd ..
		
	# Print it to file along with the XML file name
	echo $XML"	"$PID >> $PIDFILE
	
done



