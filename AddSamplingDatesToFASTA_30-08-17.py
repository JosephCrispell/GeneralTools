#!usr/usr/bin/python

###################
# Import packages #
###################

import sys # Command line arguments
import re # pattern matching

# Script to add sampling dates to FASTA file
# Author: Joseph Crispell
# Created: 30-08-17

# Command line structure:
# python AddSamplingDatesToFASTA.py

#############
# FUNCTIONS #
#############

def readSamplingInformation(fileName, headerPattern, idColumn, dateColumn, isolateSamplingDates):
	
	# Read the file in line by line
	with open(fileName) as fileLines:
		for line in fileLines:
	
			# Skip the header line
			if re.search(pattern=headerPattern, string=line):
				continue
		
			# Remove the end of line character
			line = line.rstrip()
		
			# Split the current line into its columns
			cols = line.split(",")
		
			# Store the sampling date for the current isolate
			isolateSamplingDates[cols[idColumn]] = cols[dateColumn]
	
	return isolateSamplingDates

################################
# Check command line arguments #
################################

if len(sys.argv) != 4: # First argument is the python script. Python is 0 indexed
	print "Requires three input files. Number of input arguments = ", len(sys.argv) - 1
	sys.exit()
	
###############################################
# Read in badger isolate sampling information #
###############################################

# Initialise a list to store the sampling dates for each isolate
isolateSamplingDates = {}

# Read the sampling information file and store the sampling dates
isolateSamplingDates = readSamplingInformation(fileName=sys.argv[1], headerPattern=r"WB_id(.*)", idColumn=0, dateColumn=4, 
isolateSamplingDates=isolateSamplingDates)
		
###############################################
# Read in cattle isolate sampling information #
###############################################

# Add the cattle isolate sampling dates to the sampling date list
isolateSamplingDates = isolateSamplingDates = readSamplingInformation(fileName=sys.argv[2], headerPattern=r"StringId(.*)", idColumn=51, dateColumn=3, 
isolateSamplingDates=isolateSamplingDates)

#######################
# Read the FASTA file #
#######################

# Read the file in line by line
lineNo = 0
with open(sys.argv[3]) as fileLines:

	for line in fileLines:
		
		# Increment the line number counter
		lineNo += 1
		
		# Skip the header
		if lineNo == 1:
			continue
		
		# Remove the end of line character
		line = line.rstrip()
	
		# Check if current line is sequence  - print and skip
		if not re.search(pattern=r">(.*)", string=line):
			print line
			continue
	
		# Check if reached reference sequence - if so break - not needed
		if re.search(pattern=r">Ref-(.*)", string=line):
			break
	
		# Get isolate ID from sequence
		id = line.split("_")[0][1:] # [1:] is substring function which removes first character ">"
		
		# Check if we can find date
		if id in isolateSamplingDates.keys():
			output = ">" + id + "_" + isolateSamplingDates[id]
			print output
		else:
			print "ERROR: Couldn't find date for isolate ", id



