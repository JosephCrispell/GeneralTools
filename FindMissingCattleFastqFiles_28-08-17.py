#!usr/usr/bin/python

###################
# Import packages #
###################

import sys # Command line arguments
import re # pattern matching
import os # Get files in current directory

# Script to quickly find the cattle fastq files (forward and reverse) that are missing
# Author: Joseph Crispell
# Created: 24-08-17

# Command line structure:
# python FindMissingCattleFastqFiles.py mapping.txt

################################
# Check command line arguments #
################################

if len(sys.argv) != 2: # First argument is the python script. Python is 0 indexed
	print "Requires input file. Number of input arguments = ", len(sys.argv) - 1
	sys.exit()

#######################################################
# Check the IDs of the files in the current directory #
#######################################################

# Get the files in the current directory
filesInDirectory = os.listdir(".")

# Initialise an empty dictionary
isolates = {}

# Initialise an array to store the isolate IDs from the files
for fileName in filesInDirectory:
	
	# Skip non FASTQ files and reverse files
	if not re.search(pattern=r"(.*)fastq.gz", string=fileName) or re.search(pattern=r"(.*)_R2(.*)", string=fileName):
		continue
	
	# Split the file name into its columnds
	parts = fileName.split("_")
	
	# Create the isolate ID - column1_column2
	id = parts[0] + "_" + parts[1]
	
	# Store the isolate ID in the dictionary
	isolates[id] = 1

	
############################################################################
# Read in the Isolate mapping file check which isolate's files are missing #
############################################################################

# Get the mapping file from the command line and open it
mappingFile = sys.argv[1]

# Read the file in line by line
with open(mappingFile) as fileLines:
	for line in fileLines:
		
		# Skip the header line
		if re.search(pattern=r"Isolate(.*)", string=line):
			continue
		
		# Remove the end of line character
		line = line.rstrip()
		
		# Split the current line into columns
		columns = line.split("\t")
		
		# Check if isolate ID - first column - was found in the files in the current directory
		if not columns[0] in isolates.keys():
			print "Files not present for isolate:", columns[0]
		
		
		
