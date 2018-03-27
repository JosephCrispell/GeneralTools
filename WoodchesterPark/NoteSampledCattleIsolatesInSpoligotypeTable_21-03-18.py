#!usr/usr/bin/python

###################
# Import packages #
###################

import sys # Command line arguments
import re # pattern matching

# Identify the sampled isolates from spoligotype table - create new table just for them
# Author: Joseph Crispell
# Created: 21-03-18

# Command line structure:
# python NoteSampledCattleIsolatesInSpoligotypeTable.py sequencedCattleIsolateTable.csv SpoligotypeInfo.csv

#############
# FUNCTIONS #
#############

def getSequencedSampleIDsAndEartags(sequencedCattleFile):

	# Initialise a dictionary to store the isolate IDs and associated eartags
	sampled = {}

	# Open the file and read it line by line
	with open(sequencedCattleFile) as fileLines:
		
		for line in fileLines:
			
			# Remove new line character
			line = line.rstrip()
			
			# Skip header
			if re.match(r"^Sample", line):
				continue
			
			# Skip isolates not sequenced
			if not re.match(r"(.*)yes", line):
				continue
			
			# Split the current line into its columns
			columns = re.split(",", line)
			
			# Create an ID that'll match the fastq files
			id = columns[0]
			if not re.match(r"^AF", columns[0]):
				parts = re.split("/", columns[0])
				id = "HI-" + parts[0] + "-" + parts[1] + "-" + parts[2]
			
			# Store the isolate ID and its eartag and ID to match FASTQ file
			sampled[columns[0]] = [columns[1], id]
	
	return sampled

def readAndPrintSpoligotypeInfoForSampled(spoligotypeFile, sampled, outputFile):
	
	# Open the output file
	output = open(outputFile, "w")
	
	# Open the file and read it line by line
	with open(spoligotypeFile) as fileLines:
		
		for line in fileLines:
			
			# Remove new line character
			line = line.rstrip()
			
			# Skip header
			if re.match(r"^Sample", line):
				output.write(line + ",StrainId\n")
				continue
			
			# Split the current line into its columns
			columns = re.split(",", line)
			
			# Split not sampled
			if not columns[0] in sampled:
				continue
			
			# Write the output file
			output.write(line + "," + sampled[columns[0]][1] + "\n")
			
	output.close()	
	
#################################
# Check command line arguements	#
#################################			

if len(sys.argv) != 4: # First argument is the python script. Python is 0 indexed
	print "Requires three file names in input arguments. Number of input arguments = ", len(sys.argv) - 1
	print "Command line structure:"
	print "\tpython NoteSampledCattleIsolatesInSpoligotypeTable.py sequencedCattleIsolateTable.csv SpoligotypeInfo.csv output.csv"
	sys.exit()

#############################
# Get a list of sampled IDs #
#############################

sampled = getSequencedSampleIDsAndEartags(sys.argv[1])

###############################################
# Read and print spoligotype info for sampled #
###############################################

readAndPrintSpoligotypeInfoForSampled(sys.argv[2], sampled, sys.argv[3])