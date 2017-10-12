#!usr/usr/bin/python

###################
# Import packages #
###################

import sys # Command line arguments
import re # pattern matching

# Count the number of sampled badgers associated with badger isolates in merged file
# Author: Joseph Crispell
# Created: 09-10-17

# Command line structure:
# python CountNSampledBadgersInMerged.py merged.txt badgerSamplingInfo.csv

#############
# FUNCTIONS #
#############

def getBadgerIsolateIDsFromMerged(mergedFile):

	# Initialise an array to store the badger isolates IDs
	badgerIsolates = []

	# Read the file in line by line
	with open(mergedFile) as fileLines:
		
		for line in fileLines:
	
			# Find the field header line
			if re.search(r"(^)#CHROM", line):
			
				# Get the list of isolates from the current line
				isolates = re.split(":", re.split("\t", line)[2])
			
				# Examine each of the badger isolates
				for i in range(len(isolates)):
				
					# Check if badger isolate
					if re.search(r"(^)WB", isolates[i]):
					
						# Get isolate ID
						id = re.split("_", isolates[i])[0]
						
						# Store the isolate id
						badgerIsolates.append(id)
					
				# Finish looking at merged file
				break
			
	# Print out how mnay badger isolates were found
	print("Found " + str(len(badgerIsolates)) + " badger isolates")
			
	return badgerIsolates

def returnUnique(array):

	# Initialise an array to store the unique elements
	uniqueElements = []
	
	# Examine each element in array
	for i in range(len(array)):
	
		# Add element if not already present in array
		if not array[i] in uniqueElements:
			uniqueElements.append(array[i])
			
	return(uniqueElements)

def countNumberSampledBadgers(samplingInfoFile, badgerIsolatesFromMerged):
	
	# Initialise a vector to store the badger tattoos
	tattoos = []

	# Read the file in line by line
	with open(samplingInfoFile) as fileLines:
		
		for line in fileLines:
		
			# Skip the header line
			if re.search(r"(^)WB_id", line):
				continue
		
			# Split the current line into its columns
			cols = re.split(",", line)
		
			# Check if current isolate in IDs array from merged file
			if cols[0] in badgerIsolates:
				tattoos.append(cols[3])
			
	# Find the unique tattoos
	uniqueTattoos = returnUnique(tattoos)

	print("Found " + str(len(uniqueTattoos)) + " sampled badgers")
	
#################################
# Check command line arguements	#
#################################			

if len(sys.argv) != 3: # First argument is the python script. Python is 0 indexed
	print "Requires two file names in input arguments. Number of input arguments = ", len(sys.argv) - 1
	print "Command line structure:"
	print "\tpython CountNSampledBadgersInMerged.py merged.txt badgerSamplingInfo.csv"
	sys.exit()

####################################################
# Get a list of the badger isolates in merged file #
####################################################
	
# Get the input file name from the command line
mergedFile = sys.argv[1]

# Get a list of the badger isolate IDs from merged file
badgerIsolates = getBadgerIsolateIDsFromMerged(mergedFile)

########################################################
# Count the number of badgers associated with isolates #
########################################################

# Get the input file name from the command line
samplingInfoFile = sys.argv[2]

# Count the number of sampled badgers associated with the isolates found in merged file
countNumberSampledBadgers(samplingInfoFile, badgerIsolates)