#!usr/usr/bin/python

###################
# Import packages #
###################

import sys
import re # pattern matching

# Python script to count sequenced bacteria names appearing in puplication titles
# Author: Joseph Crispell
# Created: 20-09-17

# Command line structure:
# python CountBacteriaAppearingInScholarSearchFromPoP.py ncbiBacteria.txt PoP.csv

#############
# FUNCTIONS #
#############

def getListOfBacteriaFromNCBIFile(fileName):
	
	print "Reading list of bacteria from NCBI file..."
	
	# Initialise a dictionary to store the bacteria names
	bacteria = {}

	# Read in the file line by line
	with open(fileName) as fileLines:
	
		for line in fileLines:
			
			# Skip the header line
			if line.startswith("#Organism/Name"):
				continue
			
			# Remove the end of line character
			line = line.rstrip()

			# Split the current line into its columns
			cols = re.split("\t", line)
		
			# Split the current bacteria's name into its parts
			parts = re.split(" ", cols[0])
		
			# Build name using only first two parts
			name = parts[0] + " " + parts[1]
		
			# Store the current bacteria's name, if not already present
			if name not in bacteria:
				bacteria[name] = 0

	return bacteria

def countTimesEachBacteriaFoundInPoPResults(fileName, bacteria):
	
	print "Counting number of times each bacteria is found in Publish Or Perish search results..."
	
	# Initialise a line counter
	lineNo = 0
	
	# Read the file in line by line
	with open(popFile) as fileLines:

		for line in fileLines:
		
			lineNo += 1
		
			# Skip the header line
			if line.startswith("Cites"):
				continue
		
			# Remove the end of line character
			line = line.rstrip()
		
			# Split the current line into its columns
			cols = re.split(",", line)
		
			# Get the title of the current search result
			title = cols[2]
		
			# Search for each of the bacteria
			for key in bacteria.keys():
				if re.search(key, title):
					bacteria[key] += 1	
	
			# Progress information
			if lineNo % 100 == 0:
				sys.stdout.write('.')
	
	sys.stdout.write('\n')
	
	
	return bacteria
	
###############
# Check input #
###############

if len(sys.argv) != 3:
	print "Python script to count sequenced bacteria names appearing in puplication titles\n"
	print "Command line structure:"
	print "\tpython CountBacteriaAppearingInScholarSearchFromPoP.py ncbiBacteria.txt PoP.csv"
	print "\t\tncbiBacteria.txt\tPath to file containing information for sequenced bacteria from ncbi website\n"
	print "\t\tPoP.csv\tPath to Publish Or Perish search results file in CSV format"

	sys.exit()

##############################################
# Read in the list of bacteria to search for #
##############################################

# Get the ncbi bacteria file name
ncbiFile = sys.argv[1]

# Get a list of bacteria from the file
bacteria = getListOfBacteriaFromNCBIFile(ncbiFile)

#########################################################################
# Count the number of times each bacteria appears in PoP search results #
#########################################################################

# Get the PoP from the command line
popFile = sys.argv[2]

# Count the number of times each bacteria appear in PoP search result titles
bacteria = countTimesEachBacteriaFoundInPoPResults(popFile, bacteria)
