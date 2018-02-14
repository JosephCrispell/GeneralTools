#!usr/usr/bin/python

###################
# Import packages #
###################

import sys # Command line arguments
import re # pattern matching
import random # shuffling array
import os # Creating directory

# Shuffle the years in a BASTA xml file
# Author: Joseph Crispell
# Created: 13-02-18

# Command line structure:
# python ShuffleYearsInBASTAxml.py basta.xml n

#############
# FUNCTIONS #
#############

def getXMLFileLines(xmlFileName):
	# Create a vector to store each file line
	lines = []

	# Read the file in line by line
	with open(xmlFileName) as fileLines:
		
		for line in fileLines:
		
			# Remove end of line character
			line = line.rstrip()
		
			# Store the current line
			lines.append(line)
	
	return lines

def getTipsAndDates(line):
	
	# Extract the tip and date information from the current line
	tipDateInfo = line.split("\"")[1].split(",")
	
	# Initialise a list to store the tip dates
	tipDates = {}
	
	# Examine the tip and date information
	for tipDate in tipDateInfo:
		
		# Split the tip IDs from their dates
		parts = tipDate.split("=")
		
		# Store the tip ID and date
		tipDates[parts[0]] = parts[1]
	
	return(tipDates)

def shuffleDates(tipDates):

	# Get the tip IDs
	tips = tipDates.keys()
	
	# Get the tip dates
	dates = tipDates.values()
	
	# Shuffle the values
	random.shuffle(dates)
	
	# Reassign random dates
	for i in range(0, len(tips)):
		tipDates[tips[i]] = dates[i]
		
	return(tipDates)

def buildTipDatesLine(tipDates):
	
	# Start building output line
	line = "\t\tvalue=\""
	
	# Get the tip ids
	tips = tipDates.keys()
	
	# Add the info for the first tip
	line = line + tips[0] + "=" + tipDates[tips[0]]
	
	# Add the tip IDs and dates
	for i in range(1, len(tips)):
		line = line + "," + tips[i] + "=" + tipDates[tips[i]]
		
	# Finish the line
	line = line + "\">"
	
	return(line)

def createXMLsWithRandomTipDates(xmlFileName, fileLines, nToCreate):
	
	# Create n xml files with random rates
	for i in range(0, nToCreate):
		
		# Initialise a variable to record whether found tip dates section
		foundTipDates = False
		
		# Build an output file name
		parts = xmlFileName.split("/")
		outputFileName = parts[len(parts)-1].split(".")[0] + "_randomTipDates-" + str(i)
	
		# Create a directory for that file
		if not os.path.exists(outputFileName):
			os.makedirs(outputFileName)
	
		# Open the output file
		output = open(outputFileName + "/" + outputFileName + ".xml", "w") 

		# Examine each of the file lines
		for line in fileLines:

			# Find the tip dates section
			if re.search(r"(^)\t<timeTraitSet", line):
				foundTipDates = True
				output.write(line + "\n")
	
			# If found the tip dates section - shuffle the dates
			elif foundTipDates == True:
	
				# Get the dates associated with each tip
				tipDates = getTipsAndDates(line)
			
				# Shuffle the dates
				tipDates = shuffleDates(tipDates)
		
				# Print the shuffled dates into output file
				output.write(buildTipDatesLine(tipDates) + "\n")
		
				# Note that finished with tip dates
				foundTipDates = False
	
			# Search for and change the output file names
			elif re.search(r"fileName=", line):
		
				# Get the current file name
				currentFileName = line.split("\"")[3].split(".")[0]
		
				# Replace current file name witgh new in the current line
				line = line.replace(currentFileName, outputFileName)
		
				# Print the edited line into the output file
				output.write(line + "\n")
		
			else:
				output.write(line + "\n")

		# Close the output file name
		output.close()
	
#################################
# Check command line arguements	#
#################################			

if len(sys.argv) != 3: # First argument is the python script. Python is 0 indexed
	print "Requires two input arguments. Number of input arguments = ", len(sys.argv) - 1
	print "Command line structure:"
	print "\tpython ShuffleYearsInBASTAxml.py basta.xml n"
	sys.exit()

#########################
# Read in the BASTA xml #
#########################

# Get the xml file from the command line
xmlFileName = sys.argv[1]

# Get the xml file lines
fileLines = getXMLFileLines(xmlFileName)

###########################################
# Print out XML files with shuffled years #
###########################################

# Get the number of XML files to create from the command line
nToCreate = int(sys.argv[2])

# Create n XML files with random tip dates
createXMLsWithRandomTipDates(xmlFileName, fileLines, nToCreate)