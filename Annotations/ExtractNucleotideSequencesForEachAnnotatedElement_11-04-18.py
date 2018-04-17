#!usr/usr/bin/python

###################
# Import packages #
###################

import sys # Command line arguments
import re # pattern matching
import os # Used to recognise operatibng system for file name retrieval

# Script to extract FASTA sequence from the end of complete genome annotations.gbff files
# Author: Joseph Crispell
# Created: 10-04-18

# Command line structure:
# python python ExtractFASTAFromGenbankAnnotations_10-04-18.py annotations.gbff

#############
# FUNCTIONS #
#############

def getAnnotationFilesInCurrentDirectory():
	
	# Intialise a list to store the file names
	files = []
	
	# Loop through all files and directories in current directory
	for name in os.listdir("."):
		
		# Check if current name is associated with a file
		if os.path.isfile(name) and name.endswith("gbff"):
			files.append(name)
	
	return files

def getAnnotatedElementSequencesFromFile(genbankFile, elementSequences):

	# Initialise a dictionary to store the unique coordinates for a sequence
	coords = {}

	# Initialise a variable to store the sequence for the annotations
	sequence = ""

	# Initialise a variable to record whether reached a sequence section
	foundSequence = False

	# Create a variable to track how many sequences we've looked at
	nSequences = 0

	# Read in the file line by line
	with open(genbankFile) as fileLines:

		for line in fileLines:
	
			# Remove the end of line character
			line = line.rstrip()
		
			# Split the current line into parts
			parts = re.split(" +", line)
		
			# Search for lines with coordinates and store them
			if len(parts) == 3 and re.search(r"\.\.", parts[2]):
				coords[parts[2]] = line
			
			# Check if reached sequence part
			if foundSequence == False and re.search(r"ORIGIN", line):
				foundSequence = True
				nSequences += 1
		
			# If reached sequence section store the sequence from the current line		
			if foundSequence == True and re.search(r"[a-z]", line):
		
				for i in range(2, len(parts)):
					sequence += parts[i].upper()
				
			# Check if reach the end of the sequence
			if foundSequence == True and re.search(r"//", line):
			
				# Get the nucleotide sequences of the elements
				elementSequences = getNucleotideSequences(coords, sequence, elementSequences, nSequences, genbankFile)
			
				# Reset the element and sequence information
				foundSequence = False
				coords = {}
				sequence = ""
				
	return elementSequences

def getNucleotideSequences(coords, sequence, elementSequences, nSequences, fileName):

	# Coordinates can be reverse complement (complement())
	# Coordinates can be join (join()) (AND ALSO COMPLEMENT) - convert to array of number and find range
	
	# Check if looking at first (genome) sequence or subsequent (plasmid) sequences
	tag = "GENOME"
	if nSequences > 1:
		tag = "PLASMID"
	
	# Examine each coordinate
	for coord in coords:
		
		# Get the start and end from coords and check whether compliment
		coordInfo = getStartEnd(coord)
		
		# Extract the nucleotide sequence
		elementSequence = elementSequence = sequence[(coordInfo["start"]-1):(coordInfo["end"]-1)]
		
		# Check if need the reverse compliment
		if coordInfo["compliment"] == True:
		
			elementSequence = reverseCompliment(elementSequence)
			
		# Add the current element's sequence into the dictionary
		elementSequences[fileName + ":" + str(coordInfo["start"]) + ":" + str(coordInfo["end"]) + ":" + tag] = elementSequence # Note use of dictionary will remove duplicate coordinates
		
	return elementSequences

def reverseCompliment(sequence):
	
	# Reverse the sequence
	sequence = sequence[::-1] # [start:end:step] -1 for step and blank for start and end reverses string
	
	# Get the compliment of each nucleotide
	compliments = {"A":"T", "C":"G", "G":"C", "T":"A"}
	complimentSequence = ""
	for nucleotide in sequence:
		
		# Check haven't found spurious nucleotide (found "M"s in GCF_000008585.1_ASM858v1)
		if nucleotide in compliments:
			complimentSequence += compliments[nucleotide]
		else:
			complimentSequence += "N"
			print "Found non-nucleotide character: " + nucleotide
		
	return complimentSequence
	
def getStartEnd(coord):
	
	# Initialise a dictionary to store the start, end, and whether compliment required
	output = {"start":0, "end":0, "compliment":False}
	
	# Remove ">" or "<" if present (note re.sub(a, b, string) recognises regex whereas string.replace() doesn't)
	coord = re.sub(r">|<", "", coord)
	
	# Check whether compliment required
	if re.search(r"^complement(.*)", coord):
		coord = coord[11:(len(coord) - 1)]
		output["compliment"] = True
		
	# Check whether cordinates a join() set
	if re.search(r"^join(.*)", coord):
		
		# Get all the start and end coordinates for each element being joined
		coord = coord[5:(len(coord) - 1)]
		values = convertToInt(re.split(r"\.\.|,", coord))
		
		# Calculate the range - start and end for combined
		startEnd = getRange(values)
		output["start"] = startEnd[0]
		output["end"] = startEnd[1]
		
	else:
		
		# Get the start and end coordinates for the current element
		values = convertToInt(re.split(r"\.\.", coord))
		output["start"] = values[0]
		output["end"] = values[1]
		
	return output

def getRange(array):
	
	# Initialise variables to store the min and max
	min = None
	max = None
	
	# Search for the min and max values of the array
	for value in array:
	
		# Check if found new min
		if min == None or value < min:
			min = value
			
		# Check if found new max
		if max == None or value > max:
			max = value
			
	return [min, max]
		
def convertToInt(array):

	output = [None] * len(array)
	for i in range(0, len(array)):
		output[i] = int(array[i])
		
	return output
			
#################################################################################
# Get the nucleotide sequence of each annotated element in each annotation file #
#################################################################################

# Get a list of the annotation files in the current directory
files = getAnnotationFilesInCurrentDirectory()

# Initialise a dictionary to store the nucleotide sequences of each element in each file. Key: FileName:start:end:GENOME|PLASMID
elementSequences = {}

# Examine each genbank file and get the element sequences
for genbankFile in files:
	
	# Print progress
	print genbankFile
	
	# Get the nucleotide sequence from each annotated element in file
	elementSequences = getAnnotatedElementSequencesFromFile(genbankFile, elementSequences)

############################################################################
# Print each sequence, noting its coordinates, and origin (genome/plasmid) #
############################################################################
