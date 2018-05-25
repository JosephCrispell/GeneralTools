#!usr/usr/bin/python

###################
# Import packages #
###################

import sys # Command line arguments
import re # pattern matching

# Script to transfer any elements (identified by locus tag) in text file into embl filels
# Author: Joseph Crispell
# Created: 31-08-17

# Command line structure:
# python CopyElementsIntoEMBL.py elements.txt codonTable from.embl to.embl output.embl

#############
# FUNCTIONS #
#############

def getReverseComplement(sequence):
	
	# Get the reverse
	reverse = sequence[::-1]
	
	# Note the nucleotide compliments
	complimentNucleotides = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
	
	# Get the compliment
	compliment = ""
	for nucleotide in list(reverse):
		compliment += complimentNucleotides[nucleotide]
		
	return compliment

def buildCodonTable(fileName):
	
	# Intialise a dictionary to store the codon table
	codonTable = {}
	
	# Read the file in line by line
	with open(fileName) as fileLines:
	
		for line in fileLines:
			
			# Skip header region
			if re.match(r"Codon(.*)", line):
				continue
				
			# Remove the end of line character
			line = line.rstrip()
		
			# Split the current line into parts
			parts = re.split("\t", line)

			# Store the codon information in the codon table
			codonTable[parts[0]] = parts[3]
	
	return codonTable

def takeAway(vector, number):
	
	output = [0] * len(vector)
	
	for index in range(0, len(vector)):
		
		output[index] = vector[index] - number
		
	return output

def convertToAA(sequence, codonTable, key):
	
	# Initialise a string to store the AA sequence
	aaSequence = ""
	
	# Note the start of each codon
	codonStarts = range(0, len(sequence), 3)
	
	# Remove last start if not divisible by 3
	if len(sequence) % 3 != 0:
		del codonStarts[-1]
		print "ERROR: Length of feature (coordinates: " + key + ") not divisble by 3"
	
	# Examine each codon in the sequence
	for i in codonStarts:
		
		# Get the current codon
		codon = sequence[i:i+3] # splice: start:end:by[optional]
		
		# Skip stop codon
		if codonTable[codon] == "*":
			continue
		
		# Get the current codon's AA
		aaSequence += codonTable[codon]
	
	return aaSequence

def getFeatureSequence(coords, genome, compliment):
	
	# Take 1 from the coordinates to make them 0-based
	coords = takeAway(coords, 1)
	
	# Initialise a variable to store the sequence of the feature
	sequence = sequence = genome[coords[0]:coords[1]+1:1].upper() # splice syntax: start:end:by
	
	# Check if want the compliment
	if compliment == True:

		sequence = getReverseComplement(sequence)

	return sequence

def addAASequencesToFeatures(keys, features, sequence, codonTable):
	
	# Examine each of the features
	for key in keys:
		
		# Get the coordinate information from the coords
		coords = parseCoords(key)
		
		# Get the feature information for the current feature
		featureInfo = features[key]
		
		# Skip features that aren't coding sequences
		if re.search(r"FT( +)CDS", featureInfo) == None:
			continue
			
		# Skip pseudo genes
		if re.search(r"/pseudo", featureInfo):
			continue
		
		# Get the nucleotide sequence for the current sequence
		featureSequence = getFeatureSequence(coords, sequence, re.search(r"complement(.*)", key) != None)
		
		# Convert the sequence to an amino acid sequence
		aaSequence = convertToAA(featureSequence, codonTable, key)
		
		# Add the sequence into the feature Info
		linesToAdd = "FT                   /translation=\"" + aaSequence + "\""
		features[key] += linesToAdd + "\n"
		
		print "Added AA sequence to feature (coordinates: " + key + ")"
	
	return features

def getCoordsFromFeatureInfo(featureInfo):
	
	# Initilaise a variable to store the coordinates
	coords = "NA"
	
	# Examine each line
	for line in re.split("\n", featureInfo):
		
		# Split the current line into its columns
		parts = re.split("  +", line)
		
		# Check if coords in the current line
		if re.search("\.\.", line):
			
			# Get the coordinates and finish
			coords = parts[2]
			break
	
	return coords

def findNextEmptyIndex(lines):
	index = -1
	for i in range(0, len(lines)):
		if lines[i] == None:
			index = i
			break
			
	return index

def addLine(lines, line):

	# Check if lines array still has space
	if lines[len(lines) - 1] == None:
		
		# Get the next free index
		nextIndex = findNextEmptyIndex(lines)
		
		# Store the line
		lines[nextIndex] = line

	else:
		
		# Move each element back a space
		for index in range(0, len(lines) - 1):
			lines[index] = lines[index + 1]
		
		# Add in line
		lines[len(lines) - 1] = line
	
	return lines

def getGeneInfo(lines):

	# Initialise an output string
	output = ""
	
	# Look at the previous lines available
	for index in range(len(lines) - 2, 0, -1):
	
		# Skip empty lines
		if lines[index] == None:
			continue
	
		# Split the current line into parts
		parts = re.split("  +", lines[index])
		
		# Check if reached another feature and if it has a gene tag
		if len(parts) == 3 and re.search(r"FT   gene(.*)", lines[index]):
			
			# Store the lines to come as string
			for i in range(index, len(lines)):
				
				# Don't add empty lines
				if not lines[i] == None:
					output += lines[i] + "\n"
			
			# Stop the search
			break
		
		# If reached another feature and it hasn't got a gene tag break out
		elif len(parts) == 3:
			break

	return output

def getFeatureInfoFromEMBL(fileName, tags, codonTable):
	
	# Initialise a dictionary to store the feature information
	features = {}
	
	# Initialise an array to store the feature keys
	keys = []
	
	# Initialise a variable to record that a feature was found
	foundFeature = False
	
	# Initialise a variable to record the sequence
	sequence = None
	
	# Initialise an array to store the previous 100 file lines
	previousLines = [None] * 100
	
	# Read the file in line by line
	with open(fileName) as fileLines:
	
		for line in fileLines:
			
			# Skip header region
			if re.match(r"FT(.*)", line) == None and re.match(r"SQ(.*)", line) == None and sequence == None:
				continue
				
			# Check if reached sequence
			elif re.match(r"SQ(.*)", line):
				
				# Check if we need to store a feature of interest
				if foundFeature == True:
				
					# Get the features coords from the featureInfo
					coords = getCoordsFromFeatureInfo(featureInfo)

					# Store the feature information
					keys.append(coords)

					# Store the information for the previous feature
					features[coords] = featureInfo
				
					# Reset search for features
					foundFeature = False
					
				# Note that reached sequence and skip line
				sequence = ""
				continue

					
			# Remove the end of line character
			line = line.rstrip()
		
			# Store current line
			previousLines = addLine(previousLines, line)
		
			# Split the current line into parts
			parts = re.split("  +", line)		
						
			# If haven't found a feature of interest yet - search for one!
			if foundFeature == False and sequence == None:
				
				# Check whether of the tags was found in the current line
				for tag in tags:
					if re.search(r"/locus_tag=\"" + tag, line, re.IGNORECASE) or re.search(r"/gene=\"" + tag, line, re.IGNORECASE):
					
						# Get the feature info collected from the previous lines so far for the current feature
						featureInfo = getGeneInfo(previousLines)
						foundFeature = True
						break
						
						
			# If found feature - continue building feature information until reach next gene tag
			elif foundFeature == True and parts[1] != "gene":
				featureInfo += line + "\n"
			
			# If found feature and found gene tag - store all the feature information collected
			elif foundFeature == True and parts[1] == "gene":
				
				# Get the features coords from the featureInfo
				coords = getCoordsFromFeatureInfo(featureInfo)

				# Store the feature information
				keys.append(coords)
					
				# Store the information for the previous feature
				features[coords] = featureInfo
				
				# Reset search for features
				foundFeature = False
				
			# Check if we have reached sequence
			elif sequence != None and re.match("//", line) == None:
				
				sequence += parts[1]
	
	# Remove spaces from the nucleotide sequence
	sequence = sequence.replace(" ", "")

	# Note how many of the features were found
	print "Found " + str(len(keys)) + " features of " + str(len(tags))
	
	# Add the protein sequences to the information for the features of interest
	features = addAASequencesToFeatures(keys, features, sequence, codonTable)
					
	# Return an output that has the feature dictionary and keys
	output = {"Keys":keys, "Feature Info":features}
	
	return output

def parseCoords(coords):

	output = coords

	# Check if coordinates have complement parenthesis - remove if present
	if re.search(r"complement(.*)", coords):
		output = coords.split("(")[1]
		output = output[0:(len(output) - 1)]
	
	# Remove ">" if present
	output = output.replace(">", "")
	
	# Split the coordinates into start and end
	output = output.split("..")
	
	# Convert the coordinates to numbers
	output = [int(value) for value in output]

	return output
	
def checkIfFeatureFallsBetweenCoords(currentCoords, previousCoords, features, featuresAdded, orderedKeys, outputFile):

	# Convert the current and previous coords to numbers
	currentCoords = parseCoords(currentCoords)
	previousCoords = parseCoords(previousCoords)
	
	# Calculate the difference between the end of the previous and start of the current
	difference = currentCoords[0] - previousCoords[1]
	
	if difference > 0:
		
		# Compare these coordinates to those of the features
		for key in orderedKeys:
			
			# Skip if feature already added - i.e. have found a feature that it starts after and one that it ends before
			if key in featuresAdded.keys() and featuresAdded[key][0] == True and featuresAdded[key][1] == True:
				continue
				
			# Convert the coordinates to numbers for feature
			coords = parseCoords(key)
		
			# Check if the feature starts after the previous feature's end and haven't already noted this
			if not key in featuresAdded.keys() and coords[0] >= previousCoords[1]:
				featuresAdded[key] = [True, False, None] # Found feature that starts after, Found feature that ends before, Notes
		
			# Check if the current feature ends before the end of the current coords
			if key in featuresAdded.keys() and coords[1] <= currentCoords[1]:
				
				# Add the feature information into the output file
				output.write(features[key])
			
				# Record that inserted feature
				featuresAdded[key][1] = True
				featuresAdded[key][2] = "Previous = " + str(previousCoords[0]) + ":" + str(previousCoords[1]) + "   Current = " + str(currentCoords[0]) + ":" + str(currentCoords[1]) + "\tFeature = " + str(coords[0]) + ":" + str(coords[1])
	
	return featuresAdded

def readAndPrintEMBLInsertingFeatures(EMBL, output, features, orderedKeys):
	
	# Initialise a variable to store the coords of the previous feature	
	previousCoords = None

	# Initialise a variable to record when we've found the sequence block
	foundSequence = False
	
	# Initialise a dictionary to check whether feature added
	featuresAdded = {}
	
	# Read in the file line by line
	with open(EMBL) as fileLines:
	
		for line in fileLines:
		
			# Remove the end of line character
			line = line.rstrip()
		
			# Check the format of the locus tag in the line - Damien's advice
			if re.search(r"/locus_tag=\"MB(.*)\"", line):
				line = line.replace("\"MB", "\"Mb")
			if re.search(r"/locus_tag=\"(.*)C\"", line):
				line = line.replace("C\"", "c\"")
		
			# Check if found sequence block
			if foundSequence == False and line.startswith("SQ") == True:
				foundSequence = True
			
			# Print line to file and skip if found sequence block
			if foundSequence == True:
				output.write(line + "\n")
				continue
			
			# Split the current line into parts
			parts = re.split("  +", line)
		
			# Skip lines that aren't a feature header
			if len(parts) != 3 or re.search(r"FH(.*)Key(.*)Location/Qualifiers", line):
			
				output.write(line + "\n")
				continue
			
			# Get the coordinates from the current line
			coords = parts[2]
			
			# Check if a feature fell between the current and previous coords
			if previousCoords != None and previousCoords != coords:
				featuresAdded = checkIfFeatureFallsBetweenCoords(coords, previousCoords, features, featuresAdded, orderedKeys, output)
		
			# Set the current coords to previous coords
			previousCoords = coords
			
			# Print the current line
			output.write(line + "\n")
		
	# Close the output file
	output.close()
	
	return featuresAdded

def checkIfFeaturesNotInserted(features, featuresAdded):
	
	print "Found " + str(len(features)) + " features and inserted " + str(len(featuresAdded)) + " of them."
	
	if len(features) > len(featuresAdded):
	
		# Find the features not added
		for key in features.keys():
			
			# Check if feature added
			if not key in featuresAdded.keys():
				index = index + 1
			
				print "\n-----------------------------------------------------------------------------"
				print features[key]

def checkFeaturesAdded(featuresAdded):
	
	# Initialise a count variable
	count = 0
	
	# Initialise a hashtable to store the unique locus_tags
	locusTags = {}
	
	# Examine each feature
	for key in features.keys():
		
		# Check if locus tag present and skip
		if(re.search(r"(.*)/locus_tag(.*)", features[key])):
			
			# Split the info into lines
			lines = re.split("\n", features[key])
			for line in lines:
				
				if(re.search(r"(.*)/locus_tag(.*)", line)):
					tag = re.split("=", features[key])[2]
					locusTags[tag] = 1
					break
					
			continue
		
		count = count + 1
		#print "\n-----------------------------------------------------------------------------"
		#print features[key]
	
	print "Found " + str(count) + " features without /locus_tag."
	print "Found " + str(len(locusTags)) + " unique locus tags."

def getElementTags(fileName):
	
	# Initialise a dictionary to store the feature tags
	tags = []
	
	# Read the file in line by line
	with open(fileName) as fileLines:
	
		for line in fileLines:
			
			# Skip header region
			if re.search(r"locus/gene(.*)", line):
				continue
			
			# Remove the end of line character
			line = line.rstrip()
		
			# Store the current tag
			tags.append(line)
	
	return tags

#################################
# Check command line arguements	#
#################################			

#if len(sys.argv) != 4: # First argument is the python script. Python is 0 indexed
#	print "Requires three file names in input arguments. Number of input arguments = ", len(sys.argv) - 1
#	print "Command line structure:"
#	print "\tpython CopyElementsIntoEMBL.py CopyElementsIntoEMBL.py elements.txt from.embl to.embl output.embl"
#	sys.exit()

############################################
# Get the tags of the elements of interest #
############################################

elementsTable = sys.argv[1]
tags = getElementTags(elementsTable)

###########################
# Read in the codon table #
###########################

# Get the codon table file from the command line
codonTableFile = sys.argv[2]

# Build the codon table
codonTable = buildCodonTable(codonTableFile)

#################################################
# Get feature information from input EMBL file #
#################################################
	
# Get the input file name from the command line
inputEMBL = sys.argv[3]

# Read the file and get the information for the elements of interest
featuresFound = getFeatureInfoFromEMBL(inputEMBL, tags, codonTable)
features = featuresFound["Feature Info"]
orderedKeys = featuresFound["Keys"] # Ordered by their insertion into feature table - i.e. their order in EMBL

###############################################################
# Read in EMBL file and insert features to create output EMBL #
###############################################################

# Get the input file name from the command line
EMBL = sys.argv[4]

# Open the output file
output = open(sys.argv[5],"w")

# Read in the EMBL, print out its lines and insert features at the correct locations to produce the output EMBL file
featuresAdded = readAndPrintEMBLInsertingFeatures(EMBL, output, features, orderedKeys)

# Check that all features were inserted
checkIfFeaturesNotInserted(features, featuresAdded)

# Have another look at the features added
checkFeaturesAdded(featuresAdded)