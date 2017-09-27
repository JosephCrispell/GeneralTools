#!usr/usr/bin/python

###################
# Import packages #
###################

import sys # Command line arguments
import re # pattern matching

# Script to transfer "repeat_region"s and "mobile_element"s from one EMBL file to another
# Author: Joseph Crispell
# Created: 31-08-17

# Command line structure:
# python CopyElementsIntoEMBL.py from.embl to.embl output.embl

#############
# FUNCTIONS #
#############

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
	
		# Split the current line into parts
		parts = re.split("  +", lines[index])
		
		# Check if reached another feature and if it has a gene tag
		if len(parts) == 3 and re.search(r"FT   gene(.*)", lines[index]):
			
			# Store the lines to come as string
			for i in range(index, len(lines)):
				output += lines[i] + "\n"
			
			# Stop the search
			break
		
		# If reached another feature and it hasn't got a gene tag break out
		elif len(parts) == 3:
			break

	return output

def getFeatureInfoFromEMBL(fileName):
	
	# Initialise a dictionary to store the feature information
	features = {}
	
	# Initialise an array to store the feature keys
	keys = []
	
	# Initialise a variable to record that a feature was found
	foundFeature = False
	
	# Initialise an array to store the previous 100 file lines
	previousLines = [None] * 100
	
	# Read the file in line by line
	with open(fileName) as fileLines:
	
		for line in fileLines:
			
			# Remove the end of line character
			line = line.rstrip()
		
			# Store current line
			previousLines = addLine(previousLines, line)
		
			# Split the current line into parts
			parts = re.split("  +", line)
		
			# Skip header region
			if not re.search(r"FT(.*)", line):
				
				# If reached sequence - break out
				if re.search(r"SQ(.*)", line):
					break
				else:
					continue
			
			# Check if reached annotation for repeat_region or mobile_element
			if re.search(r"FT   repeat_region(.*)", line) or re.search(r"FT   mobile_element(.*)", line):
				
				# Check if already found Feature
				if foundFeature == True:
					
					# Record the key (coords)
					keys.append(coords)
					
					# Store the information for the previous feature
					features[coords] = featureInfo
					
					# Reset the feature info
					featureInfo = line + "\n"
				else:
					
					# Note that found a feature
					foundFeature = True
				
					# Get the gene information from the previous lines
					featureInfo = getGeneInfo(previousLines)
					
					# Check if found gene - if not add current line
					if featureInfo == "":
						featureInfo = line + "\n"
				
				# Get the current features coordinates - use as key
				coords = parts[2]
			
			# Keep adding the next lines to feature Info until reach next feature
			elif foundFeature == True and len(parts) != 3:
				featureInfo += line + "\n"
			
			# Store the feature info when reached next feature
			elif foundFeature == True and (len(parts) == 3 or re.search(r"SQ(.*)", line)):
			
				# Record the key (coords)
				keys.append(coords)
			
				# Store the feature info
				features[coords] = featureInfo
				
				# Reset found feature variable
				foundFeature = False

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

#################################
# Check command line arguements	#
#################################			

if len(sys.argv) != 4: # First argument is the python script. Python is 0 indexed
	print "Requires three file names in input arguments. Number of input arguments = ", len(sys.argv) - 1
	print "Command line structure:"
	print "\tpython CopyElementsIntoEMBL.py from.embl to.embl output.embl"
	sys.exit()

#################################################
# Get feature information from input EMBL file #
#################################################
	
# Get the input file name from the command line
inputEMBL = sys.argv[1]

# Read the file and get the information for the "repeat_region"s and "mobile_element"s
featuresFound = getFeatureInfoFromEMBL(inputEMBL)
features = featuresFound["Feature Info"]
orderedKeys = featuresFound["Keys"] # Ordered by their insertion into feature table - i.e. their order in EMBL

###############################################################
# Read in EMBL file and insert features to create output EMBL #
###############################################################

# Get the input file name from the command line
EMBL = sys.argv[2]

# Open the output file
output = open(sys.argv[3],"w")

# Read in the EMBL, print out its lines and insert features at the correct locations to produce the output EMBL file
featuresAdded = readAndPrintEMBLInsertingFeatures(EMBL, output, features, orderedKeys)

# Check that all features were inserted
checkIfFeaturesNotInserted(features, featuresAdded)
