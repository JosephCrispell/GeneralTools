#!usr/usr/bin/python

###################
# Import packages #
###################

import sys # Command line arguments
import re # pattern matching
import os # Get files in current directory
import shutil # Move a file

# Get the right forward and reverse FASTQ files for the TVR sequences
# Author: Joseph Crispell
# Created: 26-03-19

# Command line structure:
# python CountNSampledBadgersInMerged.py merged.txt badgerSamplingInfo.csv

#############
# FUNCTIONS #
#############

def noteEachSamplesFASTQFiles():
	
	# Get all the FASTQ files in the current directory
	filesInDirectory = os.listdir(".")

	# Initialise an empty dictionary
	isolates = {}

	# Initialise an array to store the isolate IDs from the files
	for fileName in filesInDirectory:
	
		# Skip non FASTQ files and reverse files
		if not re.search(pattern=r"(.*)fastq.gz", string=fileName):
			continue
	
		# Split the file name into its columnds
		parts = fileName.split("_")
	
		# Note the isolate ID
		id = parts[0]
	
		# Check if encountered ID before
		if id in isolates:
			isolates[id].append(fileName)
		else:
			isolates[id] = [fileName]

	return isolates

def findSamplesWithTooManyFASTQFiles(sampleFASTQFiles):

	# Initialise an array to store the samples with too many FASTQ files
	samples = []

	# Examine the FASTQ files associated with each sample
	for id in sampleFASTQFiles:

		# Skip samples with 2 FASTQ files
		if len(sampleFASTQFiles[id]) == 2:
			continue

		# Store the current sample ID
		samples.append(id)

	return samples

def noteWhichFASTQFilesToIgnoreFromSamplesWithTooMany(sampleFASTQFiles, samplesWithTooManyFASTQs):

	# Initialise an array to store file names to ignore
	filesToIgnore = []

	# Examine each of the samples with too many FASTQ files
	for id in samplesWithTooManyFASTQs:

		# Get the FASTQ files for the current sample
		fastqs = sampleFASTQFiles[id]

		# Initialise a dictionary to store the FASTQ files IDs
		fileIdentifiers = {}

		# Examine each of the FASTQ files
		for fastq in fastqs:

			# Split the file name into its columnds
			parts = fastq.split("_")

			# Create a file identifier
			fastqIdentifier = parts[0] + "_" + parts[1]

			# Check if encountered this identifier before
			if fastqIdentifier in fileIdentifiers:
				fileIdentifiers[fastqIdentifier].append(fastq)
			else:
				fileIdentifiers[fastqIdentifier] = [fastq]

		# Examine the dictionary of FASTQ file IDs
		for fastqIdentifier in fileIdentifiers:

			# Check if less than two files found with identifier
			if len(fileIdentifiers[fastqIdentifier]) == 1:
				filesToIgnore.append(fileIdentifiers[fastqIdentifier][0])
	
	return filesToIgnore

def moveFilesToDirectory(files, directoryName):

	# Create the directory unless it already exists
	if os.path.exists(directoryName) == False:
		os.mkdir(directoryName)

	# Examine each file
	for file in files:

		# Move the file into the directory
		shutil.move(file, directoryName + file)

		# Print progress
		print "Moved " + file + " to " + directoryName

##########################################
# Find samples with too many FASTQ files #
##########################################

# Note FASTQ files for each sample 
sampleFASTQFiles = noteEachSamplesFASTQFiles()

# Identify the samples with too many FASTQ files
samplesWithTooManyFASTQs = findSamplesWithTooManyFASTQFiles(sampleFASTQFiles)

# Note the FASTQ files (from the samples with too many) to ignore
filesToIgnore = noteWhichFASTQFilesToIgnoreFromSamplesWithTooMany(sampleFASTQFiles, samplesWithTooManyFASTQs)

# Move the FASTQ files being ignored into a new directory
moveFilesToDirectory(filesToIgnore, "ExtraFASTQFilesIgnored/")
