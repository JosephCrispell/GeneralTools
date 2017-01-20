#!/usr/bin/perl
# Searching files in a directory
use warnings;
use strict;

# Author: Joseph Crisp
# Extract the sampling information needed to accompany the FASTQ files on NCBI

# Command Line Structure
# perl BuildNCBISRAUploadFile_20-01-17.pl BioSampleInfo.csv outputFile.tsv PRJN

# BioSampleInfo.csv file structure:
# 	target_db	status			msg	resp_status	accession		name			tax_id	extras	file_path	file_id	file_url
# 	0			1				2				3				4				5		6		7			8		9
# 	BioSample					processed-ok
#								Successfully loaded
#												processed-ok	SAMN06245910	AgR297	1765


# Output file structure:
# 	bioproject_accession	biosample_accession	library_ID	title	library_strategy	library_source	library_selection	
# 	0 ###					1 ###				2 ###		3 ###	4 ###				5 ###			6 ###
#
# 	library_layout	platform	instrument_model	design_description	filetype	filename	filename2	filename3	filename4
#	7 ###			8 ###		9 ###				10 ###				11 ###		12 ###		13			14			15
#
# Column Definitions:
# bioproject_accession		PRJNA[number] - the number assigned to the BioProject submission
# biosampled_accession		SAMN[number] - the numbers assigned to each sequenced from Biosampled submission - SUB[number]-BioSample-processed-ok-objects.csv
# library_ID				Unique identifier
# title						Short description
# library_strategy			Sequencing strategy - Random sequencing of the whole genome - WGS
# library-source			DNA type used - Genome DNA - GENOMIC
# library-selection			DNA selection strategy - Random selection by shearing or other method - RANDOM
# library-layout			Paired-end or single
# platform					Sequencing platform used - ILLUMINA
# instrument_model			Sequencing platform model used - Illumina MiSeq (NZ and Polyomics) and Illumina HiSeq 2500 (Sanger)
# design_description		Description of methods used to create the sequencing library
# filetype					Type of file being submitted - fastq
# filename					Name of file submitted
#
# Notes:
# PAIRED files must be put in the same line
# Save the file as a TAB delimited file - TSV

#######################################
# Get the BioSample Accession Numbers #
#######################################

# Open the file
my $biosampleFile = shift @ARGV;
open ACCESSION, $biosampleFile or die "Couldn't open file:$!";

# Initialise the necessary variables to parse the file
my $line;
my @cols;

# Initialise a hashtable to store the accession numbers link to the isolate IDs
my %isolateAccession = ();

# Begin reading the file
while(<ACCESSION>){
	$line = $_;
	chomp($line);
	
	# Skip the header lines
	next if $line =~ /^target_db/ || $line =~ /^BioSample/ || $line =~ /Successfully loaded/;
	
	# Split the current line into its columns
	@cols = split /,/, $line;
	
	# Store the accession number for the current isolate
	$isolateAccession{$cols[5]} = $cols[4];
}
close(ACCESSION);

###############################
# Get the isolates File Names #
###############################

# Get the directory containing the fastq files
my @directories = split /,/, shift @ARGV;

# Initialise a hashtable to store eahc isolates fastq file name
my %isolateFastqFiles = ();

# Initialise some necessary variables
my $dir;
my @files;
my $isolateID;
my @parts;

# Get the fastq files from each directory
foreach my $directory (@directories){

	# Get all the files in the directory
	opendir $dir, $directory or die "Cannot open directory: $!";
	@files = readdir $dir;

	# Add the fastq files present to the array
	foreach my $file (@files){
		
		if($file =~ /fastq.gz$/){
			
			# Get the isolate ID from the file name
			# 12754_7#2_1.fastq.gz
			# N11_S8_L001_R1_001.fastq.gz
			# AgR111_S8_L001_R1_001.fastq.gz
			@parts = split /_/, $file;
			$isolateID = $parts[0];
			if($file =~ /#/){
				$isolateID = $parts[0] . "_" . $parts[1];
				$isolateID =~ s/#/-/g;
			}
			
			# Store the file names (forward and reverse)
			if(exists($isolateFastqFiles{$isolateID})){
				$isolateFastqFiles{$isolateID} = $isolateFastqFiles{$isolateID} . "\t" . $file;
			}else{
				$isolateFastqFiles{$isolateID} = $file;
			}
		}
	}
}

######################
# Create Output File #
######################

# Open the output file
my $outputFile = shift @ARGV;
open OUTPUT, ">$outputFile" or die "Couldn't open file:$!";

# Print out the header
print OUTPUT "bioproject_accession\tbiosample_accession\tlibrary_ID\ttitle\tlibrary_strategy\tlibrary_source\tlibrary_selection\tlibrary_layout	platform\tinstrument_model\tdesign_description\tfiletype\tfilename\tfilename2\tfilename3\tfilename4\n";

# Get a list of the isolate ids
my @isolates = keys(%isolateAccession);

# Get the PRNJ number
my $bioprojectNumber = shift @ARGV;

# Initialise some necessary variables
my $outputLine;
my $platform;

# Print the necessary information for each isolate
foreach $isolateID (@isolates){
	
	# All have the same PRNJ - bioproject - number
	$outputLine = $bioprojectNumber;
	
	# Add the isolates biosample accession number
	$outputLine = $outputLine . "\t" . $isolateAccession{$isolateID};
	
	# Add isolate ID and other information
	$outputLine = $outputLine . "\t" . $isolateID . "\tNew Zealand Isolate\tWGS\tGENOMIC\tRANDOM\tPaired-end\tILLUMINA";
	
	# Note which platform the isolate was sequenced on
	$platform = "Illumina MiSeq";
	$platform = "Illumina HiSeq 2500" if $isolateID =~ /-/;
	$outputLine = $outputLine . "\t" . $platform;
	
	# Add further information and file type
	$outputLine = $outputLine . "\tLibrary preparation available on request\tfastq";
	
	# Add the file names
	$outputLine = $outputLine . "\t" . $isolateFastqFiles{$isolateID};
	
	print OUTPUT "$outputLine\n";
	
}

# Close the output file
close(OUTPUT);