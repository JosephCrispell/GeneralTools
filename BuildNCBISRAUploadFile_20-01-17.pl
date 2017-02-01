#!/usr/bin/perl
# Searching files in a directory
use warnings;
use strict;

# Author: Joseph Crisp
# Extract the sampling information needed to accompany the FASTQ files on NCBI

# Command Line Structure
# perl BuildNCBISRAUploadFile_20-01-17.pl BioSampleInfo.csv FASTQDirectory1,FastQDirectory2 outputFile.tsv PRJN

# BioSampleInfo.csv file structure:
# 	target_db	status			msg	resp_status	accession		name			tax_id	extras	file_path	file_id	file_url
# 	0			1				2				3				4				5		6		7			8		9
# 	BioSample					processed-ok
#								Successfully loaded
#												processed-ok	SAMN06245910	AgR297	1765

# Output file structure:
# 	bioproject_accession	biosample_accession	library_ID	title/short description	library_strategy
#	0 ###					1 ###				2 ###		3 ###					4 ###
#			
#	library_source	library_selection	library_layout	platform	instrument_model	design_description
#	5 ###			6 ###				7 ###			8 ###		9 ###				10 ###
#
#	reference_genome_assembly (or accession)	alignment_software	forward_read_length	reverse_read_length	filetype
#	11											12					13 ###				14 ###				15 ###
#
#	filename	MD5_checksum	filetype	filename	MD5_checksum
#	16 ###		17 ###			18 ###		19 ###		20 ###

# Column Definitions:
# bioproject_accession		PRJNA[number] - the number assigned to the BioProject submission
# biosampled_accession		SAMN[number] - the numbers assigned to each sequenced from Biosampled submission - SUB[number]-BioSample-processed-ok-objects.csv
# library_ID				Unique identifier
# title/short description	Short description
# library_strategy			Sequencing strategy - Random sequencing of the whole genome - WGS
# library-source			DNA type used - Genome DNA - GENOMIC
# library-selection			DNA selection strategy - Random selection by shearing or other method - RANDOM
# library-layout			paired or single
# platform					Sequencing platform used - ILLUMINA
# instrument_model			Sequencing platform model used - Illumina MiSeq (NZ and Polyomics) and Illumina HiSeq 2500 (Sanger)
# design_description		Description of methods used to create the sequencing library
# reference_genome_assembly For aligned data (accession number)
# alignment_software		For aligned data
# forward_read_length		Length of the Forward Reads
# reverse_read_length		Length of the Reverse Reads
# filetype					Type of file being submitted - fastq	FORWARD
# filename					Name of file submitted					FORWARD
# MD5_checksum				MD5 check sum value for file			FORWARD
# filetype					Type of file being submitted - fastq	REVERSE
# filename					Name of file submitted					REVERSE
# MD5_checksum				MD5 check sum value for file			REVERSE

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

# Initialise a hashtable to store each isolates fastq file name
my %isolateFastqFiles = ();
my @fastqFilesInfo;
my $md5Checksum = 15;

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
			}
			
			# Get the md5 checksum value for the current file
			$md5Checksum = `md5sum $directory"/"$file`;
			$md5Checksum = (split / +/, $md5Checksum)[0];
			print "md5 checksum for $file: $md5Checksum\n";
			
			# Store the file names (forward and reverse)
			if(exists($isolateFastqFiles{$isolateID})){
				@fastqFilesInfo = @{$isolateFastqFiles{$isolateID}};
				$fastqFilesInfo[2] = $file;
				$fastqFilesInfo[3] = $md5Checksum;
				$fastqFilesInfo[4] = $directory;
				@{$isolateFastqFiles{$isolateID}} = @fastqFilesInfo;
			}else{
				@fastqFilesInfo = {};
				$fastqFilesInfo[0] = $file;
				$fastqFilesInfo[1] = $md5Checksum;
				@{$isolateFastqFiles{$isolateID}} = @fastqFilesInfo;
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
print OUTPUT "bioproject_accession\tbiosample_accession\tlibrary_ID\ttitle/short description\tlibrary_strategy\tlibrary_source\tlibrary_selection\tlibrary_layout\tplatform\tinstrument_model\tdesign_description\treference_genome_assembly\talignment_software\tforward_read_length\treverse_read_length\tfiletype\tfilename\tMD5_checksum\tfiletype2\tfilename2\tMD5_checksum2\n";

# Get a list of the isolate ids
my @isolates = keys(%isolateAccession);

# Get the PRNJ number
my $bioprojectNumber = shift @ARGV;

# Initialise some necessary variables
my $outputLine;
my $platform;
my $readLength;
my $error;

# Print the necessary information for each isolate
foreach $isolateID (@isolates){
	
	# All have the same PRNJ - bioproject - number
	$outputLine = $bioprojectNumber;
	
	# Add the isolates biosample accession number
	$outputLine = $outputLine . "\t" . $isolateAccession{$isolateID};
	
	# Add isolate ID and other information
	$outputLine = $outputLine . "\t" . $isolateID . "\tNew Zealand Isolate\tWGS\tGENOMIC\tRANDOM\tpaired\tILLUMINA";
	
	# Note which platform the isolate was sequenced on
	$platform = "Illumina MiSeq";
	$platform = "Illumina HiSeq 2500" if $isolateID =~ /#/;
	$outputLine = $outputLine . "\t" . $platform;
	
	# Add further information and file type
	$outputLine = $outputLine . "\tLibrary preparation available on request\t\t";
	
	# Add the length of the forward and reverse reads
	$readLength = 300;
	$readLength = 250 if $isolateID =~ /^N/;
	$readLength = 100 if $isolateID =~ /#/;
	$outputLine = $outputLine . "\t" . $readLength . "\t" . $readLength;
	
	# Add the forward and reverse file information
	if(exists($isolateFastqFiles{$isolateID})){
		@fastqFilesInfo = @{$isolateFastqFiles{$isolateID}};
		$outputLine = $outputLine . "\tfastq\t" . $fastqFilesInfo[0] . "\t" . $fastqFilesInfo[1] . "\tfastq\t" . $fastqFilesInfo[2] . "\t" . $fastqFilesInfo[3];
	
		print OUTPUT "$outputLine\n";

		# Copy the fastq files into an upload folder
		print "Copying files for $isolateID...\n";
		$error = `cp $fastqFilesInfo[4]"/"$fastqFilesInfo[0] "FilesToUpLoad/"`;
		$error = `cp $fastqFilesInfo[4]"/"$fastqFilesInfo[2] "FilesToUpLoad/"`;
	}else{
		print "ERROR!: Couldn't find the file for isolate: $isolateID\n";
	}
	
	
	
	
}

# Close the output file
close(OUTPUT);