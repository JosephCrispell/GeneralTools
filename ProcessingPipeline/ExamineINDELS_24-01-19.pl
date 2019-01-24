#!/usr/bin/perl
# Searching files in a directory
use warnings;
use strict;
use Term::ANSIColor; # For Coloured Print Statements

# Author: Joseph Crispell
# Examining merged vcf file and creating a INDEL summary

# Command Line Run Structure:
# perl ExamineINDELS.pl filteredMerged.txt

# Merged VCF file structure
# For each Variant Position within the filtered file the following information is available:
# #CHROM	POS	Sample 1:Sample 2:Sample 3:...      \t
# 0			1	2
#				|
#				 ---> 	DP	HQDP	MQ	AlleleSupport	AlleleCalled	QUAL	FQ		Result		Ref	Alt	;
# 						0 	1  		2 	3				4				5		6		7			8	9

#############
# FUNCTIONS #
#############

sub getCurrentDate{

	# Get the information about the current date
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();

	# Build a date string
	my $date = $mday . "-" . ($mon + 1) . "-" . (1900 + $year);

	return $date;
}

sub checkInputArguments{
	
	# Get the array of arguments - first argument
	my @arguments = @{$_[0]};

	# Build the help statement
	my $helpStatement = "Perl Script to rescue variant position site information\n\n" .
						"Command Line Structure:\n" .
						"\tperl ExamineINDELS.pl filteredMerged.txt presenceAbsenceMatrix.txt indelIndex.txt\n" .
						"\t\tfilteredMerged.txt\t\tFiltered merged vcf file returned from FilterVariants.pl or RescureVariantPositionInfo.pl\n\n" .
						"This script generates the following output file: \n" .
						"\t\tpresenceAbsenceMatrix.txt\tPresence [1 or F (if failed)] absence [0] matrix where each sample represents a column and each INDEL a row.\n";
	
	# Check if the correct number of arguments was provided
	if(scalar @arguments != 1){

		print $helpStatement;
		exit 0;
	}

	# Check if help requested
	if($arguments[0] eq "-help" || $arguments[0] eq "" || $arguments[0] eq "-h" || $arguments[0] eq "help" || $arguments[0] eq ""){

		print $helpStatement;
		exit 0;
	}
}

sub indexSamples{

	# Get the array of sample IDs from the first argument
	my @sampleIDs = @{$_[0]};

	# Initialise a hashtable to store the IDs
	my %indexedSamples = ();

	# Examine each sample ID
	for(my $i = 0; $i < scalar(@sampleIDs); $i++){

		# Check already encountered ID
		if(exists($indexedSamples{$sampleIDs[$i]})){
			
			print color("red"), "ERROR: Sample ID: $sampleIDs[$i] is none unique!\n", color("reset");
			exit 0;
		}else{
			$indexedSamples{$sampleIDs[$i]} = $i;
		}
	}

	return %indexedSamples;
}

sub getINDEL{

	# Get the sample information - first argument
	my $info = $_[0];

	# Split the information into the quality scores
	my @parts = split /;/, $info;

	# Initialise a variable to store INDEL
	my $indel = "";

	# Check if INDEL present - either ref or alt are more than 1 character
	if(length($parts[8]) > 1 or length($parts[9]) > 1){
		$indel = $parts[8] . ":" . $parts[9];
	}

	return $indel;
}

sub recordINDELsInSamples{

	# Get the name of the merged file from the arguments
	my $filteredFile = $_[0];

	# Open the filtered file
	open FILTERED, $filteredFile or die "Couldn't open $filteredFile:$!";

	# Initialise variables to parse the filtered.txt file
	my $line;
	my @parts;
	my @sampleInfo;
	my $position;
	my $indel;

	# Initialise a hashtable to store the column index of each sample
	my %sampleIndices;
	my @sampleIDs;

	# Initialise a hashtable to store the samples that have each INDEL
	my %indelSamples = ();

	# Read the filtered file line by line
	while(<FILTERED>){
		$line = $_;
		
		# Remove the end of line character
		chomp($line);
		
		# Skip the header region
		if($line =~ /^##/){
			next;
		}

		# Get the information available for each sample
		@parts = split /\t/, $line;
		@sampleInfo = split /:/, $parts[2];


		# Get the sample IDs from the fields line
		if($line =~ /^#/){
			
			# Store and index the sample IDs in the sample info array
			@sampleIDs = @sampleInfo;
			%sampleIndices = indexSamples(\@sampleInfo);
			next;
		}

		# Note the current position on the genome
		$position = $parts[1];

		# Examine the quality information associated with each sample
		for(my $i = 0; $i < scalar(@sampleInfo); $i++){
		
			# Skip isolates that have no data available
			next if $sampleInfo[$i] eq "-------------------";
		
			# Get INDEL information
			$indel = getINDEL($sampleInfo[$i]);

			# Skip if INDEL found
			next if $indel eq "";

			# Add current position to INDEL found
			$indel = $position . ":" . $indel;

			# Check if already encountered INDEL
			if(exists($indelSamples{$indel})){

				# Add the current sample's ID to the array of IDs associated with the current INDEL - check if failed
				if ($sampleInfo[$i] =~ /Fail/){
					push (@{$indelSamples{$indel}}, $sampleIDs[$i] . ":FAIL");
				}else{
					push (@{$indelSamples{$indel}}, $sampleIDs[$i] . ":PASS");
				}
				
			}else{

				# Start an array recording the samples associated with the current INDEL - check if failed
				if ($sampleInfo[$i] =~ /Fail/){
					@{$indelSamples{$indel}} = ($sampleIDs[$i] . ":FAIL");
				}else{
					@{$indelSamples{$indel}} = ($sampleIDs[$i] . ":PASS");
				}
				
			}
		}
	}

	# Close the filtered file
	close(FILTERED);

	# Build the presence absence matrix file
	printPresenceAbsenceMatrix(\%indelSamples, \%sampleIndices, \@sampleIDs);
}

sub printPresenceAbsenceMatrix{

	# Get the input arguments
	my %indelSamples = %{$_[0]};
	my %sampleIndices = %{$_[1]};
	my @sampleIDs = @{$_[2]};

	# Create the output files to store the presence absence matrix and INDEL index
	my $date = getCurrentDate();
	my $presenceAbsenceFile = "indelInfo_presenceAbsence_" . $date . ".txt";

	# Open the output files
	open MATRIX, ">$presenceAbsenceFile" or die "Couldn't open $presenceAbsenceFile:$!";

	# Print headers into each output file
	print MATRIX "INDEL(Pos:Ref:Alt)\t" . join("\t", @sampleIDs) . "\n";

	# Get an array of each of the indels
	my @indels = keys %indelSamples;

	# Initialise a variable to note each sample's ID and filter result
	my @parts;

	# Initialise an array to act as a row in the presence absence matrix
	my @presenceAbsence;

	# Examine each of the INDELs
	for(my $i = 0; $i < scalar(@indels); $i++){

		# Create an array of zeros - length equals number of samples
		@presenceAbsence = (0) x scalar(@sampleIDs);

		# Examine each of the samples with the current INDEL
		foreach my $sample (@{$indelSamples{$indels[$i]}}){

			# Get the sample ID - without filter result
			@parts = split /:/, $sample;

			# Note that the current sample has the current INDEL - check if failed
			if($parts[1] eq "FAIL"){
				$presenceAbsence[$sampleIndices{$parts[0]}] = "F";
			}else{
				$presenceAbsence[$sampleIndices{$parts[0]}] = 1;
			}
			
		}

		# Print the current row in the presence absence matrix (for the current INDEL)
		print MATRIX $indels[$i] . "\t" . join("\t", @presenceAbsence) . "\n";
	}

	# Close the output files
	close(MATRIX);
}

################################
# Examining INDELs in filtered #
################################

# Check the input arguments
checkInputArguments(\@ARGV);

# Get the input argument
my $filtered = $ARGV[0];

# Record which INDELs are present in each isolate
recordINDELsInSamples($filtered);