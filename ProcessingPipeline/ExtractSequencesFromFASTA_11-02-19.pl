#!/usr/bin/perl
# Searching files in a directory
use warnings;
use strict;
use Term::ANSIColor; # For Coloured Print Statements

# Author: Joseph Crisp
# Extract's sequence - by index - from FASTA

# Command Line Run Structure: Should be ran in Directory with VCF Files
# perl ExtractSequenceFromFASTA.pl index sequences.fasta

#############
# Functions #
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
	my $helpStatement = "Perl Script to extract single sequence from FASTA file.\n\n" .
						"Command Line Structure:\n" .
						"\tperl ExtractSequenceFromFASTA.pl index sequences.fasta\n" .
						"\t\tindex\t\t\tA comma separated list of indices of sequences to extract (0-indexed). No spaces.\n" .
						"\t\tsequences.fasta\t\tA FASTA formatted file containing nucleotide sequences\n\n" .
						"This script prints the selected sequence into a file named: extractedSequences_[DATE].fasta\n";
	
	# Check if the correct number of arguments was provided
	if(scalar @arguments != 2){

		print $helpStatement;
		exit 0;
	}

	# Check if help requested
	if($arguments[0] eq "-help" || $arguments[0] eq "" || $arguments[0] eq "-h" || $arguments[0] eq "help" || $arguments[0] eq ""){

		print $helpStatement;
		exit 0;
	}
}

sub getSequenceIndices{

	# Split the argument string into array
	my @parts = split /,/, $_[0];

	# Convert the strings to integers
	my @sequenceIndices = ();
	for(my $i = 0; $i < scalar(@parts); $i++){

		$sequenceIndices[$i] = scalar($parts[$i]);
	}

	# Convert the array of integers to a hash
	my %indices = ();
	foreach my $index (@sequenceIndices){
		$indices{$index} = 1;
	}

	return %indices;
}

sub extractSequences{

	# Get the FASTA file name
	my $fastaFile = $_[0];

	# Get the indices of the sequences to be extracted
	my %indices = %{$_[1]};
	my $nSequencesToExtract = keys %indices;
	my $nSequencesExtracted = 0;

	# Open the Reference Sequence file
	open FASTA, $fastaFile or die "Couldn't open $fastaFile:$!";

	# Create an output file
	my $date = getCurrentDate();
	my $outputFile = "extractedSequences_" . $date . ".fasta";
	open OUTPUT, ">$outputFile" or die "Couldn't open $outputFile:$!";

	# Initialise a variable to note the number of sequences examined
	my $sequenceIndex = -1;

	# Initialise a variable to store the sequence and its ID
	my $name = "";
	my $sequence = "";

	# Initialise the necessary variables to parse the fasta file
	my $line;
	my $lineNo = 0;

	# Begin reading the Reference sequence file
	while(<FASTA>){
	
		# Get the current line and remove the end of line character
		$line = $_;
		chomp($line);

		# Increment the line number
		$lineNo++;

		# Check for FASTA info in first line - skip it if so
		next if $lineNo == 1 && $line !~ /^>/;
	
		# Check if found new sequence
		if($line =~ /^>/){

			# Either print or discard the previous sequence
			if(exists($indices{$sequenceIndex})){
				
				# Print the sequence to file
				print OUTPUT ">$name\n";
				print OUTPUT "$sequence\n";
				print "Extracted sequence for: $name\n";

				# Increment extracted sequence counter
				$nSequencesExtracted++;

				# Check if extracted all sequences requested
				if($nSequencesExtracted == $nSequencesToExtract){
					last;
				}
			}

			# Increment the sequence counter
			$sequenceIndex++;

			# Get the name of the current sequence
			$name = substr($line, 1); # offset by 1

			# Reset the sequence
			$sequence = "";
		}else{
		
			# Store the line as the next part of the sequence
			$sequence = $sequence . $line . "\n";
		}
	}

	# Check the last sequence
	if($nSequencesExtracted != $nSequencesToExtract && exists($indices{$sequenceIndex})){
		print OUTPUT ">$name\n";
		print OUTPUT "$sequence\n";
		print "Extracted sequence for: $name\n";
	}

	# Close the input and output files
	close(FASTA);
	close(OUTPUT);

	# Note the output file created
	print "Finished extracting sequences. Created output file: $outputFile\n";
}

#############
# Arguments #
#############

# Check the input arguments
checkInputArguments(\@ARGV);

# Note the index of the sequence to be extracted
my %indices = getSequenceIndices($ARGV[0]);

# Note the input FASTA file
my $fastaFile = $ARGV[1];

#####################
# Reading sequences #
#####################

# Extract and print selected sequences from FASTA
extractSequences($fastaFile, \%indices);