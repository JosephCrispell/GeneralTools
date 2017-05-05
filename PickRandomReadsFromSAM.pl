#!/usr/bin/perl
# Searching files in a directory
use warnings;
use strict;

# Import package for shuffling array
use List::Util qw(shuffle);

# Author: Joseph Crisp
# Pick a random number of reads from a SAM file of unmapped reads

# Command Line Structure:
# perl PickRandomReadsFromSAM.pl n unmapped.sam selectedReads.txt

# SAM file structure:
# M01808:112:000000000-AHA0P:1:1101:11099:9551	77	*	0	0	*	*	0	0	CTGACATCGG	GGGGGGGGGG	AS:i:0	XS:i:0
# 0												1	2	3	4	5	6	7	8	9			10			11		12

# Search BLAST databases for query sequences command:
# sudo blastn -query query.fasta -out blast.txt -db nr -remote


###########################
# Command Line Arguements #
###########################

# Get the number of reads to select from the SAM file
my $nReads = shift @ARGV;

# Open the SAM file containing the unmapped reads
my $samFile = shift @ARGV;
open SAM, $samFile or die "Couldn't open $samFile:$!";

################################
# Choose the Random file lines #
################################

# Find the number of lines in the file
my $nLines = (split /\ /, `wc -l $samFile`)[0];

# Pick n random lines to extract the reads from
my @shuffledIndices = shuffle(0 .. $nLines - 1);
my @selectedLines = @shuffledIndices[0 .. $nReads - 1];

# Store the picked indices in a hashtable
my %pickedLines = ();
for(my $i = 0; $i < scalar(@selectedLines); $i++){
	$pickedLines{$selectedLines[$i]} = 1;
}

#######################################
# Extract the Reads from the SAM file #
#######################################

# Initialise the necessary variables for parsing the SAM file
my $line;
my @cols;
my $lineNo = -1;
my $nFound = 0;

# Begin reading the SAM file
while(<SAM>){
	$line = $_;
	chomp($line);
	$lineNo++;
	
	# Was the current line selected?
	if(exists($pickedLines{$lineNo})){
		
		# Store the read sequence
		@cols = split /\t/, $line;
		$pickedLines{$lineNo} = $cols[9];
		$nFound++;
	}
	
	# Break out once we have found all the lines selected
	if($nFound == $nReads){
		last;
	}
}
close(SAM);

########################################
# Print the selected reads out to file #
########################################

my $outputFile = shift @ARGV;
open OUTPUT, ">$outputFile" or die "Couldn't open $outputFile:$!";

my @keys = keys %pickedLines;

for(my $i = 0; $i < scalar(@keys); $i++){
	print OUTPUT ">Query_$i\n$pickedLines{$keys[$i]}\n";
}

close(OUTPUT);
