#!/usr/bin/perl
# Searching files in a directory
use warnings;
use strict;

# Author: Joseph Crisp
# Convert paired FASTQ files into FASTA (VERY CRUDE!)

# Command Line Structure
# perl ConvertPairedFastqsIntoFasta.pl forward.fastq reverse.fastq

###############################################
# Open the Forward FASTQ file and store Reads #
###############################################

# Forward
my $forward = shift @ARGV;
open FORWARD, $forward or die "Couldn't open $forward:$!";

# Initialise variables to parse file
my $line;
my $sequenceInNextLine = 0;

# Initialise an array to store the sequence of each read
my @forwardReads = ();
my $pos = -1;

# Begin reading file and storing reads
print "Reading sequences from Forward FASTQ file...\n";
while(<FORWARD>){
	$line = $_;
	chomp($line);

	# Check if found new read tag
	if($line =~ /^@/){
		$sequenceInNextLine = 1;
		
	}elsif($sequenceInNextLine == 1){
		
		$pos++;
		$forwardReads[$pos] = $line;
		$sequenceInNextLine = 0;
	}
}
close(FORWARD);

###############################################
# Open the Reverse FASTQ file and store Reads #
###############################################

# Reverse
my $reverse = shift @ARGV;
open REVERSE, $reverse or die "Couldn't open $reverse:$!";

# Reset variables to parse file
$sequenceInNextLine = 0;

# Initialise an array to store the sequence of each read
my @reverseReads = ();
$pos = -1;

# Begin reading file and stroing reads
print "Reading sequences from Reverse FASTQ file...\n";
while(<REVERSE>){
	$line = $_;
	chomp($line);

	# Check if found new read tag
	if($line =~ /^@/){
		$sequenceInNextLine = 1;
		
	}elsif($sequenceInNextLine == 1){
		
		$pos++;
		$reverseReads[$pos] = $line;
		$sequenceInNextLine = 0;
	}
}
close(REVERSE);

##################################################################
# Print out the Forward and Reverse Reads into FASTA file format #
##################################################################

# Name and open output FASTA file
my @fileNameParts = split /_/, $forward;
my $fasta = $fileNameParts[0] . "_" . $fileNameParts[1] . "_RawReads.fasta";
open OUTPUT, ">$fasta" or die "Coudln't open $fasta:$!";

# Print forward reads
print "Printing forward reads to FASTA file...\n";
for($pos = 0; $pos < scalar(@forwardReads); $pos++){
	print OUTPUT ">F-$pos\n$forwardReads[$pos]\n" if $forwardReads[$pos] !~ /^(.)\1{1,}$/;
}

# Print forward reads
print "Printing reverse reads to FASTA file...\n";
for($pos = 0; $pos < scalar(@reverseReads); $pos++){
	print OUTPUT ">R-$pos\n$reverseReads[$pos]\n" if $reverseReads[$pos] !~ /^(.)\1{1,}$/;
}

# Close the FASTA file
close(OUTPUT);