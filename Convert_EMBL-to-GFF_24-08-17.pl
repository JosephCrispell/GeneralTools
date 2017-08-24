#!/usr/local/bin/perl -w
use strict;
use Bio::SeqIO;

# Script to convert EMBL format to GFF format
# Command line structure:
# perl Convert_EMBL-to-GFF.pl file.embl > output.gff

## Script edited from: http://ratt.sourceforge.net/transform.html

# Check that input arguements were given correctly
if (scalar(@ARGV) != 1) {
    die "USAGE: Convert_EMBL-to-GFF.pl file.embl > output.gff\n"; 
}

###########################################
# Convert input file into gff output file #
###########################################

# Note that using bioperl functions

# Initialise variables for creating the output file
my $line;
my @cols;

# Print header

# Read in the input file
my $embl = shift @ARGV;
my $input = Bio::SeqIO->new(-file=>$embl,-format=>'EMBL');
while (my $seq = $input->next_seq) {
  for my $feature ($seq->top_SeqFeatures) {
  
	# Get the new line
	$line = $feature->gff_string;
	
	# Split the line into its columns
	@cols = split /\t/, $line;
	
	# Replace the second column with "feature"
	$cols[1] = "feature";
	
	# Remake the line
	$line = join("\t", @cols);
	
	print "$line\n";
  }
}