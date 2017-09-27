#!/usr/bin/perl
# Searching files in a directory
use warnings;
use strict;

# Author: Joseph Crisp
# Convert the Sequence FASTA File to a NEXUS File format and Create an accompanying TRAITS File

# Command Line Structure
# perl ConvertFasta2Phylip.pl sequences.fasta sequences.phy 

# FASTA File Format:
#	4 20				noSamples	noSites
#	>AgR111_S8_1.vcf
#	GGTCGCGCTCCGACACGCCC
#	>AgR121_S1_2.vcf
#	GGTCGCGCTCCGNCANGCCC
#	>12755_8#21_202.vcf
#	GGTCGCGCTCCGACANGCCC
#	>12754_7#93_92.vcf
#	GGTCGCGCTCCGACGNGCCC

# PHYLIP File Format:
# 		3	20
#	isolateID	ATGCTCGATAGCTAGATCGA
#	isolateID	ATGCTCGATAGCTAGATCGA
#	isolateID	ATGCTCGATAGCTAGATCGA
#	
# Note that the isolateID can only be 10 characters

# Open in the Sample FASTA File
my $fastaFile = shift @ARGV; 
open FASTA, $fastaFile or die "Couldn't open $fastaFile:$!";

# Open the Output File
my $phylipFile = shift @ARGV; # Takes Output File name from Command Line Input
open PHYLIP, ">$phylipFile" or die "Couldn't open $phylipFile:$!";

# Create variables to store the isolateIDs and sequence
my $sequenceName;
my $sequence;
my $outputLine;
my $isolateID;

# Create the necessary variables for parsing the FASTA file
my $line;
my @cols;
my $lineNo = 0;
my @parts;

# Begin reading the FASTA file
while(<FASTA>){
	$lineNo++;
	$line = $_;
	chomp($line);
	
	# Get the number of isolates and sequence length from the first line
	if($lineNo == 1){
		@cols = split /\ +/, $line;
		
		print PHYLIP "\t$cols[0]\t$cols[1]\n";
		next;
	}
	
	# Have we found a new isolate?
	if($line =~ /^>/){
		
		# Print out the information for the previous isolate
		if($lineNo != 2){
			
			# Get the isolates ID
			@parts = split /\_/, $sequenceName;
			$isolateID = $parts[0];
			
			# Pad the isolateId with spaces to match 10 character rule
			$outputLine = sprintf "%-10s %-s\n", $isolateID, $sequence;
			
			print PHYLIP "$outputLine";			
		}
		
		# Get the isolate name
		$sequenceName = substr($line, 1);
		
		# Reset the sequence
		$sequence = "";
	}else{
		$sequence = $sequence . "" . $line;
	}
}

# Print out the information for the last isolate
# Get the isolates ID
@parts = split /\_/, $sequenceName;
$isolateID = $parts[0];
	
# Pad the isolateId with spaces to match 10 character rule
$outputLine = sprintf "%-10s %-s\n", $isolateID, $sequence;
print PHYLIP "$outputLine";

# Close the input and output files
close(FASTA);
close(PHYLIP);

