#!/usr/bin/perl
use warnings;
use strict;

# Author: Joseph Crisp
# Get part of FASTA sequence

# Command Line Structure
# perl GetPartOfFASTASequence.pl start end sequence.fasta label
#
# Note only reads first sequence!!

###############################
# Get command line arguements #
###############################

# Get start and end of region of interest
my $start = shift @ARGV;
my $end = shift @ARGV;

# Open the input fasta file
my $fastaFile = shift @ARGV;
open FASTA, $fastaFile or die "Couldn't open $fastaFile:$!";

# Get region name
my $name = shift @ARGV;

####################################
# Read in sequence from FASTA file #
####################################

# Initialise a variable to store the sequence
my $sequence = "";

# Initialise the necessary variables for parsing the file
my $line;
my $lineNo = 0;

# Begin reading the file
while(<FASTA>){
	$lineNo++;
	$line = $_;
	chomp($line);	

	# Skip the first line
	if($line =~ /^>/){
		
		next if $lineNo == 1;
		last if $lineNo > 1;
	}
	
	# Store the sequence
	$sequence = $sequence . "" . $line;
}

# Close the input FASTA file
close(FASTA);

##################
# Get the region #
##################

# Get the region
my $region = substr($sequence, $start - 1, $end - $start);

# Convert to array of nucleotides
my @nucleotides = split //, $region;

# Get the reverse sequence
my @reverse = reverse @nucleotides;

# Get the compliments of the forward and reverse
my @compliment = GetCompliment(\@nucleotides);
my @reverseCompliment = GetCompliment(\@reverse);

# Convert the nucleotide sequences to strings
my $region_reverse = join("", @reverse);
my $region_compliment = join("", @compliment);
my $region_reverseCompliment = join("", @reverseCompliment);

####################################
# Print the region into FASTA file #
####################################

print ">$name\n$region\n";
print ">$name-compliment\n$region_compliment\n";
print ">$name-reverse\n$region_reverse\n";
print ">$name-reverseCompliment\n$region_reverseCompliment\n";


#############
# FUNCTIONS #
#############

sub GetCompliment{
   
	# Get the sequence as the input arguement
	my @nucleotides = @{$_[0]};
   
	# Define hashtable of nucleotide compliments
	my %compliments = (
		A => 'T',
		C => 'G',
		G => 'C',
		T => 'A'
	);
	
	# Initialise an array to store the compliment
	my @compliment = ();
	
	# Build the compliment
	for(my $i = 0; $i < scalar(@nucleotides); $i++){
		$compliment[$i] = $compliments{$nucleotides[$i]};
	}
	
	return @compliment;
}