#!/usr/bin/perl
# Searching files in a directory
use warnings;
use strict;
use Term::ANSIColor; # For Coloured Print Statements

# Author: Joseph Crisp
# Calculates the porportion of isolates that passed, or were rescued, at each site in the filtered.txt file

# Command Line Run Structure:
# perl CalculateIsolateCoverageFromFiltered.pl filtered.txt

# For each Variant Position within the filtered file the following information is available:
# #CHROM	POS	Sample 1:Sample 2:Sample 3:...      \t
# 0			1	2
#				|
#				 ---> 	DP	HQDP	MQ	AlleleSupport	AlleleCalled	QUAL	FQ		Result		RefAlt	;
# 						0 	1  		2 	3				4				5		6		7			8

# Fasta File Structure
# 	2 18					noSamples sequenceLength
#	>Sample1ID
#	GCTGATCGTGNGTAGGCC
#	>Sample2ID
#	GCTGATCGTGNGTAGGCT


##################
# Check For Help #
##################

my $firstArgument = shift @ARGV; # Open merged VCF File

# Print help information if needed
if($firstArgument eq "-help" || $firstArgument eq "" || $firstArgument eq "-h" || $firstArgument eq "help"){
	
	print color("blue"), "Perl Script to calculate the porportion of isolates that passed, or were rescued, at each site in the filtered.txt file\n";
	print "\nCommand Line Structure:\n";
	print "\tperl CalculateIsolateCoverageFromFiltered.pl filtered.txt\n", color("reset");
	print "\t\tfiltered.txt\t\tPath to the filtered VCF file.\n";
	
	exit 0;
	
}

#####################################
# Open filtered.txt and output file #
#####################################

# Open the filtered VCF file
my $filtered = $firstArgument;
open FILTERED, $filtered or die "Couldn't open $filtered:$!";

# Create an output file
my $date = substr((split /_/, $filtered)[2], 0, -4);
my $output = "IsolateVariantPositionCoverage_" . $date . ".txt";

# Print a check of the input information
print "\nInput Settings:\n";
print "Filtered VCF file: $filtered\n";
print "\nThe following output file is produced: $output\n\n";

####################################
# Calculate each isolates coverage #
####################################

# Keep track of what position in FASTA sequence investigating
my @isolates;
my @isolateCov;
my @isolateInfo;
my $nSites = 0;

# Initialise variables for reading the file
my $line;
my @cols;

# Begin reading the filtered file
while(<FILTERED>){
	my $line = $_;
	chomp($line);
	
	# Skip the header lines
	next if $line =~ /^##/;
	
	# Split the current line into its columns
	@cols = split /\t/, $line; 
	
	# Extract Sample IDs from Header Line
	#	#CHROM	POS	Sample 1:Sample 2:Sample 3:...
	# 	0 		1	2
	if($line =~ /^#/){
		
		# Get the isolate names
		@isolates = split /\:/, $cols[2]; # Extract the Sample Names into Array
		
		# Initialise the coverage array to store each isolates coverage
		@isolateCov = (0) x scalar(@isolates);
		next;
	}	
	
	# Increment the site count
	$nSites++;
	
	# Get each isolates information
	@isolateInfo = split /\:/, $cols[2];
	
	# Check whether each isolate's allele passed, failed or was rescued at the current position
	for(my $pos = 0; $pos < scalar(@isolateInfo); $pos++){
	
		# Skip isolates that failed
		next if $isolateInfo[$pos] =~ /Fail/ || $isolateInfo[$pos] =~ /----------/;
		
		# Add to current isolates count of sites that passed/were rescued
		$isolateCov[$pos]++;
	}

	# Progress information
	print "Finished examining site: $nSites\n" if $nSites % 1000 == 0;
	
}
close(FILTERED);

############################################
# Print each isolates coverage out to file #
############################################

# Open an output file
open OUTPUT, ">$output" or die "Couldn't open $output:$!";

# Print a header into the file
print OUTPUT "Isolate\tCoverage\n";

# Print the coverage of each isolate
for(my $pos = 0; $pos < scalar(@isolates); $pos++){
	
	# Convert count to proportion
	$isolateCov[$pos] = $isolateCov[$pos] / $nSites;
	
	# Print the coverage out
	print OUTPUT "$isolates[$pos]\t$isolateCov[$pos]\n";
}
close(OUTPUT);