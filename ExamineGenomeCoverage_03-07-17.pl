#!/usr/bin/perl
# Searching files in a directory
use warnings;
use strict;
use Term::ANSIColor; # For Coloured Print Statements

# Author: Joseph Crisp
# Investigate the coverage of the sequenced isolates on the reference genome

# Command Line Structure:
# perl ExamineGenomeCoverage.pl threshold genomeCoverage.txt genomeCoverageSummary.txt missingPositions.txt isolateCoverageSummary.txt 

# Genome coverage file structure:
# POS	filenameA:filenameB:filenameC
# 1		20:2:5
# 		DPA:DPB:DPC
# 0		1

# Get first variable from command line
my $firstArg = shift @ARGV;

# Print help information if requested
if($firstArg eq "-help" || $firstArg eq "" || $firstArg eq "-h"){
	print color("green"), "Perl Script to remove isolates with poor site coverage from FASTA.\n\nCommand Line Structure:\n", color("reset");
	print "\tperl ExamineGenomeCoverage.pl threshold genomeCoverage.txt\n";
	print "\t\tthreshold\t\tThreshold for number mapped reads aligned for site to be considered as covered.\n";
	print "\t\tgenomeCoverage_date.txt\tPath to table noting the read depth of each isolate for each position on reference genome. File should have date as last delimited part of name\n";

	exit 0;	
}

#######################################################################
# Examine the isolate coverage at each position of the Reference Genome

# Get the Read Depth threshold from the command line
my $threshold = $firstArg;

# Open the IsolateCoverage file
my $coverage = shift @ARGV;
open COVERAGE, $coverage or die "Couldn't open $coverage:$!";

# Get the date from the file name
my $date = (split /\./, (split /_/, $coverage)[-1])[0];

# Open an output genome coverage file
my $coverageOutput = "genomeCoverageSummary_DP-$threshold". "_$date.txt";
open GENOMECOV, ">$coverageOutput" or die "Couldn't open $coverageOutput:$!";

# Open a file to record any gaps with no coveraeg in any isolates
my $missingPositions = "missingPositions_DP-$threshold". "_$date.txt";
open MISSING, ">$missingPositions" or die "Couldn't open $missingPositions:$!";

# Create a name for the output file that reads the coverage of each isolate
my $isolateCov = "isolateCoverageSummary_DP-$threshold". "_$date.txt";

# Check the input information
print "Input settings:\n";
print "\tRead Depth threshold = $threshold\n";
print "\tGenome coverage file: $coverage\n";
print "\nThe following output files will be produced:\n";
print "\t$coverageOutput\tTable recording the average read depth, across the isolates, for each position on the reference genome\n";
print "\t$missingPositions\t\tFile recording the postions on the reference genome that no isolates had sufficient aligned reads\n";
print "\t$isolateCov\tTable recording the average read depth, across the reference genome, of each isolate\n\n";

# Initialise a hashtable to store information for each isolate: Index: IsolateID
my @isolateInfoAtIndex = ();
my @info = ();
my @isolateNames = ();

# Initialise variables necessary to fill in positions where there was no coverage in any isolates
my $nPositionsMissed;
my $outputLine;

# Initialise two arrays for calculating the mean coverage and percentage coverage
my @meanDepth = ();
my @percentageCov = ();

# Initialise a variable to record the average coverage at each position on the genome
my $averageDepth;

# Initialise the necessary variables to parse the file
my $line;
my @cols;
my @parts = ();
my $lineNo = 0;
my $pos;

print color("green"), "Reading genome coverage file...\n", color("reset");

# Initialise an array to count the number of cattle and badgers with coverage above th threshold
my @nSpecies = ();

# Begin reading the genome coverage file
while(<COVERAGE>){
	$line = $_;
	chomp($line);
	@cols = split /\t/, $line;
	
	# Examine the header
	if($line =~ /^POS/){
		
		print GENOMECOV "POS\tMeanDepth\n";
		
		# Note the column index of each isolate and the species
		# **** VCF file name begins with TB for cattle isolate and WB for badger
		for(my $i = 1; $i < scalar(@cols); $i++){
		
			# Note the isolate id at the current column
			$isolateInfoAtIndex[$i - 1] = $cols[$i];
		}
		
		# Initialise the arrays for calculating mean depth and percentage coverage
		@meanDepth = (0) x (scalar(@cols) - 1);
		@percentageCov = (0) x (scalar(@cols) - 1);
		next;
	}	
	
	# Record the current line
	$lineNo++;
	
	# Have we skip any genome positions since the last line?
	if($lineNo != $cols[0]){
		
		# Calculate how many positions have been skipped
		$nPositionsMissed = $cols[0] - $lineNo;
		print MISSING "Missing Positions Found: $lineNo\t$cols[0]\t$nPositionsMissed\n";
		
		# Insert a line in the output table
		for(my $i = 0; $i < $nPositionsMissed; $i++){
			
			$pos = $lineNo + $i;
			print GENOMECOV "$pos\t0\t0\n";
		}
		
		# Update the line number
		$lineNo = $cols[0];
	}
	
	# Reset the average depth variable
	$averageDepth = 0;
	
	for(my $i = 1; $i < scalar(@cols); $i++){
		
		# Calculating the average read depth at the current position
		$averageDepth += $cols[$i];
		
		if($cols[$i] >= $threshold){
		
			# Make the calculations for the percentage coverage
			$percentageCov[$i - 1]++ if $cols[$i] >= $threshold;
		}
		
		# Make the calculations for the mean depth for the current isolate
		$meanDepth[$i - 1] += $cols[$i];
	}
	
	# Calculate the average read depth at the current position across the isolates
	$averageDepth = $averageDepth / scalar(@isolateInfoAtIndex);
	
	# Print out the summary information for the current genome position
	$outputLine = $cols[0] . "\t" . $averageDepth . "\t";
	
	print GENOMECOV "$outputLine\n";
	
	# Print progress information
	if($lineNo % 100000 == 0){
		print "Finished reading current line:\t$lineNo\n";
	}
}	
close(COVERAGE);
close(GENOMECOV);
close(MISSING);

# Print the mean depth and percentage coverage for each isolate
open ISOLATECOV, ">$isolateCov" or die "Couldn't open $isolateCov:$!";

# Print a Header
print ISOLATECOV "IsolateID\tMeanDepth\tPercentageCoverage\n";

# Initialise a variable to store the isolate ID
my $isolate;

for(my $i = 0; $i < scalar(@isolateInfoAtIndex); $i++){
	
	# Calculate the mean depth for the current isolate
	my $mean = $meanDepth[$i] / $lineNo;
	
	# Calculate the percentage coverage
	my $percentage = $percentageCov[$i] / $lineNo;
	
	# Get the isolate info
	$isolate = $isolateInfoAtIndex[$i];
	
	# Print out the isolate's mean depth and percentage coverage
	print ISOLATECOV "$isolate\t$mean\t$percentage\n";
}
close(ISOLATECOV);

print color("blue"), "Finished.\n", color("reset");