#!/usr/bin/perl
# Searching files in a directory
use warnings;
use strict;
use Term::ANSIColor; # For Coloured Print Statements

# Author: Joseph Crisp
# For sites, of isolates, that fail filtering - this script looks back at data available to rescue an allele where possible

# Command Line Run Structure:
# perl RescueVariantPositionInfo.pl nIsolatesAlleleMustBePresentIn highQualityBaseDepthThreshold alleleSupportThreshold filtered.txt

# For each Variant Position within the filtered file the following information is available:
# #CHROM	POS	Sample 1:Sample 2:Sample 3:...      \t
# 0			1	2
#				|
#				 ---> 	DP	HQDP	MQ	AlleleSupport	AlleleCalled	QUAL	FQ		Result		Ref	Alt	;
# 						0 	1  		2 	3				4				5		6		7			8	9

#############
# FUNCTIONS #
#############

sub countAllelesAtCurrentSite{
	
	# Get the isolate info from the input
	my @isolateInfo = @{$_[0]};
	
	# Initialise a hashtable to store the alelle counts
	my %alleles = ();
	
	# Initialise a variable to store each isolates quality
	my @qualityInfo;
	
	# Count the alleles present, and how many isolates they were found in at the current site
	for(my $i = 0; $i < scalar(@isolateInfo); $i++){
		
		# Skip isolates with no information ior that failed filtering
		next if $isolateInfo[$i] =~ /^-----/ || $isolateInfo[$i] =~ /Fail/;
		
		# Split the information available for the current isolate at the current site
		@qualityInfo = split /;/, $isolateInfo[$i];
		
		# Check which allele was called
		if(exists($alleles{$qualityInfo[4]})){
			$alleles{$qualityInfo[4]}++;
		}else{
			$alleles{$qualityInfo[4]} = 1;
		}
	}
	
	return %alleles;
}

sub returnSupportedAllele{
	
	# Get the allele depths from the input
	my @alleleHighQualityBaseDepth = @{$_[0]};
	
	# Initialise an index variable
	my $index = 0;
	
	# Check if alternate allele has more high quality bases supporting it
	if($alleleHighQualityBaseDepth[1] > $alleleHighQualityBaseDepth[0]){
		$index = 1;
	}
	
	return $index;
}

sub checkForAlleleToRescue{
	
	# Get the input arguements
	my $isolateInfo = $_[0];
	my $nIsolatesAlleleMustBePresentIn = $_[1];
	my $depthThreshold = $_[2];
	my $supportThreshold = $_[3];
	my %alleles = %{$_[4]};
	
	# Get the quality information for the current isolate at the current position
	my @qualityInfo = split /;/, $isolateInfo;
		
	# Note alleles present
	my @allelesPresent = ($qualityInfo[8], $qualityInfo[9]);
	my @highQualityBaseDepth = split /,/, $qualityInfo[1];
		
	# Count high quality base depth
	my @depthForAllelesPresent = ();
	$depthForAllelesPresent[0] = $highQualityBaseDepth[0] + $highQualityBaseDepth[1];
	$depthForAllelesPresent[1] = $highQualityBaseDepth[2] + $highQualityBaseDepth[3];
		
	# Note which allele is best supported by the high quality
	my $chosen = returnSupportedAllele(\@depthForAllelesPresent);
		
	# Calculate the proportion of high quality bases supporting allele
	my $proportionHighQualityBasesSupportingAllele = 0;
	if($depthForAllelesPresent[0] + $depthForAllelesPresent[1] > 0){
		$proportionHighQualityBasesSupportingAllele = $depthForAllelesPresent[$chosen] / ($depthForAllelesPresent[0] + $depthForAllelesPresent[1]);
	}
			
	# Check if chosen allele was observed in a sufficient number of isolates and is supported with a sufficient number of high quality bases
	if(exists($alleles{$allelesPresent[$chosen]}) && 
	$alleles{$allelesPresent[$chosen]} >= $nIsolatesAlleleMustBePresentIn && $depthForAllelesPresent[$chosen] >= $depthThreshold && $proportionHighQualityBasesSupportingAllele >= $supportThreshold){
		
		# Re-build the current isolates information at the current site
		$isolateInfo = $qualityInfo[0] . ";" .	$qualityInfo[1] . ";" .	$qualityInfo[2] . ";" . $proportionHighQualityBasesSupportingAllele;
		$isolateInfo .= ";" . $allelesPresent[$chosen] . ";" . 	$qualityInfo[5] . ";" .	$qualityInfo[6] . ";Rescue;" . $qualityInfo[8] . ";" . $qualityInfo[9];
	}

	return $isolateInfo;
}

################################
# Check command line arguments #
################################

# Get first argument
my $argument = shift @ARGV;

# Check if asking for help
if($argument eq "-help" || $argument eq "" || $argument eq "-h" || $argument eq "help"){
	print color("green"), "Perl Script to rescue variant position site information\n\nCommand Line Structure:\n", color("reset");
	print "\tperl RescueVariantPositionInfo.pl nIsolateAlleleMustBePresentIn highQualityBaseDepthThreshold alleleSupportThreshold filtered.txt\n";
	print "\t\tnIsolateAlleleMustBePresentIn\tN. isolates that allele must be found in to be considered for rescuing\n";
	print "\t\thighQualityBaseDepthThreshold\tN. high quality forward/reverse bases necessary to support a rescued allele\n";
	print "\t\talleleSupportThreshold\t\tProportion of high quality bases necessary to support a rescued allele\n";
	print "\t\tfiltered.txt\t\t\tFull path to output from FilterVariants.pl\n";

	exit 0;
}

##################################
# Get the command line arguments #
##################################

# Get thresholds
my $nIsolatesAlleleMustBePresentIn = $argument; # Number of isolates allele must have been observed in
my $depthThreshold = shift @ARGV; # Number of high quality reads supporting allele
my $supportThreshold = shift @ARGV; # Proportion of high quality reads supporting allele

# Reassign filtered.txt file argument
my $filteredFile = shift @ARGV;

# Build output file name
my $date = substr((split /_/, $filteredFile)[2], 0, -4);
my $output = "filtered-rescued_" . $nIsolatesAlleleMustBePresentIn . "-" . $depthThreshold . "-" . $supportThreshold . "_" . $date . ".txt";

# Print arguments back to user
print color("green"), "Input argument settings:\n", color("reset");
print "IsolateAlleleMustBePresentIn: \t$nIsolatesAlleleMustBePresentIn\n";
print "highQualityBaseDepthThreshold: \t$depthThreshold\n";
print "alleleSupportThreshold: \t$supportThreshold\n";
print "filtered.txt: \t\t\t$filteredFile\n";
print color("green"), "\nProduces the following output file: ", color("reset");
print "$output\n";

###############################################
# Examine each site in the filtered text file #
###############################################

# Open the file
open FILTERED, $filteredFile or die "Couldn't open $filteredFile:$!";

# Open an output file
open OUTPUT, ">$output" or die "Couldn't open $output:$!";

# Initialise variables to parse the filtered.txt file
my $line;
my $lineNo = 0;
my @parts;

# Initialise variables to handle each isolates quality information at each site
my @isolateInfo;

# Initialise variables to store each site summary
my %alleles;

# Read the filtered file line by line
while(<FILTERED>){
	$line = $_;
	chomp($line);
	$lineNo++;
	
	# Skip the header lines
	if($line =~ /^#/){
		print OUTPUT "$line\n";
		next;
	}
	
	# Get the information available for each isolate
	@parts = split /\t/, $line;
	@isolateInfo = split /:/, $parts[2];
	
	# Count the alleles present, and how many isolates they were found in at the current site
	%alleles = countAllelesAtCurrentSite(\@isolateInfo);
	
	# Examine the quality information available for each isolate
	for(my $i = 0; $i < scalar(@isolateInfo); $i++){
		
		# Skip isolates that didn't fail
		next unless $isolateInfo[$i] =~ /Fail/;
		
		# Check whether there is an allele to rescued for the current isolate at the current position
		$isolateInfo[$i] = checkForAlleleToRescue($isolateInfo[$i], $nIsolatesAlleleMustBePresentIn, $depthThreshold,
		$supportThreshold, \%alleles);
	}
	
	# Print the isolate information into the output file
	$line = $parts[0] . "\t" . $parts[1] . "\t";
	$line .= join(":", @isolateInfo);
	print OUTPUT "$line\n";
	
	# Progress information
	if($lineNo % 1000 == 0){
		print "Finished reading line: $lineNo\n";
	}
}

# Close the input and output files
close(FILTERED);
close(OUTPUT);