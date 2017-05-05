#!/usr/bin/perl
# Searching files in a directory
use warnings;
use strict;
use Term::ANSIColor; # For Coloured Print Statements

# Author: Joseph Crisp
# Creates a FASTA sequence file from the Filtered Variant positions information

# Command Line Run Structure:
# perl CreateFastaFromFiltered.pl verbose proximityFilter filtered.txt output.fasta

# For each Variant Position within the filtered file the following information is available:
# #CHROM	POS	Sample 1:Sample 2:Sample 3:...      \t
# 0			1	2
#				|
#				 ---> 	DP	HQDP	MQ	AlleleSupport	AlleleCalled	QUAL	FQ		Result		;
# 						0 	1  		2 	3				5				6		7		8	

# Fasta File Structure
# 	2 18					noSamples sequenceLength
#	>Sample1ID
#	GCTGATCGTGNGTAGGCC
#	>Sample2ID
#	GCTGATCGTGNGTAGGCT


#########
# Input #
#########

my $verbose = shift @ARGV; # Open merged VCF File

# Print help information if needed
if($verbose eq "-help" || $verbose eq ""){
	
	print color("blue"), "Perl Script to create a Sequence FASTA file from the filtered Variant positions.\n";
	print "\nCommand Line Structure:\n";
	print "\tperl CreateFastaFromFiltered.pl verbose proximityFilter filtered.txt output.fasta\n", color("reset");
	print "Proximity Filter\tIf any SNPs found within X positions of one another they are removed."
	
}else{

	# Get the proximity filter value
	my $proximityFilt = shift @ARGV;
	print color("blue"), "Proximity Filter = $proximityFilt\n" if $verbose == 1;
	
	
	# Open the Filtered Variants file
	my $inputFile = shift @ARGV;
	print color("blue"), "Attempting to open Filtered File...\n" if $verbose == 1;
	open INPUT, $inputFile or die "Couldn't open $inputFile:$!";

	##################
	# Read Sequences #
	##################

	# Keep track of what position in FASTA sequence investigating
	my $sequencePos = -1; 
	
	# Initialise an array to store the position for each SNP
	my @positions = ();

	# Initialise the necessary variables to recording the sequence information
	my @isolateNames = ();
	my @isolateInfo = ();
	my $noIsolates;
	my @isolateSnpDetails = ();
	my @isolateSequences; # Matrix (Array of Sequence Arrays)
	

	# Begin reading the filtered file
	print "Reading sequences from Filtered File...\n" if $verbose == 1;	

	my $line;
	my @cols;
	while(<INPUT>){
		my $line = $_;
		chomp($line);
		next if $line =~ /^##/;
		
		# Split SNP Position Information Line into an Array
		@cols = split /\t/, $line; 
		
		####### Extract Sample IDs from Header Line #######
		#	#CHROM	POS	Sample 1:Sample 2:Sample 3:...
		# 	0 		1	2
		if($line =~ /^#/){
			@isolateNames = split /\:/, $cols[2]; # Extract the Sample Names into Array
			$noIsolates = @isolateNames;
			
			next;
		}	

		$sequencePos++; # Move to next Position in FASTA Sequence Array
		
		# Note the position for the current line
		$positions[$sequencePos] = $cols[1];
		
		####### Build each Sample Sequence #######

		@isolateInfo = split /\:/, $cols[2];
	
		# Investigate the SNP Position information from each sample
		for(my $pos = 0; $pos < $noIsolates; $pos++){
		
			# Insert an N where isolates have no information available
			if($isolateInfo[$pos] =~ /^----/){ # If no Information Available from Current isolate...
				$isolateSequences[$pos][$sequencePos] = "N";
				next; # Move to next Sample
			}	

			# Extract the Current Sample SNP Position Details:
			# DP	HQDP	MQ	AlleleSupport	AlleleCalled	QUAL	FQ		Result
			# 0		1		2	3				4				5		6		7
			@isolateSnpDetails = split /\;/, $isolateInfo[$pos];
		
			if($isolateSnpDetails[7] eq "Pass"){
				$isolateSequences[$pos][$sequencePos] = $isolateSnpDetails[4]; # If passed print out allele called
			}else{
				$isolateSequences[$pos][$sequencePos] = "N"; # If failed insert an "N" - no information available for this particular site
			}
		}
	}
	close(INPUT);

	print "Finished Reading Sequences.\n" if $verbose == 1;

	######################################
	# Note which sites are uninformative #
	######################################

	print "Finding Informative Sites...\n" if $verbose == 1;
	
	# Find the number of sites
	my $nSites = ($sequencePos + 1);

	# Initialise an array to note the alleles removed
	my @removedCounts = (0) x 4; # ACGT
	
	# Keep count of the sites removed by proximity filter
	my $nRemovedByProximity = 0;
	my $remove = 0;

	# Initialise an array to record which sites are informative
	my @informativeSites = (0) x $nSites;
	my $nInformativeSites = 0;

	# Examine each site
	for(my $siteIndex = 0; $siteIndex < $nSites; $siteIndex++){
		
		# Reset proximity check
		$remove = 0;
		
		# Check if site will be removed by proximity filter
		if(($siteIndex != 0 && $positions[$siteIndex] - $positions[$siteIndex - 1] <= $proximityFilt) ||
		($siteIndex != $sequencePos && $positions[$siteIndex + 1] - $positions[$siteIndex] <= $proximityFilt)){
			
			$nRemovedByProximity++;
			$remove = 1;
		}
		
		# Find the first site that isn't and N and compare to all others
		for(my $isolateI = 0; $isolateI < $noIsolates; $isolateI++){
	
			# Is the current isolate's site an N?
			if($isolateSequences[$isolateI][$siteIndex] ne "N"){
			
				# Only look if informative if site not removed by proximity filter
				if($remove == 0){
				
					# Compare the current isolate to all the other isolates at the current site, are the same or different?
					for(my $isolateJ = 0; $isolateJ < $noIsolates; $isolateJ++){
				
						# Are the sites different? - INFORMATIVE!
						if($isolateSequences[$isolateJ][$siteIndex] ne "N" && 
						$isolateSequences[$isolateI][$siteIndex] ne $isolateSequences[$isolateJ][$siteIndex]){
							$informativeSites[$siteIndex] = 1;
							$nInformativeSites++;
					
							# Finish comparing
							last;
						}				
					}
				}
			
				# If site is uninformative, record the allele removed - ACGT
				if($informativeSites[$siteIndex] == 0){
					if($isolateSequences[$isolateI][$siteIndex] eq "A"){
						$removedCounts[0]++;
					}elsif($isolateSequences[$isolateI][$siteIndex] eq "C"){
						$removedCounts[1]++;
					}elsif($isolateSequences[$isolateI][$siteIndex] eq "G"){
						$removedCounts[2]++;
					}elsif($isolateSequences[$isolateI][$siteIndex] eq "T"){
						$removedCounts[3]++;
					}else{
						print "ERROR: Site allele note recognised: $isolateSequences[$isolateI][$siteIndex]\n";
					}
				}
			
				# Finish looking at the current site
				last;
			}
		}
	}
	
	###########################
	# Print out the sequences #
	###########################
	
	# Open the output file
	my $outputFile = shift @ARGV;
	open OUTPUT, ">$outputFile" or die "Couldn't open $outputFile:$!";
	
	# Print out a header = NumberIsolates	SequenceLength
	print OUTPUT "$noIsolates $nInformativeSites\n";
	
	# Initialise a variable to store a sequence
	my $sequence = "";

	# Print out the sequences
	for(my $isolateIndex = 0; $isolateIndex < $noIsolates; $isolateIndex++){
	
		print OUTPUT ">$isolateNames[$isolateIndex]\n";
		$sequence = "";

		for(my $siteIndex = 0; $siteIndex < $nSites; $siteIndex++){
		
			# Skip non-informative sites
			next if $informativeSites[$siteIndex] == 0;
		
			# Add the site allele to the current individuals sequence
			$sequence = $sequence . $isolateSequences[$isolateIndex][$siteIndex];
		}
	
		print OUTPUT "$sequence\n";
	}
	close(OUTPUT);

	#################################################
	# Print out information about the sites removed #
	#################################################

	print "$nRemovedByProximity of $nSites were removed by proximity filtering\n" if $verbose == 1;
	print "$nInformativeSites informative sites were retained.\n" if $verbose == 1;
	print "Allele Counts of removed sites:\n" if $verbose == 1;
	print "A\tC\tG\tT\n" if $verbose == 1;
	print "$removedCounts[0]\t$removedCounts[1]\t$removedCounts[2]\t$removedCounts[3]\n" if $verbose == 1;
	
	print color("green"), "Finished Creating Sequence Fasta File.\n", color ("reset") if $verbose == 1;
}
