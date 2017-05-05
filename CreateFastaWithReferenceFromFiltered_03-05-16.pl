#!/usr/bin/perl
# Searching files in a directory
use warnings;
use strict;
use Term::ANSIColor; # For Coloured Print Statements

# Author: Joseph Crisp
# Creates a FASTA sequence file from the Filtered Variant positions information

# Command Line Run Structure:
# perl CreateFastaFromFiltered.pl verbose proximityFilter referenceSeq.fasta filtered.txt output.fasta

# Reference Genome FASTA file structure
# 	>gi|31791177|ref|NC_002945.3| Mycobacterium bovis AF2122/97 chromosome, complete genome
# 	TTGACCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACGCGGTCGTCTCCGAACTTAACGGCGACC
#	CTAAGGTTGACGACGGACCCAGCAGTGATGCTAATCTCAGCGCTCCGCTGACCCCTCAGCAAAGGGCTTG
#	...

# NOTE - doesn't remove non-informative sites - as these will differ against the reference!!

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
	print "\tperl CreateFastaWithReferenceFromFiltered.pl verbose proximityFilter referenceSeq.fasta filtered.txt output.fasta\n", color("reset");
	print "Proximity Filter\tIf any SNPs found within X positions of one another they are removed."
	
}else{

	# Get the proximity filter value
	my $proximityFilt = shift @ARGV;

	# Open the Reference Sequence file
	my $reference = shift @ARGV;
	print color("blue"), "Attempting to open Reference Sequence File...\n" if $verbose == 1;
	open REF, $reference or die "Couldn't open $reference:$!";
	
	# Initialise a variable to store the sequence
	my $sequence = "";
	
	# Initialise the necessary variables to parse the fasta file
	my $line;
	my @cols;
	
	# Begin reading the Reference sequence file
	print color("blue"), "Reading the Reference Genome File: $reference...\n" if $verbose == 1;
	while(<REF>){
		$line = $_;
		chomp($line);
		
		# Skip the header line
		next if $line =~ /^>gi/;
		
		# Store the line as the next part of the sequence
		$sequence = $sequence . $line;
	}
	close(REF);
	
	# Create an array of nucleotides from the fasta sequence
	my @ReferenceNucleotides = split //, $sequence;
	print color("blue"), "Finished Reading the Reference Genome File: $reference. Stored sequence as array of nucleotides.\n" if $verbose == 1;
	
	
	# Open the Filtered file as input
	my $filtered = shift @ARGV;
	print color("blue"), "Attempting to open Filtered Variants File...\n" if $verbose == 1;
	open FILTERED, $filtered or die "Couldn't open $filtered:$!";

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

	while(<FILTERED>){
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
			
			# Add an identifier for the Reference Sequence
			$isolateNames[scalar(@isolateNames)] = "Ref-1997";
			
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
		
		# Add the allele from the reference sequence
		$isolateSequences[$noIsolates][$sequencePos] = $ReferenceNucleotides[$cols[1] - 1];
	}
	close(FILTERED);

	print "Finished Reading Sequences.\n" if $verbose == 1;

	#######################################
	# Remove Sites Using Proximity filter #
	#######################################
	
	# Initialise an array to record wich sites to remove
	my @remove = (0) x ($sequencePos + 1);
	
	# Create a variable to count the number of sites removed using the proximity filter
	my $nRemovedByProximity = 0;
	
	for(my $i = 0; $i < ($sequencePos + 1); $i++){

		# Apply the proximity filtered
		if(($i != 0 && $positions[$i] - $positions[$i - 1] <= $proximityFilt) ||
		($i != $sequencePos && $positions[$i + 1] - $positions[$i] <= $proximityFilt)){
		
			$remove[$i] = 1;
			$nRemovedByProximity++;
		}
	}
		
	#####################
	# Create Fasta File #
	#####################

	####### Print out each Samples Sequence #######

	# Fasta File Structure
	# 	2 18					noSamples sequenceLength
	#	>Sample1ID
	#	GCTGATCGTGNGTAGGCC
	#	>Sample2ID
	#	GCTGATCGTGNGTAGGCT
	
	print "Opening Output File...\n" if $verbose == 1;

	# Open the Output File
	my $outputFastaFile = shift @ARGV; # Open Output File to Print out to
	open OUTPUT, ">$outputFastaFile" or die "Couldn't open $outputFastaFile:$!";
	
	print "Creating Fasta File...\n" if $verbose == 1;
	
	# Print the Sequences out to File - note two slightly different formats
	$noIsolates = $noIsolates + 1; # Remember we have added the reference
	my $sequenceLength = ($sequencePos + 1) - $nRemovedByProximity;
	print OUTPUT "$noIsolates $sequenceLength\n";
		
	for(my $pos = 0; $pos < $noIsolates; $pos++){
	
		# Print isolates file name
		print OUTPUT ">$isolateNames[$pos]\n";
	
		# Print isolate Sequence - only the informative sites
		for(my $i = 0; $i < ($sequencePos + 1); $i++){
			
			# Print if not filtered by proximity
			print OUTPUT "$isolateSequences[$pos][$i]" if $remove[$i] == 0; 
		}
		
		print OUTPUT "\n";
	}
	close(OUTPUT);
	
	print color("green"), "\t$nRemovedByProximity of ", $sequencePos + 1, " were removed by proximity filtering\n" if $verbose == 1;
	print color("green"), "Finished Creating Sequence Fasta File.\n", color ("reset") if $verbose == 1;
}
