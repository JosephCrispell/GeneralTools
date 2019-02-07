#!/usr/bin/perl
# Searching files in a directory
use warnings;
use strict;
use Term::ANSIColor; # For Coloured Print Statements

# Author: Joseph Crisp
# Creates a FASTA sequence file from the Filtered Variant positions information

# Command Line Run Structure:
# perl CreateFastaFromFiltered.pl verbose proximityFilter referenceSeq.fasta filtered.txt

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
#				 ---> 	DP	HQDP	MQ	AlleleSupport	AlleleCalled	QUAL	FQ		Result		Ref	Alt	;
# 						0 	1  		2 	3				4				5		6		7			8	9

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
	print "\tperl CreateFastaWithReferenceFromFiltered.pl verbose proximityFilter reference.fasta filtered.txt\n", color("reset");
	print "\t\tverbose\t\t\tProvide 1 if want detailed output to terminal, 0 if not\n";
	print "\t\tProximity Filter\tIf any SNPs found within X positions of one another they are removed.\n";
	print "\t\treference.fasta\t\tPath to the M. bovis reference genome FASTA file\n";
	print "\t\tfiltered.txt\t\tPath to the filtered VCF file.\n";
	print "\t\t*** NOTE *** Currently inserts \"N\" for an INDEL\n";

	exit;
}
	

# Get the input settings
my $proximityFilt = shift @ARGV;
my $reference = shift @ARGV;
my $filtered = shift @ARGV;

# Get the date from the filtered file
my $date = substr((split /_/, $filtered)[2], 0, -4);

# Create the name for the output FASTA file
my $outputFastaFile = "sequences_Prox-" . $proximityFilt . "_" . $date . ".fasta"; # Open Output File to Print out to
my $outputFastaPositionsFile = "fastaPositions_Prox-" . $proximityFilt . "_" . $date . ".txt"; # Open Output File to Print out to

# Print a check of the input information
print "\nInput Settings:\n" if $verbose == 1;
print "Proximity filter = $proximityFilt\n" if $verbose == 1;
print "Reference sequence file: $reference\n" if $verbose == 1;
print "Filtered VCF file: $filtered\n" if $verbose == 1;
print "\nThe following output files are produced: \n" if $verbose == 1;
print "\t$outputFastaFile\n" if $verbose == 1;
print "\t$outputFastaPositionsFile\n" if $verbose == 1;

# Open the Reference Sequence file
open REF, $reference or die "Couldn't open $reference:$!";

# Initialise a variable to store the sequence
my $sequence = "";

# Initialise the necessary variables to parse the fasta file
my $line;
my @cols;

# Begin reading the Reference sequence file
print color("blue"), "Reading the Reference Genome File...\n" if $verbose == 1;
while(<REF>){
	
	$line = $_;
	chomp($line);
	
	# Skip the header line
	next if $line =~ /^>/;
	
	# Store the line as the next part of the sequence
	$sequence = $sequence . $line;
}
close(REF);

# Create an array of nucleotides from the fasta sequence
my @referenceNucleotides = split //, $sequence;
	
# Open the Filtered file as input
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
my @indelCounts = ();

# Begin reading the filtered file
print "Reading sequences from Filtered File...\n", color ("reset") if $verbose == 1;	

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

	# Initialise a counter for INDELS or Ns
	$indelCounts[$sequencePos] = 0;

	# Investigate the SNP Position information from each sample
	for(my $pos = 0; $pos < $noIsolates; $pos++){
	
		# Insert an N where isolates have no information available
		if($isolateInfo[$pos] =~ /^----/){ # If no Information Available from Current isolate...
			$isolateSequences[$pos][$sequencePos] = "N";
			next; # Move to next Sample
		}	

		# Extract the Current Sample SNP Position Details:
		# DP	HQDP	MQ	AlleleSupport	AlleleCalled	QUAL	FQ		Result	Ref	Alt
		# 0		1		2	3				4				5		6		7		8	9
		@isolateSnpDetails = split /\;/, $isolateInfo[$pos];
		
		# Insert an N if the current isolate has an INDEL at the this position. Will be an INDEL if the allele called has more than 1 nucleotide
		if(length($isolateSnpDetails[4]) > 1){
			$isolateSequences[$pos][$sequencePos] = "N";
			$indelCounts[$sequencePos]++;
			next; # Move to next Sample
		}

		# Insert the current sample's allele
		if($isolateSnpDetails[7] ne "Fail"){
			$isolateSequences[$pos][$sequencePos] = $isolateSnpDetails[4]; # If passed print out allele called
		}else{
			$isolateSequences[$pos][$sequencePos] = "N"; # If failed insert an "N" - no information available for this particular site
		}
	}
		
	# Add the allele from the reference sequence
	$isolateSequences[$noIsolates][$sequencePos] = $referenceNucleotides[$cols[1] - 1];
}
close(FILTERED);

#######################################
# Remove Sites Using Proximity filter #
#######################################
	
# Initialise an array to record wich sites to remove
my @remove = (0) x ($sequencePos + 1);
	
# Create a variable to count the number of sites removed using the proximity filter
my $nRemovedByProximity = 0;
my $nRemovedAsINDELs = 0;
my @nAlellesOfSitesRemoved = (0) x 4;
	
for(my $i = 0; $i < ($sequencePos + 1); $i++){

	# Check if only INDELs found at the current position
	if($indelCounts[$i] == $noIsolates){
		$remove[$i] = 1;
		print "$indelCounts[$i] INDELs found at $i\n" if $verbose == 1;
		$nRemovedAsINDELs++;
		next
	}elsif($indelCounts[$i] > 1){
		print "$indelCounts[$i] INDELs found at $i\n" if $verbose == 1;
	}

	# Apply the proximity filtered
	if(($i != 0 && $remove[$i - 1] == 0 && $positions[$i] - $positions[$i - 1] <= $proximityFilt) ||
	($i != $sequencePos && $positions[$i + 1] - $positions[$i] <= $proximityFilt)){
		
		$remove[$i] = 1;
		$nRemovedByProximity++;
			
		if($referenceNucleotides[$i] eq 'A'){
			$nAlellesOfSitesRemoved[0]++;
		}elsif($referenceNucleotides[$i] eq 'C'){
			$nAlellesOfSitesRemoved[1]++;
		}elsif($referenceNucleotides[$i] eq 'G'){
			$nAlellesOfSitesRemoved[2]++;
		}elsif($referenceNucleotides[$i] eq 'T'){
			$nAlellesOfSitesRemoved[3]++;
		}else{
			print "ERROR!: Reference allele not recognised: $referenceNucleotides[$i]\n";
		}
	}
}
	
##################################
# Report the sites used in FASTA #
##################################
	
# Open the Output File
open POSITIONS, ">$outputFastaPositionsFile" or die "Couldn't open $outputFastaPositionsFile:$!";
print color ("blue"), "Recording FASTA positions ...\n" if $verbose == 1;
	
# Print the sequence positions used in FASTA
print POSITIONS "Position\n";
for(my $i = 0; $i < ($sequencePos + 1); $i++){
		
	# Print position if not filtered by proximity
	print POSITIONS "$positions[$i]\n" if $remove[$i] == 0; 
}
close(POSITIONS);
	
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

# Open the Output File
open FASTA, ">$outputFastaFile" or die "Couldn't open $outputFastaFile:$!";
	
print "Creating Fasta File...\n" if $verbose == 1;
	
# Print the Sequences out to File - note two slightly different formats
$noIsolates = $noIsolates + 1; # Remember we have added the reference
my $sequenceLength = ($sequencePos + 1) - ($nRemovedByProximity + $nRemovedAsINDELs);
print FASTA "$noIsolates $sequenceLength\n";
	
for(my $pos = 0; $pos < $noIsolates; $pos++){
	
	# Print isolates file name
	print FASTA ">$isolateNames[$pos]\n";
	
	# Print isolate Sequence - only the informative sites
	for(my $i = 0; $i < ($sequencePos + 1); $i++){
		
		# Print if not filtered by proximity
		print FASTA "$isolateSequences[$pos][$i]" if $remove[$i] == 0; 
	}
		
	print FASTA "\n";
}
close(FASTA);
	
print color("green"), "\t$nRemovedByProximity of ", $sequencePos + 1, " were removed by proximity filtering\n" if $verbose == 1;
print color("green"), "\t$nRemovedAsINDELs of ", $sequencePos + 1, " were removed as all isolates were INDELs\n" if $verbose == 1;
print color("green"), "Finished Creating Sequence Fasta File.\n", color ("reset") if $verbose == 1;

print "\nAllele Counts of removed sites:\n" if $verbose == 1;
print "A\tC\tG\tT\n" if $verbose == 1;
print "$nAlellesOfSitesRemoved[0]\t$nAlellesOfSitesRemoved[1]\t$nAlellesOfSitesRemoved[2]\t$nAlellesOfSitesRemoved[3]\n" if $verbose == 1;

