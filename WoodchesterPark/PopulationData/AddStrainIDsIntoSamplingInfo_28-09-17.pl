#!/usr/bin/perl
# Searching files in a directory
use warnings;
use strict;
use Term::ANSIColor; # For Coloured Print Statements

# Author: Joseph Crisp
# Add strain IDs - associated with FASTA sequences - into sampling information file

# Command Line Run Structure:
# perl AddStrainIDsIntoSamplingInfo.pl sequences.fasta linkFile.csv cattleSamplingInfo.csv

#############
# FUNCTIONS #
#############

sub getSequenceIDs{
	
	# Get the FASTA file name
	my $fastaFile = $_[0];
	open FASTA, $fastaFile or die "Couldn't open $fastaFile:$!";

	# Initialise a hashtable to store the cattle sequence IDs
	my %sequenceIDs = ();
	my $sequenceID;

	# Initialise the necessary variables for parsing the file
	my $line;

	# Read the FASTA file line by line
	while(<FASTA>){
		$line = $_;
		chomp($line);
	
		# Skip none sequence header lines and badgers and reference
		next if $line !~ /^>/ || $line =~ /^>WB/ || $line =~ /^>Ref/;
	
		# Get the sequence name from the current line
		$sequenceID = substr((split /_/, $line)[0], 1);
		$sequenceIDs{$sequenceID} = 1;
	}
	close(FASTA);
	
	return(%sequenceIDs)
}

sub noteSequenceIDsAssociatedWithStrainIDs{

	# Get the input file
	my $linkFile = $_[0];

	# Open the link file
	open LINK, $linkFile or die "Couldn't open $linkFile:$!";

	# Initialise a hashtable to store the sequence IDs associated with each strain ID
	my %strainIDs = ();

	# Initialise the necessary variables for parsing the file
	my $line;
	my @cols;

	# Read the file line by line
	while(<LINK>){
		$line = $_;
		chomp($line);
		
		# Remove any special characters
		$line =~ s/\R//g;
	
		# Skip the header line
		next if $line =~ /^Seq/;
	
		# Split the current line into it's columns
		@cols = split /,/, $line;
	
		# Check if encountered strain ID already
		if(exists($strainIDs{$cols[1]})){
			push(@{$strainIDs{$cols[1]}}, $cols[0]);
		}else{
			@{$strainIDs{$cols[1]}} = ($cols[0]);
		}
	}
	close(LINK);
	
	return(%strainIDs);
}

sub readAndPrintSampleInfoWithSequenceIDsPresentInfFASTA{

	# Get the sampling information file
	my $sampleInfo = $_[0];
	
	# Get the strain IDs hashtable
	my %strainIDs = %{$_[1]};
	
	# Get the sequence IDs hashtable
	my %sequenceIDs = %{$_[2]};

	# Open the sampling information file
	open INFO, "$sampleInfo" or die "Couldn't open $sampleInfo:$!";

	# Create and and open an output file
	my $output = substr($sampleInfo, 0, -4) . "_AddedStrainIDs.csv";
	print "Creates the following output file: $output\n";
	open OUTPUT, ">$output" or die "Couldn't open $output:$!";

	# Initialise the necessary variables for parsing the file
	my $line;
	my @cols;

	# Initialise a variable to store the sequence IDs associated with a strain ID
	my @seqIDs;

	# Read the file in line by line
	while(<INFO>){
		$line = $_;
		chomp($line);
	
		# Skip the header line but print it out
		if($line =~ /^CPH/){
			print OUTPUT "$line\n";
			next;
		}
		
		# Split the current line inot its columns
		@cols = split /,/, $line;
	
		# Get the sequence IDs associated with the current string ID
		if(exists($strainIDs{$cols[1]})){
			@seqIDs = @{$strainIDs{$cols[1]}};
		
			# Select sequence ID which was present in FASTA file and insert into correct column
			for(my $i = 0; $i < scalar(@seqIDs); $i++){
				
				if(exists($sequenceIDs{$seqIDs[$i]})){
					$cols[-3] = $seqIDs[$i];
				}
			}
		}
	
		# Print the current line back out to file
		$line = join ",", @cols;
		print OUTPUT "$line\n";
	}
	close(INFO);
	close(OUTPUT);	
}

##################
# Check for help #
##################

my $firstArg = shift @ARGV;

if($firstArg eq "-help" || $firstArg eq "-" || $firstArg eq "help" || $firstArg eq "h" || $firstArg eq ""){
	print "Perl Script to add strain IDs into cattle sample information file\n";
	print "\tCommand line structure:\n";
	print "\t\tperl AddStrainIDsIntoSamplingInfo.pl sequences.fasta linkFile.csv cattleSamplingInfo.csv\n";
	print "\t\tsequences.fasta\tPath to file containing cattle sequences\n";
	print "\t\tlinkFile.csv\tPath to file listing the Strain IDs associated with each Sequence ID\n";
	print "\t\tcattleSamplingInfo.csv\tPath to cattle sampling information CSV\n";
	
	exit 0
}

###################################################
# Get list of cattle sequence IDs from FASTA file #
###################################################

# Get the FASTA file name
my $fastaFile = $firstArg;

# Get the cattle strain IDs and store in a hashtable
my %sequenceIDs = getSequenceIDs($fastaFile);

##########################################
# Note the Strain IDs for each Strain ID #
##########################################

# Get the link file name
my $linkFile = shift @ARGV;

# Store the sequence IDs associated with each strain ID in a hashtable
my %strainIDs = noteSequenceIDsAssociatedWithStrainIDs($linkFile);

###############################################################
# Add appropriate sequence IDs into sampling information file #
###############################################################

# Get the sampling information file
my $sampleInfo = shift @ARGV;

# Read the file and print it out adding in the appropriate sequence ID where necessary
readAndPrintSampleInfoWithSequenceIDsPresentInfFASTA($sampleInfo, \%strainIDs, \%sequenceIDs);