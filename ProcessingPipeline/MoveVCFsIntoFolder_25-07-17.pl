#!/usr/bin/perl
# Searching files in a directory
use warnings;
use strict;

# Author: Joseph Crisp
# Moves VCF files (identified by input file) into specified directory

# Command Line Structure
# perl MoveVCFsIntoFolder.pl isolatesToRemove.txt directory

# Isolates to remove file structure:
# Isolate		Reason
# Identifier	PoorVariantPositionCoverage

################################
# Get the command line options #
################################

# Open the isolates to remove file
my $isolatesToRemoveFile = shift @ARGV;
open INPUT, $isolatesToRemoveFile or die "Couldn't open $isolatesToRemoveFile:$!";

# Note the directory name where the VCF files are to be move to
my $directory = shift @ARGV;

#####################################################
# Note the identifiers of each VCF file to be moved #
#####################################################

# Initialise an array to store the identifiers
my @vcfIdentifiers = ();
my $pos = -1;

# Initialise the necessary variables for parsing the file
my $line;

# Begin reading the file
while(<INPUT>){
	$line = $_;
	chomp($line);
	
	# Skip the header
	next if $line =~ /^Isolate/;
	
	# Store the identifier
	$pos++;
	$vcfIdentifiers[$pos] = (split /\t/, $line)[0];
}
close(INPUT);

###################################################
# move the VCF files into the specified directory #
###################################################

# Initialise a variable to store the VCF file name associated with each identifier
my $file;

# Initialise a variabel to catch any errors
my $error;

foreach my $identifier (@vcfIdentifiers){

	# Get the associated file
	$file = `ls | grep $identifier`;
	chomp($file);
	
	# Move the file to the specified directory
	$error = `mv $file $directory`;
		
	# Check if there was an error
	if($error ne ""){
		print "ERROR: There was an error when move the file: $file\n\t$error\n";
	}else{
		print "Successfully moved $file to $directory\n";
	}
}
