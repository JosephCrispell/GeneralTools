#!/usr/bin/perl
# Searching files in a directory
use warnings;
use strict;
use List::MoreUtils qw(uniq);

# Author: Joseph Crisp
# Examine the BLAST output resulting from running a FASTA file against a local database
# Output from the following command:
# blastn -query example.fasta -out results.txt -db Database/databaseName -outfmt 0

# Command Line Structure
# perl ExamineBLASTOutput.pl nRecords outputBLAST.txt

# BLAST output file structure:
# 	BLASTN 2.2.25+
#	\n 
# 	Reference: Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb
# 	Miller (2000), "A greedy algorithm for aligning DNA sequences", J
# 	Comput Biol 2000; 7(1-2):203-14.
# 	\n\n 
# 	Database: proteinCodingSequences.fasta
#   	         2,127 sequences; 1,868,325 total letters
# 	\n\n 
# 	Query= .15056_6_10.1
# 
# 	Length=387784
#                                                                       Score     E
# 	Sequences producing significant alignments:                          (Bits)  Value
# 
# 	lcl|1013528                                                           6536    0.0   
# 	lcl|1013288                                                           5145    0.0   
# 	...
#	\n
#	>lcl|1013528 
#	Length=3539
#
#	 Score = 6536 bits (3539),  Expect = 0.0
#	 Identities = 3539/3539 (100%), Gaps = 0/3539 (0%)
#	 Strand=Plus/Plus
#
#	Query  308962  ATGTTTTTAAAAGAAATTGAAATGCAAGGCTTTAAGTCCTTTGCGGATAAAACTAAGGTA  309021
#	               ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#	Sbjct  1       ATGTTTTTAAAAGAAATTGAAATGCAAGGCTTTAAGTCCTTTGCGGATAAAACTAAGGTA  60
#	...
#
#	>lcl|1013288 
#	...
#	\n\n	
#	Lambda     K      H
#	    1.33    0.621     1.12 
#
#	Gapped
#	Lambda     K      H
#	    1.28    0.460    0.850 
#
#	Effective search space used: 2761440996
#	\n
#	  Database: proteinCodingSequences.fasta
#	    Posted date:  Aug 24, 2015  2:48 PM
#	  Number of letters in database: 1,868,325
#	  Number of sequences in database:  2,127
#	\n\n
#	Matrix: blastn matrix 1 -2
#	Gap Penalties: Existence: 0, Extension: 2.5

# Output line structure:
#FileName:GeneId(Score,E),GeneId(Score,E),...

# Get the number of records to be retained
my $nRecords = shift @ARGV;
my $recordsRetained = 0;

# Open the BLAST output file
my $blastFile = shift @ARGV;
open BLAST, $blastFile or die "Couldn't open $blastFile:$!";

# Get the query FASTA file name from the command line
my $query = shift @ARGV;

# Initialise an array to store the IDs of all the significant hits
my @hits = ();
my $index = -1;

# Initialise additional variables for tracking position in file
my $foundTable = 0;

# Initialise necessary variables for parsing a file
my $line;
my $output;
my @parts;
my $hitInfo;

# Begin reading the file
while(<BLAST>){
	$line = $_;
	chomp($line);
	
	if($line =~ /Value$/){
		$foundTable = 1;
		
	}elsif($foundTable == 1 && $line =~ /^>/){
		$foundTable = 0;
		$recordsRetained = 0;
	}
	elsif($foundTable == 1 && $line !~ /^>/ && $line ne ""){
		
		$recordsRetained++;
		print "$line\n" if $recordsRetained <= $nRecords;	
	}
}

