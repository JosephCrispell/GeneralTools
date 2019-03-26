#!/usr/bin/perl
# Searching files in a directory
use warnings;
use strict;
use Term::ANSIColor; # For Coloured Print Statements

# Author: Joseph Crisp
# Filtering the Recorded nucleotides according to Quality Metrics:
#		Read Depth
#		High Quality Base Depth
#		Mapping Quality
#		Allele Support
#		Proportion Isolates with Sufficient Coverage
#		Quality Score (QUAL)
#		FQ

# Command Line Run Structure: Should be ran in Directory with VCF Files
# perl FilterVariants.pl verbose DP HQDP MQ PRA COV QUAL FQ merged.txt

# MergedVCFs File Sample Information:
# #CHROM	Pos	Sample 1:Sample 2:Sample 3:...					\t
# 0			1 	2  		
# 				|		 
#				 ---->	DP	HQDP	MQ	QUAL	FQ	Ref	Alt		;
# 						0 	1  		2 	3		4	5	6

# Creates file where information for each SNP position is recorded (For each of the samples available):
# #CHROM	POS	Sample 1:Sample 2:Sample 3:...
# 0			1	2
#				|
#				 ---> 	DP	HQDP	MQ	AlleleSupport	AlleleCalled	QUAL	FQ		Result		Ref	Alt	;
# 						0 	1  		2 	3				4				5		6		7			8	9

###########
# Filters #
###########

# Get the Filter thresholds from the Command Line options
my $verbose = shift @ARGV;

# Print out Help instructions if needed
if($verbose eq "-help" || $verbose eq ""){
	print color("green"), "Perl Script to filter the Variant Positions in the Merged VCF File.\n\nCommand Line Structure:\n", color("reset");
	print "\tperl FilterVariants.pl verbose DP HQDP MQ PRA COV QUAL FQ merged.txt\n";
	print "\t\tverbose\t\tProvide 1 if want detailed output to terminal, 0 if not\n";
	print "\t\tDP\t\tRead Depth\n";
	print "\t\tHQDP\t\tHigh Quality Base Depth\n";
	print "\t\tMQ\t\tMapping Quality\n";
	print "\t\tPRA\t\tProportion of aligned Reads supporting Allele called\n";
	print "\t\tCOV\t\tProportion sites with sufficient coverage across isolates\n";
	print "\t\tQUAL\t\tGeneral Quality score\n";
	print "\t\tFQ\t\tScore that relates to the probability of multiple alleles being present across samples (irrelevant for single samples?)\n";
	print "\t\tmerged.txt\t\tPath to merged VCFs file\n";
	
}else{

	my $dpFilter = shift @ARGV;
	my $hqdpFilter = shift @ARGV;
	my $mqFilter = shift @ARGV;
	my $supportFilter = shift @ARGV;
	my $propCovFilter = shift @ARGV;
	my $qualFilter = shift @ARGV;
	my $fqFilter = shift @ARGV;
	
	print "Filter Selection:\n" if $verbose == 1;
	print "\tDP = $dpFilter\n\tHQDP = $hqdpFilter\n\tMQ = $mqFilter\n\tSUP = $supportFilter\n\tCOV = $propCovFilter\n\tQUAL = $qualFilter\n\tFQ = $fqFilter\n" if $verbose == 1;

	###############
	# Preparation #
	###############

	####### Get Output File for Filtered Nucleotide Positions #######
	
	my $mergedFile = shift @ARGV; 
	open INPUT, $mergedFile or die "Couldn't open $mergedFile:$!"; # Attempt to open Current VCF File
	
	print "\nInput merged VCF file: $mergedFile\n" if $verbose == 1;
	
	# Get the date from the merged VCF file name
	my $date = substr((split /_/, $mergedFile)[1], 0, -4);
	
	# Create and open the output files
	my $outputFilterFile = "filtered_" . $dpFilter . "-" . $hqdpFilter . "-" . $mqFilter . "-" . $supportFilter . "-" . $propCovFilter . "-" . $qualFilter . "-" . $fqFilter . "_" . $date . ".txt"; 
	open OUTPUTFILT, ">$outputFilterFile" or die "Couldn't open $outputFilterFile:$!"; # Assigns OUTPUTFILT File to Label and checks it can be Opened
	my $outputCoverageFile = "snpCoverage_" . $date . ".txt";
	open OUTPUTCOV, ">$outputCoverageFile" or die "Couldn't open $outputCoverageFile:$!"; # Assigns OUTPUTCOV File to Label and checks it can be Opened

	print "\nThis script produces the following output files:\n" if $verbose == 1;
	print "\t$outputFilterFile\t\tA filtered version of the input merged VCFs file.\n" if $verbose == 1;
	print "\t$outputCoverageFile\t\t\t\tA file detailing the coverage of each site in the merged VCF file across the isolates\n\n" if $verbose == 1;
	
	print OUTPUTCOV "#CHROM\tPOS\tPropIsolatesWithSuffCov\tPropIsolateWithPosWithSuffCov\tNoIsolatesWithPos\n";
	
	print color("blue"), "Output File Checked. \n" if $verbose == 1;

	###################
	# Begin Filtering #
	###################

	print "Beginning Variant Filtering...\n" if $verbose == 1;
	
	# DP	HQDP	MQ	QUAL	FQ	Ref 	Alt		;
	# 0 	1  		2 	3		4	5		6

	# Initialise the necessary variables
	my @snpInfo;
	my @isolateInfo;
	my $readDepth;
	my @hqBaseDepth;
	my $mappingQuality;
	my $qual;
	my $fq;
	my $ref;
	my $alt;
	
	# Initialise variables needed for further calculations
	my $refProportion; # Proportion of HQ bases supporting Reference allele
	my $altProportion; 
	my $call;
	my @alleleHQDepth = (0) x 2;
	
	# Initialise variables to keep track of which isolates have coverage
	my $noIsolatesWithSufficientCoverage;
	my $noIsolatesWithPosNotPresent;
	my $noIsolatesWithPosPresent;
	my $proportionIsolatesWithSufficientCoverage;
	my $proportionIsolatesWithSufficientCoverageWherePresent;
	
	# Initialise an array to note the alleles of the sites removed
	my @nAlellesOfSitesRemoved = (0) x 4;
	
	# Initialise variables to store output
	my $result;
	my $resultInfo;
	my $outputLine;

	# Keep a count of the number of passed Variant Positions
	my $noPassedVPs = 0; # Note this counts VPs where at least one isolate passed the filtering criteria
	my $noVPs = 0;
	
	# Begin Reading the MergedVCFs file
	my $line;
	my @cols = ();
	while(<INPUT>){ # Read in each Line of Current VCF File
		$line = $_;
		chomp($line);
	
		# Print the VCF Header and Fields to File (Fields don't change)
		if($line =~ /^#/){
			print OUTPUTFILT "$line\n";
			next;
		}

		# Split the current line: #CHROM	POS	Sample 1:Sample 2:Sample 3:...
		@cols = split /\t/, $line;
		$outputLine = $cols[0] . "\t" . $cols[1] . "\t";
		$noVPs++;
	
		# Extract into Array SNP information for each Sample	
		@snpInfo = split /\:/, $cols[2];
	
		# Examine the SNP information from each of the isolates
		$noIsolatesWithSufficientCoverage = 0;
		$noIsolatesWithPosNotPresent = 0;
		for(my $pos = 0; $pos < scalar(@snpInfo); $pos++){
		
			# If there is no information available for the current sample then can't filter
			if($snpInfo[$pos] =~ /^----/){
				$outputLine = $outputLine . $snpInfo[$pos];
				$outputLine = $outputLine . "\:" if $pos < scalar(@snpInfo) - 1;
				$noIsolatesWithPosNotPresent++;
				next;
			}
			
			# Extract the Current Samples SNP Information
			# DP	HQDP	MQ	QUAL	FQ	Ref	Alt		;
			# 0 	1  		2 	3		4	5	6
			@isolateInfo = split /\;/, $snpInfo[$pos];
					
			# Store the quality Metrics
			$readDepth = $isolateInfo[0];
			@hqBaseDepth = split /\,/, $isolateInfo[1]; # ReferenceForward, ReferenceReverse, AlternateForward, AlternateReverse
			$mappingQuality = $isolateInfo[2];
			$qual = $isolateInfo[3];
			$fq = $isolateInfo[4];
			$ref = $isolateInfo[5];
			$alt = $isolateInfo[6];
	
			# Calculate the Proportion of High Quality Bases supporting the Reference Allele
			$refProportion = 0;
			if($hqBaseDepth[0] + $hqBaseDepth[1] != 0){
				$refProportion = ($hqBaseDepth[0] + $hqBaseDepth[1]) / ($hqBaseDepth[0] + $hqBaseDepth[1] + $hqBaseDepth[2] + $hqBaseDepth[3]);
			}
		
			# Calculate the Proportion of High Quality Bases supporting the Alternate Allele
			$altProportion = 0;
			if($hqBaseDepth[2] + $hqBaseDepth[3] != 0){
				$altProportion = ($hqBaseDepth[2] + $hqBaseDepth[3]) / ($hqBaseDepth[0] + $hqBaseDepth[1] + $hqBaseDepth[2] + $hqBaseDepth[3]);
			}
		
			# Allocate Nucleotide Position - Either Ref or Alt Position
			$call = "NA"; 
			if($refProportion >= $altProportion){
				$call = "Ref";
				$alleleHQDepth[0] = $hqBaseDepth[0];
				$alleleHQDepth[1] = $hqBaseDepth[1];
			}else{
				$call = "Alt";
				$alleleHQDepth[0] = $hqBaseDepth[2];
				$alleleHQDepth[1] = $hqBaseDepth[3];
			}
			
		
			####### Apply the Filters #######
			$result = "Fail";
			
			# READ DEPTH
			if($readDepth >= $dpFilter){
				
				# HIGH QUALITY BASE DEPTH
				if($alleleHQDepth[0] >= $hqdpFilter && $alleleHQDepth[1] >= $hqdpFilter){
					
					# MAPPING QUALITY
					if($mappingQuality >= $mqFilter){
						
						# ALLELE SUPPORT
						if(($call eq "Ref" && $refProportion >= $supportFilter) || ($call eq "Alt" && $altProportion >= $supportFilter)){
							
							# QUALITY SCORE 
							if($qual >= $qualFilter){
									
								# FQ (FQ is negative)
								if($fq <= $fqFilter){
									
									$result = "Pass";
									$noIsolatesWithSufficientCoverage++;
								}
							}
						}
					}					
				}
			}
			
			# Where an Alternate position has passed but there is no alternate allele present set result to FAIL
			if($result eq "Pass" && $call eq "Alt" && $alt eq "."){
				$result = "Fail";
			}
		
			# Prepare the SNP filter information from the current sample DP	HQDP	MQ	AlleleSupport	AlleleCalled	QUAL	FQ		Result	Ref	Alt
			$resultInfo =  $readDepth . "\;" . $isolateInfo[1] . "\;" . $mappingQuality . "\;";
			if($call eq "Ref" && $result eq "Pass"){
				$resultInfo = $resultInfo . $refProportion . "\;" . $ref . "\;";
				
			}elsif($call eq "Alt" && $result eq "Pass"){
				$resultInfo = $resultInfo . $altProportion . "\;" . $alt . "\;";
				
			}else{
				$resultInfo = $resultInfo . $altProportion . "\;" . "N" . "\;";
			}
			$resultInfo = $resultInfo . $qual . "\;" . $fq . "\;" . $result . ";" . $ref . ";" . $alt;

			# Store the Filtering Results for the current isolate
			$outputLine = $outputLine . "$resultInfo";
		
			# Add in a separator if haven't reached the last isolate
			$outputLine = $outputLine . "\:" if $pos < scalar(@snpInfo) - 1;
		}	
	
		# SNP COVERAGE ACROSS ISOLATES
		$proportionIsolatesWithSufficientCoverage = $noIsolatesWithSufficientCoverage / scalar(@snpInfo);
		if( $proportionIsolatesWithSufficientCoverage >= $propCovFilter){
			print OUTPUTFILT "$outputLine\n";
			
			$noPassedVPs++ if $outputLine =~ /Pass/;
		
		# Record the allele of the site removed
		}else{
			
			if($ref eq 'A'){
				$nAlellesOfSitesRemoved[0]++;
			}elsif($ref eq 'C'){
				$nAlellesOfSitesRemoved[1]++;
			}elsif($ref eq 'G'){
				$nAlellesOfSitesRemoved[2]++;
			}elsif($ref eq 'T'){
				$nAlellesOfSitesRemoved[3]++;
			}elsif(length($ref) == 1){ # Check that we haven't found indel
				print "ERROR!: Reference allele not recognised: $ref\n";
			}
		}
		$noIsolatesWithPosPresent = scalar(@snpInfo) - $noIsolatesWithPosNotPresent;
		$proportionIsolatesWithSufficientCoverageWherePresent = 0;
		$proportionIsolatesWithSufficientCoverageWherePresent = $noIsolatesWithSufficientCoverage / $noIsolatesWithPosPresent if $noIsolatesWithPosPresent != 0;
		
		
		# Print the Proportion of Isolates that have coverage at the current SNP
		print OUTPUTCOV "$cols[0]\t$cols[1]\t$proportionIsolatesWithSufficientCoverage\t$proportionIsolatesWithSufficientCoverageWherePresent\t$noIsolatesWithPosPresent\n";

	} # END
	close(INPUT);
	close(OUTPUTFILT);
	close(OUTPUTCOV);
	
	print color("green"), "Variant Filtering Complete.\n", color("reset") if $verbose == 1;
	my $noFiltered = $noVPs - $noPassedVPs;
	print "\t$noPassedVPs Variant Positions of $noVPs passed the filtering criteria ($noFiltered).\n";
	my $nFailedSites = $noVPs - $noPassedVPs;
	print "\t$nFailedSites sites that failed across the isolates were retained in merged file since coverage filter was set to zero.\n" if $verbose == 1 && $propCovFilter == 0;
	
	print "\nAllele Counts of removed sites:\n" if $verbose == 1;
	print "A\tC\tG\tT\n" if $verbose == 1;
	print "$nAlellesOfSitesRemoved[0]\t$nAlellesOfSitesRemoved[1]\t$nAlellesOfSitesRemoved[2]\t$nAlellesOfSitesRemoved[3]\n" if $verbose == 1;
}

