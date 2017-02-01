#!/usr/bin/perl
# Searching files in a directory
use warnings;
use strict;

# Author: Joseph Crisp
# Extract the sampling information needed to accompany the FASTQ files on NCBI

# Command Line Structure
# perl ExtractMetadataForNCBIUpload.pl isolatesUsed.csv output.txt

# IsolatesUsed.csv file structure:
# IsolateID	Year	Host	Region	Area	REA
# 0			1		2		3		4		5

# NCBI Metadata file:
# *sample_name	sample_title	bioproject_accession	*organism	strain	isolate	host	isolation_source	*collection_date
# 0 ###			1				2						3 ###		4		5 ###	6 ###	7					8 ###

# *geo_loc_name	*sample_type	altitude	biomaterial_provider	collected_by	culture_collection	depth	env_biome	genotype	
# 9	###			10 ###			11			12						13				14					15		16			17

# host_tissue_sampled	identified_by	lab_host	lat_lon	mating_type	passage_history	samp_size	serotype	serovar	specimen_voucher	
# 18					19				20			21					22				23			24			25		26

# temp	description
# 27	28

#############################################################
# Read in the Isolates Used File and Store Necessary Fields #
#############################################################

# Open the sampling information file
my $isolatesUsed = shift @ARGV;
open USED, $isolatesUsed or die "Couldn't open file:$!";

# Initialise the necessary variables for parsing the file
my $line;
my @cols;

# Create a Hashtable to store the isolate metadata - removes duplicates
my %isolateInfo = ();
my $info;

# Begin reading the file
while(<USED>){
	$line = $_;
	chomp($line);
	
	# Skip the header line
	next if $line =~ /^IsolateID/;
	
	# Split the current line into its columns
	@cols = split /,/, $line;
	
	# Print out the columns into the output files if we don't want to ignore isolate
	$info =  "$cols[0]" . "\t\t\tMycobacterium bovis\t\t$cols[0]\t". $cols[2] . "\t\t". $cols[1] . "\tNew Zealand: " . $cols[3] . "\tCulture\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t";
	$isolateInfo{$cols[0]} = $info;
}

# Close the input file
close(USED);

#############################
# Print the Metadata needed #
#############################

# Open an output file
my $outfile = shift @ARGV;
open OUTPUT, ">$outfile" or die "Couldn't open file:$!";

print OUTPUT "# This is a submission template for batch deposit of \'Microbe; version 1.0\' samples to the NCBI BioSample database (http://www.ncbi.nlm.nih.gov/biosample/).\n";
print OUTPUT "# See SAMN02911891 for an example record of this type of BioSample (http://www.ncbi.nlm.nih.gov/biosample/SAMN02911891).\n";
print OUTPUT "# Fields with an asterisk (*) are mandatory. Your submission will fail if any mandatory fields are not completed. If information is unavailable for any mandatory field, please enter \'not collected\', \'not applicable\' or \'missing\' as appropriate.\n";
print OUTPUT "# All other fields are optional. Leave optional fields empty if no information is available.\n";
print OUTPUT "# You can add any number of custom fields to fully describe your BioSamples, simply include them in the table.\n";
print OUTPUT "# CAUTION: Be aware that Excel may automatically apply formatting to your data. In particular, take care with dates, incrementing autofills and special characters like / or -. Doublecheck that your text file is accurate before uploading to BioSample.\n";
print OUTPUT "# TO MAKE A SUBMISSION:\n";
print OUTPUT "#     1. Complete the template table (typically in Excel, or another spreadsheet application)\n";
print OUTPUT "#     2. Save the worksheet as a Text (Tab-delimited) file - use \'File, Save as, Save as type: Text (Tab-delimited)\'\n";
print OUTPUT "#     3. Upload the file on the \'Attributes\' tab of the BioSample Submission Portal at https://submit.ncbi.nlm.nih.gov/subs/biosample/.\n";
print OUTPUT "#     4. If you have any questions, please contact us at biosamplehelp\@ncbi.nlm.nih.gov.\n";
print OUTPUT "*sample_name	sample_title	bioproject_accession	*organism	strain	isolate	host	isolation_source	*collection_date	*geo_loc_name	*sample_type	altitude	biomaterial_provider	collected_by	culture_collection	depth	env_biome	genotype	host_tissue_sampled	identified_by	lab_host	lat_lon	mating_type	passage_history	samp_size	serotype	serovar	specimen_voucher	temp	description\n";

# Print out the Isolate metadata
my @keys = keys(%isolateInfo);

for(my $i = 0; $i < scalar(@keys); $i++){

	$line = $isolateInfo{$keys[$i]};
	print OUTPUT "$line\n";
}

# Close the output file
close(OUTPUT);