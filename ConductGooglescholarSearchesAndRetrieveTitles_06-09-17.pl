#!/usr/bin/perl
use warnings;
use strict;

# Author: Joseph Crisp
# Conducting a Google Scholar search and retrieving the titles of each of the results

# Command line structure:
# perl ConductGoogleScholarSearchesAndRetieveTitles.pl output.txt

# Searching google scholar:
# https://scholar.google.co.uk/scholar?as_q=&as_epq=&as_oq=&as_eq=&as_occt=any&as_sauthors=&as_publication=&as_ylo=&as_yhi=&btnG=&hl=en&as_sdt=0%2C5
# url can be broken up into:
# 	https://scholar.google.co.uk/scholar?
#	&as_vis=1					# Turn off patents and citations - in combination with &as_std=1,5 at end
#	&as_q=					 	# All of words
#	&as_epq=				 	# Exact phrase
#	&as_oq=						# At least one of
#	&as_eq=						# Without
#	&as_occt=					# any=anywhere in article title=search in title
#	&as_sauthors=				# Authors
#	&as_publication=			# Journal
#	&as_ylo=					# Min Year
#	&as_yhi=					# Max Year
#	&start=						# Result number to start with - google presents 10 results per page - CHANGE TO GWET TO NEXT SET
#	&btnG=						# ?
#	&hl=en						# Language
#	&as_sdt=1,5					# Used in combination with &as_viz to turn off patents and citations

# Notes:
# <h3 class="gs_rt"> is the class for each result. Title is found between >...</a>
# Think about changing ip - system('ipconfig /release & ipconfig /renew');
# Think about changing user or cycling through random users

###########################
# Check if help requested #
###########################

my $firstArg = shift @ARGV;
if($firstArg eq "-help" || $firstArg eq "" || $firstArg eq "-h"){
	print color("green"), "Perl Script to conduct a Google Scholar search and retrieve the titles of each of the results:\n", color("reset");
	print "\tperl ConductGoogleScholarSearchesAndRetieveTitles.pl output.txt\n";
	print "\t\toutput.txt\t\tPath to output file where the result titles are stored\n";
	print "\nPlease note that search parameters are currently specified in script. Sorry!\n";

	exit 0;	
}

#############
# FUNCTIONS #
#############

sub getArticleContainers{

	# Get the page text from the input argument
	my $pageText = $_[0];

	# Intialise an array to store the article container text
	my @containers = ();

	# Check if there are any present
	if($pageText =~ /h3 class="gs_rt"/){

		# Split the text into characters
		my @pageCharacters = split //, $pageText;
	
		# Initialise variables to parse the webpage
		my $subString;
		my $containerStart = -1;
		my $index = -1;
	
		# Examine each character - search for the <h3 class="gs_rt"> ... </h3> container
		for(my $i = 0; $i < scalar(@pageCharacters); $i++){
			
			# If not yet found container start - <h3 class="gs_rt"> - look for it
			if($containerStart == -1){
			
				# Build a substring
				$subString = substr $pageText, $i, 18;
			
				# Check if matches start of container
				if($subString eq "<h3 class=\"gs_rt\">"){
				
					# Record the start of the article container
					$containerStart = $i;
				}	
			
			# If found container look for end - </h3>
			}else{
			
				# Build a substring
				$subString = substr $pageText, $i, 5;
			
				# Check if matches end of container
				if($subString eq "</h3>"){
				
					# Store the container substring
					$index++;
					$containers[$index] = substr $pageText, $containerStart, ($i + 5) - $containerStart;
					
					# Reset the container start
					$containerStart = -1;
				}	
			}
		}
	}

	return @containers;
}

sub getArticleTitles{
	
	# Get the containers from the input arguments
	my @containers = @{$_[0]};
	
	# Initialise an array to store the article titles
	my @titles;
	my $index = -1;
	my $title;

	# Examine each article container
	for(my $i = 0; $i < scalar(@containers); $i++){
	
		# Check if hyperlink is present and retrieve title
		if($containers[$i] =~ /<\/a>/){
		
			# Split the container text and retrieve title
			$title = (split /">/, $containers[$i])[-1];
		
			# Remove the characters at the end of the title
			$title = substr $title, 0, -9;
		}else{
			# Split the container text and retrieve title
			$title = (split /<\/span>/, $containers[$i])[-1];
		
			# Remove the characters at the end of the title
			$title = substr $title, 0, -5;
		}
	
		# Remove bold indicators
		$title =~ s/<b>//g;
		$title =~ s/<\/b>//g;
	
		# Remove formatting indicators of the form: "&xxx;" - curly brackets allows 1-10 characters to be between "&" and ";"
		$title =~ s/\&.{1,10};//g;
	
		# Remove spaces at start
		while($title =~ /^ /){
			$title = substr $title, 1;
		}
	
		# and end if present
		while($title =~ / $/){
			$title = substr $title, 0, -1;
		}
	
		# Store the title
		$index++;
		$titles[$index] = $title;
	}

	return @titles;
}

sub range{
	
	# Get the input arguments
	my $start = $_[0];
	my $end = $_[1];
	my $by = $_[2];
	
	# Initialise an array to store the output
	my @array = ();
	my $index = -1;
	
	# Fill the array
	for(my $value = $start; $value <= $end; $value += $by){
		$index++;
		$array[$index] = $value;
	}
	
	return @array;
}

sub printTitlesIntoOutputFile{

	# Get the output file name from the input arguments
	my $output = $_[0];
	open OUTPUT, ">$output" or die "Couldn't open $output:$!";
	
	# Get the titles from the input arguments
	my @titles = @{$_[1]};

	# Print the titles into the file
	foreach my $title (@titles){
		print OUTPUT "$title\n";
	}

	# Close the output file
	close(OUTPUT);
}

#########################
# Set search parameters #
#########################

# Set the user - needed to access google scholar
my $user = "Joseph Crispell"

# Set the search criteria
my $all = "";
my $exactPhrase = "Whole Genome Sequencing";
my $atLeastOne = "bacteria bacterium";
my $without = "";
my $where = "any"; # "any" or "title"
my $minYear = 2007;
my $maxYear = 2017;
my $start = 0;

##############################
# Get google scholar results # DOESN'T WORK FOR LARGE NUMBERS OF REQUESTS - LEADS TO BLOCKING OF IP???
##############################

# Loop through each page - until no results found
my @starts = range(0, 1000, 10);

# Create arrays to store the overall list of titles and a list for each google scholar page
my @currentPageTitles = ();
my @titles = ();

# Initialise variables for parsing the search results
my $url;
my $page;
my @containers;

# Intialise variables for sleeping
my $secondsToSleep;

foreach my $start(@starts){
	
	# Wait random number of seconds - between 10 and 60 seconds - trying to not get stopped by google
	$secondsToSleep = int(rand(51)) + 10;
	print "Sleeping $secondsToSleep seconds...\n";
	sleep($secondsToSleep);
	
	# Build the URL for the google scholar webpage
	$url = "\"https://scholar.google.co.uk/scholar?as_viz=1&as_q=$all&as_epq=$exactPhrase&as_oq=$atLeastOne&as_eq=$without&as_occt=$where";
	$url .= "&as_sauthors=&as_publication=&as_ylo=$minYear&as_yhi=$maxYear&start=$start&btnG=&hl=en&&as_std=1,5%2C5\"";
	
	# Get the html formatted text from the URL
	$page = `wget -U $user $url -O -`;
	
	# Get the article containers
	@containers = getArticleContainers($page);

	# Check that there were results on the current page
	if(scalar(@containers) == 0){
		print "Found the last page of results! Started at $start result\n";
		last;
	}
	
	# Get the titles from the article containers
	@currentPageTitles = getArticleTitles(\@containers);
	
	# Combine the titles from the current page with the growing list for all the google scholar pages
	@titles = (@titles, @currentPageTitles);
}

########################################
# Print the titles into an output file #
########################################

# Open the output file
my $output = $firstArg;
printTitlesIntoOutputFile($output, \@titles);