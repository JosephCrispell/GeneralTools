#!/usr/bin/perl
use warnings;
use strict;

# Author: Joseph Crisp
# Conducting and Parsing Google Scholar searches

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
	}else{
		print "No article containers found in page!\n";
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

#########################
# Set search parameters #
#########################

# Set the user - needed to access google scholar
my $user = "JosephCrispell";

# Set the search criteria
my $all = "";
my $exactPhrase = "Whole Genome Sequencing";
my $atLeastOne = "tuberculosis mycobacterium";
my $without = "";
my $where = "title"; # "any" or "title"
my $minYear = 2000;
my $maxYear = 2017;
my $start = 0;

##############################
# Get google scholar results #
##############################

# Loop through each page - until no results found
my @starts = range(0, 10, 10);

# Create arrays to store the overall list of titles and a list for each google scholar page
my @currentPageTitles = ();
my @titles = ();

foreach my $start(@starts){
	
	# Build the URL for the google scholar webpage
	my $url = "\"https://scholar.google.co.uk/scholar?as_viz=1&as_q=$all&as_epq=$exactPhrase&as_oq=$atLeastOne&as_eq=$without&as_occt=$where";
	$url .= "&as_sauthors=&as_publication=&as_ylo=$minYear&as_yhi=$maxYear&start=$start&btnG=&hl=en&&as_std=1,5%2C5\"";
	
	# Get the html formatted text from the URL
	my $page = `wget -U $user $url -O -`;
	
	# Get the article containers
	my @containers = getArticleContainers($page);

	# Check that there were results on the current page
	if(scalar(@containers) == 0){
		last;
	}
	
	# Get the titles from the article containers
	my @currentPageTitles = getArticleTitles(\@containers);
	
	# Combine the titles from the current page with the growing list for all the google scholar pages
	@titles = (@titles, @currentPageTitles);
}

for(my $i = 0; $i < scalar(@titles); $i++){
	print "$titles[$i]\n\n";
}
