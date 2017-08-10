#!/usr/bin/perl
use warnings;
use strict;

# Author: Joseph Crisp
# Conducting and Parsing Google Scholar searches

# Searching google scholar:
# https://scholar.google.co.uk/scholar?as_q=&as_epq=&as_oq=&as_eq=&as_occt=any&as_sauthors=&as_publication=&as_ylo=&as_yhi=&btnG=&hl=en&as_sdt=0%2C5
# url can be broken up into:
# 	https://scholar.google.co.uk/scholar?
#	as_q=					 	# All of words
#	&as_epq=				 	# Exact phrase
#	&as_oq=						# At least one of
#	&as_eq=						# Without
#	&as_occt=					# any=anywhere in article title=search in title
#	&as_sauthors=				# Authors
#	&as_publication=			# Journal
#	&as_ylo=					# Min Year
#	&as_yhi=					# Max Year
#	&btnG=						# ?
#	&hl=en						# Language
#	&as_sdt=0%2C5				# ?

# Notes:
# The href immediately preceeding: "gs_ico gs_ico_nav_next" will take you to the next page as long as "aria-label="Next" disabled" isn't present
# I think "class="gs_ri"><h3 class="gs_rt"><span class="gs_ctu"><span class="gs_ct1">" means the 
# "&start=10" can be added onto the url to move to next 10 - e.g. next page - google limits to 10 entries per page

#########################
# Set search parameters #
#########################

my $user = "JosephCrispell";

my $all = "all+these+words";
my $exactPhrase = "";
my $atLeastOne = "";
my $without = "";
my $where = "title"; # "any" or "title"
my $minYear = 2000;
my $maxYear = 2017;

my $url = "\"https://scholar.google.co.uk/scholar?as_q=$all&as_epq=$exactPhrase&as_oq=$atLeastOne&as_eq=$without&as_occt=$where";
$url .= "&as_sauthors=&as_publication=&as_ylo=$minYear&as_yhi=$maxYear&btnG=&hl=en&as_sdt=0%2C5\"";

##############################
# Get google scholar results #
##############################

my $page = `wget -U $user $url -O -`;

print "$page\n";

