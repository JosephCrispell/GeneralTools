#!/usr/bin/perl
# Searching files in a directory
use warnings;
use strict;
use List::Util qw(shuffle); # Shuffle function

# Author: Joseph Crisp
# Randomly shuffle the years associated with each isolate - in nexus and traits file

# Command Line Structure
# perl RandomlyShuffleYearsInXML.pl betweenYear input.xml output.xml

# Nexus File Structure:
#<?xml version="1.0" encoding="UTF-8" standalone="no"?><beast beautitemplate='Standard' beautistatus='' namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood" required="" version="2.4">
#
#	<data id="sequences" name="alignment">
#		<sequence taxon="TB1445_0.96">
#			NGCNGAGNNNNNNNNGNNNNNCTCCGCNGCGAGGANNGTGGAGCGGACNCNNNGGGAACNACGCCCGGCNTGNNNNNGGCGCACATCCNCTGCTCCCTTTANNNNNGCGAGCNNGGNACCGACTNNCCGACCATCTTCNTTGNCCTAGCCTGCNNCCGATCCCACNNCCGANCNCNCNNNCNCCCCANTCNAGTCCCGNCCGTCTGAGCNGNNCNNNNTTGGTNNNCTAGNACGCCNNCCNTNNNGTCNNNCCCGCCGCAGACATGCCGGNTTTGNTCANNCGANNCCNTNTCNNNNNNNNCCCGNCGGGATGCGCGGTGGCTGTAACATCAGGCCGNCGNCCCCNTGCTGCGCTCGAGCACAGTCCCGCGTAGCTCGCCGGCCCCGCGAGCTCTGGTNTCGAGGGCCCGNCNCGGNATNNGNGGTCTGTANNTGCGGAGGCANCCNNNNTGCTGNCGGGGGNACCGANNNNCNGTGNGGTGCCTAGGGCGGCNTGCGGGGGCGTGNNNNCGCCGGCGGCAGTCNGCGTCGCGGTTCGCCCCATCAGGCCGNNNANNGCGGCNNNNNCNTCNGNTGCGGGCGGCNNGNGCNACCCANGCGCCCCNCGGTNNCGCCCNNGCCANCCCCNCCGCNNGNANNCANNNNGNNNCTACGCGAGCTGCAGGCGATGTCGTGNGCCTCCGTGCGNNCGCGCAGCTGCAGGCCTCTNGNTCCGGCGGAGCNNCTNNNNNCAGGGCCCAGGCCTATGTACTNATCTGAGCCNNNCCGANCCATACCAGCTCTNGCNAGCNNANATTNNCGTGGCCCCCTCGAAGGGGCGCGACCAGGACCGATACCCCACCCGTCGTATGCNTCNGCCGCNCNGATGAGCTGNGNNNNGGCTCCGCNNNNNNNNATGTTTNGCACGTGCCCGCGGTGAGGTGGNGNAGACCNNGGTGCNNGACNCNNNNAGNCTGGACACATGACGACTTCNNGGGNTCCCTCTGGTCGGTGGGNCCCCTTCNGCGTAGGCCCCGTCGNGTCGTCGGTGAGGCGTGATGCNGGGCNGNGTGCCGNNNCGCAAGGNCGCTCGCCTCNGGCGGTNCGGCNCGTGCATGNNGTTGTGGTTTACTGTATGCNATCTAGCATGCGGGAGGNAGCGTCNCGCGTAGNGCATAACGGCCGCGCGGTNANGGGCNGAGTCTCCCATGACCTCNGGTCGACNGCNNNCTNACNNNGTGNNGCAGGCGGCCCATCCGGGNGAGTCCGTNAANNNNACCCCTCACTCTGTCTACCGANNGCTCCGGGNCNAGCTGGGGGACNGNNNNNCGCTGCGCCCGGACANCCGGNNNCGATTCNGGCCNNNACAACCCCCGGGNCNCTAGCCCAANGNGTCCCCGCACTGCCGGNCCCCGGTNGGNGGAGGCCCNNNGNCGGCCTTCTGCGCCCGGGGGATGTGNCCGCACCAAGGGAACCCGCCNNGCGCCGCACGGAAAGGGGGGGGCGCCGCTCTGNTGNGCGNCAGATACCTAGCGCNCCCNGNANNAGGCCCCANNNCTCCGACATCNCGCCGTCCNCGTTCGTGTCGGGGGNGCCGGACAGGATGCCGTGNNNACCCACGCGGNCCNNGGNGNNNGCAANNGANNNNNCCNNNGCCCCGTAGGACGCGCCATGANAGNNCCAAATCTNGGAAATNCTNNNGCCGCCCAGCTGGACANACGTCAGNNCNNNNNGGGCCAGCTGGAAANNCGGGCGGNGANNNNGCGCCAGCCGTACGCNCTCCTNTGGAGCTCGGATCTNCCCCAGTCNNNNNTGCCCNNNGATGNNNNCCANCTGTGGCGCGACGTCGGANGCNNCCNNNGCTATACACACCCNTCGACGNACGGCATCGGGGATCGCGNNGNTCTGCTGCNCATTNTCCAGCGCTCTGTNNCGCGGANATNCNNCCNNCGGGTCCCNNTCCCCGGGGCCNNNCNGNCCGGGCNGCGCCCTGCACGGCAGCCNGNNGGGGCANCCAGGAGGATCCACCNCACGNGACCCCGCCNNCGCANCCANNGCCCCATCACGCGCGGAGGGGACCGGTCTNCCNCCTTTCNNTNNACGCGGNNNNNCNNNNCCCGTCNGGATAGTTCAGACAGTGCGTAGGGACANGGGCTNNNNAGGCGCTGTCGCNNTCTNNGGGCAGCCTTNNGCTNGGGCGNTGGNANCCCTTGCGGCTCACGNTCGCNAGGAACCGCAGCAACCAGNNGGGAGGTNCGCATGCGCAACNNATTCGGCCCNNCTTGCCTCGACGNAACNGTNCGNACTCNCGTCAGCGCNNGCCNNNCNNNCTNNCCCCCANNCCCCCANCACGGGGCGGACGCGCACTNNNCCNNCNAGANNNNNNACCCCGGNTGGNNNNGNNNCGGTAGGGTCAGCGCACCCTGACCCGCCGTGNNGCNNNNNNNGTAGCTCNNCGTCGCNNTCCTCGTCTCTGTGGGGCAGGCTANNGNTCTGGCGGCATCGNNNTGNACGTTGGGCNCCTACCCGCAAACGCCTGAGCGANCTCCCACTCTGCCCNNNGNCCCCGCGGGTGCGCCAACCCTTCGNCCNAGNNNNCTNNGCCNGGGGCANNTTCANGCCNTNNGGTGGTCGGCNTCNNNCNCCTTAGCGTGNNCTCTCNAGAGGACCCGCANCTGGGGGAGNNCGAGCTCNCNGGNCCANNGGGNTNCAAANNCAGNCTNCAGNNNGNNAGGCGGTTCGCTNACCNNNTCTTCCACATGNCCCCGGAGGCGTCAGTGGCTCTGAACNANNCTCANNACGCCCCGNCGCCATGTACGCTTANCCCAGACNCCNGGCGGATCGCGTGNGGNTTGCCGTGGCNCGGNGNNTGGNAAGTGNNGCGAGNCANNNGGNNGTGCNNCNGTCAGCTTGGTAGCACGGCCTGNTANACTTGCAGCACNGGACAGCCGAACCCCCCGNGCCGCACCGGCGCGTNNGCATGGGCGGGCCGAATCGGNNCNNTCANNGGNNCGCACNNGGCTTNCATACCNCTGNCGTACCCGAGGANTGCGGANTCNGNCCGGGANCGGCNGCGCGCCCTNGCNGGNCNNNNNNCGGNNGNGCCAGGGCCCNGCNNTGTANANNNNGCNNNCCANTCANCCNCNCGGCGCGTGNAACCGACNNNNNGNTNNNNAGGANNNGCNNGNCNCCTATANTCAGCNNNGGNNNNNNNNNCTNNNNNNNTNNGGNNCNNNNCCTANNNNATTANCCGGCATTCGATNCCGCGCCCGCTGCNGCGCNACAGCGTCATNNGNGGGTAACANNGNNNACGCTGGCCCCCCCNNGGGCCNCNNNGGTGGCGNNNNNGCACGTGGGNNGCGNCGCAGCAGGCCGANNCGNATCTGACGNTGNTGCCACCGCGNGCTGCTCCACCGNGGNNGCCCGCNNNCCTACGAGCGATCCNNNCNNCGATTTGCCNTNGCCACTACCGGANGAGGTCGGTCNCCCGGTNCANACAANGNTCGCCCACNCGCGCANTNTCGNTCCTTCGNNGCCGCCATNCGACNGAACACCTCGTATCCNCCGGGNANCGACGGCGGCCCCTTTACAGNTCGGCGCGTTCNNNCCNCTCCGGGGAGCGGNCCCCAGGCCAGNNCCTCNNNAAGGCTNCGTGCCNCATCTTNGANAGCGNNGNCTTGTNGGCNGNATCGCAGTGTCCCCNGAGGGGAANNGCCGNNCNCNNGCCTCTGAGTTGGGGCGTGCGNTTTAGNNNNNGNNNGGAGTANCGGACTGTTAANGGGCCCNNGGNCNNANCTCACGCNNTTTGCANGTGNCCTANNNNTCGCGAANNANAGGCNNNTGCTCCCGATCGCCTCTTCTCCNGGCCTGGACTACGGNNGTNTGNCACAGGGCGGNGGGGCGGGTGGCCTCNCGAAGCGCACCTNCNGANGTGCAGCAACGAGCGCNCTNANCTAGACGGNCTGACCNGAGNNCANGCTNANNCCTTCGGCGCTCTAAGGNTACCCTNCCGGCGGCGGTTGGCGCCCGCCNNNGGGNGANNNNNNNNNNTNNNNGNNGGCTCNNNNNNCCCCGAGNNNNATGGCNCTGNATNACGCACTCCGTACGAGGAGGCNCCNCGCCNNACGACCGCTCACGGGGTGGGGNGNNGCCGCNNNGATAGCNTCNNNNGGGGGACCCCATCGGGTCGCCCCTCCCCGAGCGCGCGTTCTGCCNGNCNATGGGGATCGACNCGGCNCGGTCCGAGGGGNTNCNNCGGGNGCNNTGGGTTTGGTACTCCNACAGTNCTGNNNNNNNNCCCNNNNNNNNNNNNNNTTCTAANNNGATNNGGCGGGTNNGGAGAGCGGGNNNGCTNNNNTTGCNNNNGGNNCTGGCGCCCCGNNNACCGGGGGNCTCTTNNACCCGCCGCNNCTCCCGTTTGCANGCCTTTCGANAGGACACTCCAGACTGGCCCNNNANTCAGCNACNNNGGGGNNNGGGTCCNGGNGTCNAGCGGGGGNTGGGGGCGTCCAAGTACANCGAATNNNATNAGCGTGGCCGGGCNNNAANNAGCGNTNTNNGNCTCNANNNTNCGGNGGAGNNGGNCAGCCGGNAGGGGCAGNGGGNNNNACCNNNTGCTAGCNGGCAACCCTCAGCTNGGGTATGCGNNNNGTNGGGACCTGGNCNNNTCCCGTNACCCGANANNNAGNNNNNNNNACAGTANNNNNNNNNCTGCCNGTGACTNGTCGGNACGTTCGNAAGGGTNNNNNNANGTGNGCNNNCNNNCACGNNNGANNNNTCGGNGNNNNNGCGCCATCCGCGGNNCCCCTGANNTCAGCTGCTGNAGNNAGGGTAGAGCCCCGCNGGGAGNGNCTGCNNCNCNNNNCNGNGGTTTGTNNNNACNCGCGGAGNNCCTNGGGANCCCTGGGTCNNTCATANTNANNNATGNGGTCAGGCCCNTTGGACGGAGAGNGTGNNACNGGNTNTCCGCAGGCCCACGCCCACGCCGACAGGGCGCACGCCGGTCCGCTANGCNCGCGGCGNNNNNCNNNNGNNTCACCNNNNNNNNNGACTGCCATCGCGAATTAANNGGCAGNNCACTCTCAGCCCTGNCNNNTACTCGCTGGNGGNCNTNNNCAGTGGGCGACNNNCCCTCCCCGNGCNGGCCCCGCGCTTCGGCTANCACCGAGGCCCCANCCGGTATGTNNAGNCCACCCGCNNCGAGCCTGGCCGTCNAAGTCCCGGGCGCGGCNACGTAGGGCNNNNNCGGCTGGGAAATCNNGGNNAGAGNGGNNNNNCCCNTGTNGCTCTGGCGCGNNNCGCNNCNCACNCNNTNNGNNNNNGNGATAACCCNNCANNCGCGCCGTACCGGGGTGCCCNGGCCTGGCAAGNGGAGCCATGNCCGATCCGCCACNGGAGCACCTGGAGCAGAGCTCNNCGCATAGAACACCCGCNNCAGNGAGCGATGCCGCNACAGCGCAGGGGCCGTNGCCGNTGGGGGCGTTCGNGCCCCNNNNAGCGNACNCCCGCNCGGNNNAGNNNNACTCGGTCGCCGCCNCCAGGCGATTCGCCGTGATGGTCNCNCGGCGAGCNTGNANGGGTGCAGNNTAGGNGCGCNNCGGNNNCNCNGCGGAACGTCTGGNCCCCCGCANNGCNNGGNAGGAGGGNNTGGGACNCAGGCGTCTGCGGGGTCCNCCGCNGGCATNGGTCCCTGTCCATGCGGCGNNGCGCCTNNCNCNNNCGGCCTAGNTACCACANNNNNGGNCACGTCNCCNCCCCAGACNGTNNTCGTTGGGCGCGGGGGCTCAGCCATCGGCAGGTANNCCCCCTACGGTCGNCCATNNNGATGCCAGTCAGGCANCNCCANACNGNGTCGCGGGGCGGCGGGCGGNGNCGGCGACCNTGTNACGGATAGGNNGNCCNNGNCACGNCCNNGGNGCTGANCCNCAGTCATCGGCCGCTAACGCTCGNNCNNGNCNNNNTNATNCCGGNNTAGNCGGCGNCNANANCCGTANACGTNNCGNTNGCCNGCCGATTGGCCGGNTGNNCCGNGCGTGGTTGAGAGGGANCAGGACCNAGCCCGACNTGNNGNNTCANNNACAGNACGGTTNCGCTGNCANNNNNCNAGGGTGGCAGGGCTCGACTCCNCGAGGGNCANGGCTGTCGAGACCCCTCGNGAGCCCGGACACTGTCNCTNAGNCCTCCCNGACGCNAGCNNNATNCGAGNNNCCGTACTNCTCNTCCGGCCGGGCGTTGTTNCCCCGTTAGANNNNGTGGACGGAGGNCNNGACCNANCCAACCCTAGGNNNCGTGCGGTCGGNGACCGGGTGCGCATGTATGCNGCCCGACGGNNTGGGGGGGGGACCCTATNGTGGCNNGNNGTGCGCGACACNNCTATGAACNCACGCGAGACNNCNTCCCTTGCCGGCGTGGTGCCCCGACACGCNCGTCANCGCGGCCGCCGNGTACGCCCGTGCGCTTTTGCNAGNANGGNNCNATNACCGGGTCACCGGCCNNTGGGGNNNNATCTGGAGNNNGNNNNNANNNNGGGCCGNGATTTANTCGTGAATCGTAGCTCTGCCAGGNCACCANTCTCGNTGTCCNGGAGGCCCATGGCGGNGATTNTCTGGTCGGTATCGTTNNCCAANANTCACNNNACGCCCNTGCCCTTGNNCAGCCAANNNCGAGGGGTCCGCCAGAGCCCTTCNCTCGGGNNCNGGCTCTCCACGCGGGTNCGCNGTATNACCGNTTTNGNGGCGCNGGCGCGGTGCNGNCGNGCTTNTGGTCCGCTGGCACCNNNCATCAGCGGCCGGCGNCCGGTGGGGGANACGNCCGGGGCNNTTGCTCGCGACAGGCNCCCACCAGCTCGTCANTGTGNTNCANGGANCGANNTNCCGGTTCTCTNACANNNNANNCNGCCGCGGGGGGGCAGGANGNGCCGGTGTGAGGGCCNGCNCGGNCCCNNCCGNGCCGAGTCGTTTACGGNCNGGAACNNCGTGCCCACCGGCGGGCCGNGNGNNTGGCGGTNGTNACGNNNGTGCGANCANGCATNNNNCGACGTTTACNNGATCACGCGNNGNNGTCGGCCCCCTGCTGNNCANNGNGTCNACNTGCACNGACCCAAGNCCCACTCNGNTGCCNNAATNTCTGGGCACGCNNNCNCGGTNNACACGGGTGCCCNTGGANACCTCCGGNNACNNCTNATCGNCNTTACACCGGTTGGCGNGGGNNNNTNTGGCGAACNAGGGCCTGANNGNNTTGCGCNGNGNGCNGCCGCCGGGGNNAGGTCCCCGTGGTGGCNGATGGNCNCGAGANTGAGTCGNNGCCCCTNNCGCGGCGCATCGGCCCGTCACTACNCCCCCTNCGANNCTGNNGATAACGGGCCTATTCTTNCNTTTGCACCNNCANNNNCCCGNTNNNGGGTTNGTNNNNNNCNNCGNNGGCCCANGGTCNCCNGNTGTCGNCGTTGGCGACCACNGGGCCANNGAGGGGTCGGCGCCACCGGCNCGCCCTCGAAANGTGTCTCTGGCNCNNNGCGCCNTGGCCGCNGTGTGGNNTGGGNGCCGCGCTAAATNNGGNCGNNNGGAANNNGGGAANNCGCACCCNCACNCGGCAGCNGGCGNNNNGGCCNGGCGNTGCCTCNANGCCAGGAGGCTAAAGGCCATACGCCCGANCGNNNNNGNNNNNNTGGNNNCNNNNNNNCNCTAGNCACTGTCNCTCTANCGNNNNNCTGCCGNACCCTGACNTATATCTAGAGGGAGNCCGNGNCCGNGGTCGTCTCCGCGANCANTGGTNGCGNNCTNGNGGGACTGAGCCGGCGGGTCGCCCGTGTCNNNNGTTAGTGNCCCCNGCCGGCGCATGCNGGTCNGTGATTGCTGTCCNCCTNCGCNGACTCCGCGCCGATCCGCGGCATNNNGATCTGGCCGCGACTACACCGTGCNCCNNNTGNCCGTNCGGNGCGCTNNANGGACCCGCCAGCCGNANGTGGGCCCGNGCCACGNNNNNANNGCGGAGCNGNCCNTCCCANGCTGGNNCCCGGGTCGCTGCGATNCNGCGCGCGNTNATCGGCNACGNCGNCAGGCANTTCGCANCGGNGAGGCGACCCCTCGTCCCGCACNNNGNNGTAGGCNNNGNAAGCCCTGTTGCCNGTGNGCCNTANNGGGGGGGAACGCTNCGGTTGNTTTGCTGGCCAATGNGGGCCCGCCTAACCGTTTCGTCCACCTCTCNCGCCGGNCGTGGAGCCNAGTNCCGGTCGAGCGGATCCTGGCTGTNTTTGCGANNNAGNNGTANNNGCGCCCCCCTGNNNAGGTCAGCTGNCACGCCTTNGCTCGCTCGGGGACGCCGTTAACCGGCCTTGCCGGAACCNCGGGCCCCATATCGCAGCCTACGGCCTNGTACGNNCNTGTGTAGGGGNNCCGCAGCGACCGGAGGCNNTGCNNNNNNNGNCTTGNTNNTGCCCGCGTCGGGAGGGATTCAATCCCCACGGCNNNAGGAAACAGCCCCGCCAGGGAGTTAACCCNCGGGGGCCGANANATTCTGGAGCCGGCGTTNNCGATNCGCATCNNCCGAGCCGGCNCCAGGGCGGGTTGGGTNCAGTCTGGGGTCCGGGCCGAGCCGTCGCGCCCTATCNNNCGAGGCGCGGGACAGCTGGCTGGGCGGCGNGGGCNTACGCTGAACGCCGGGCGTGCTGGGGTGGGGNGGATGCNGAGNACGNGCCNNCGCGNATGCCNCNGGTCTGCNCCGTCNAGCGCCAGCGGCGGTCGCTGCNNNNNNNNNGCNGGGGTGTACGGGTGGGNGCGTCCGGGNGGAGGTCAGGCCGGNCGNTNTNNGNNNCNGGNCGTCACNAGCCNAGCNNNACTGAGGNGNNNNNNGTCNTTCCGACGAAGCCCGCNAANACCTCAGGTGCACTCGCCGGCNNGCCATNGANNGNNCCNGAGATCTGNCCGGGANANTAGCAGGCGGGNNNNNCGCNTNCNNCGCCCCNNCNNCGNCAGCCGCGNNNCACGGGGCNANGGGCATATNGGCTCNCACGGCNCCGATGACGACTCCNGGGGTGAGTNCACCCNGAAAACAAGTGTTCGANGGCCCCTCCTCGGGGACTGCGCCCTGTGGGTGCCCNCNGGCCANAAGGCGCCGCGCGGGGCACGNNTGTTNGNCATGNNCGCGNGGGGACCTCAGNTGCTTGGCNANTGCNGGTANGACGGTCCACGGAGCTCNCNNCNNNGGNCCTGCAACANNNNANCAATNCACCCACGTNNNGTGCTGTGACGATCNNNNNNCGNANGCGTCAGGGGGTCNNGCCNGNATCNNGCAGGCCGCANNCAGCCCCNNGCCGNCAGATTGNGTGNNNNCTNNNTCANANTATCCANNNNNNNTAGATCCAGTNCGAGGGGGGTGCTTACCAACTNNCNNNNNTCNGATGTCGNCCGACNCGGCCTGGNNGNGCCCGGCCGACGGCGCCCGAAGGCCCGNNNCCCGGGAGTGNCCCGCGNNNNGCGGGGGGNNNN
#		</sequence>
#		...
#	</data>
#    
#	...
#
#	<run id="mcmc" spec="MCMC" chainLength="100000000" storeEvery="100000">
#		<state id="state" storeEvery="5000">
#			<tree id="Tree.t:sequences" name="stateNode">
#				<trait id="dateTrait.t:sequences" spec="beast.evolution.tree.TraitSet" traitname="date">
#					TB1445_0.96=2008.85792349727,
#					...
#					TB1390_0.65=2005.18082191781                
#					
#					<taxa id="TaxonSet.sequences" spec="TaxonSet">
#				...
#		</state>
#		
#		...
#
#		<logger id="tracelog" fileName="TipRandomisation_1.log" logEvery="10000" model="@posterior" sanitiseHeaders="true" sort="smart">
#			...
#		</logger>
#
#		...
#
#		<logger id="treelog.t:sequences" fileName="TipRandomisation_1.trees" logEvery="10000" mode="tree">
#			<log id="TreeWithMetaDataLogger.t:sequences" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@Tree.t:sequences"/>
#		</logger>
#
#	</run>
#
#</beast>

## Note the option for random shuffling
my $betweenYear = shift @ARGV; # Do you want shuffling to only be done between years? [0,1]

## Open XML File
my $xmlFile = shift @ARGV;
open XML, $xmlFile or die "Couldn't open $xmlFile:$!";

# Initialise two arrays - storing the isolate IDs and sampling times
my @isolates = ();
my @samplingDates = ();
my $isolateIndex = -1;
my $foundDatesSection = 0;

# Initialise an array to store all of the XML file lines
my @fileLines = ();
my $lineNo = 0;

# Create necessary variables for parsing the file
my $line;
my @parts;
my $part;

# Begin reading the Nexus File
while(<XML>){
	$line = $_;
	chomp($line);
	
	# Remove any special end of line characters
	$line =~ s/\v//g;
	
	# Store each of the file lines
	$lineNo++;
	$fileLines[$lineNo - 1] = $line;
	
	# Search for the date assignment section
	if($line =~ /beast.evolution.tree.TraitSet/){
		$foundDatesSection = 1;
	}elsif($foundDatesSection == 1 && $line !~ /TaxonSet.sequences/ && $line =~ /=/){
	
		# Get the isolate ID and sampling date from the current line
		$part = (split /\t/, $line)[-1];
		
		# Remove comma if present
		$part = substr($part, 0, -2) if $part =~ /,/;
		
		# Get the isolate name and sampling date
		@parts = split /=/, $part;
		
		# Store the isolate name and sampling date
		$isolateIndex++;
		$isolates[$isolateIndex] = $parts[0];
		$samplingDates[$isolateIndex] = $parts[1];
	
	}elsif($line =~ /TaxonSet.sequences/){
		$foundDatesSection = 0;
	}
}

# Close the input XML file
close(XML);

## Shuffle the sampling dates

# Initialise an array to store the new random sampling dates
my @shuffledSamplingDates;

# Check if shuffling is to be done between years
if($betweenYear == 1){
	# Initialise an array to store the randomly chosen indices for the dates
	my @randomIndices = ();

	# Initialise variables to store information about random selections
	my $randomIndex;
	my $year;
	my $yearOfRandomChoice;

	# Randomly choose a date for each isolate that doesn't fall in same year
	for(my $i = 0; $i < scalar(@samplingDates); $i++){
	
		# Get the sampling year of the current year
		$year = (split /\./, $samplingDates[$i])[0];
		$yearOfRandomChoice = $year;
	
		while($year == $yearOfRandomChoice){
	
			# Chose a random index
			$randomIndex = rand @samplingDates;
		
			# Get the year of the randomly chosen isolate
			$yearOfRandomChoice = (split /\./, $samplingDates[$randomIndex])[0];	
		}
	
		# Store the random index
		$randomIndices[$i] = $randomIndex;
	}
	
	# Store the randomly chosen dates
	@shuffledSamplingDates = @samplingDates[@randomIndices];
}else{
	@shuffledSamplingDates = shuffle @samplingDates;
}

## Print out the XML file with shuffled dates

# Reset the variable noting whether the dates section has been found
$foundDatesSection = 0;

# Open the output file
my $outputXML = shift @ARGV;
open OUTPUT, ">$outputXML" or die "Couldn't open $outputXML:$!";

# Get the input and output file names minus suffixes (and remove path from input)
my $inputPrefix = (split /\//, substr($xmlFile, 0, -4))[-1];
my $outputPrefix = substr($outputXML, 0, -4);

# Print each of file lines from the original file - skip smapling dates block
foreach $line (@fileLines){

	# Replace input file name with output
	$line =~ s/$inputPrefix/$outputPrefix/g;

	# Print the XML file line out to file as long as it isn't part of the sequence block
	print OUTPUT "$line\n" if $foundDatesSection == 0;

	# Check whether we have found the sampling dates block
	if($line =~ /beast.evolution.tree.TraitSet/){
		$foundDatesSection = 1;
		
		# Print out the sampling dates
		for($isolateIndex = 0; $isolateIndex < scalar(@isolates) - 1; $isolateIndex++){
			print OUTPUT "\t\t\t\t\t$isolates[$isolateIndex]=$shuffledSamplingDates[$isolateIndex],\n";
		}
		print OUTPUT "\t\t\t\t\t$isolates[-1]=$shuffledSamplingDates[-1]\n";
		print OUTPUT "\t\t\t\t\t\n";
	}
	
	# Check if we have reached the end of the sampling dates block
	if($foundDatesSection == 1 && $line =~ /TaxonSet.sequences/){
		$foundDatesSection = 0;
		print OUTPUT "$line\n";
	}
}

# Close the output file
close(OUTPUT)
