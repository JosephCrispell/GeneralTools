#########
# About #
#########

# Processing pipeline for Whole Genome Sequenced Mycobacterium bovis FASTQ files
# Created by Joseph Crispell
# Email: crispelljoseph@gmail.com
# This file created: 12-10-2017
# Last updated: 12-10-2017

################
# Requirements #
################

# Running on Ubuntu v16.04
# fastqc - tool to examine quality of FASTQ files
# cutadapt - used by trimGalore
# bwa - alignment algorithm
# samtools - package of tools to handle sequencing files
# blastn - tool to compare a sequence to a database of sequences
# bcftools - package of tools to handle sequencing files
# java jdk - programming language necessary to run some of the scripts
# perl - programming language necessary to run some of the scripts

###############################################
# Downloading and Indexing Reference Sequence #
###############################################

The recently updated M. bovis reference sequence is available here: https://www.ncbi.nlm.nih.gov/nuccore/LT708304
NOTE - the updated annotation file doesn't contain the repeat_regions or PE/PPE feature annotations
I have an updated version of this annotation file and am working towards replacing these annotations
The old sequence with working annotation file can be found here: https://www.ncbi.nlm.nih.gov/nuccore/BX248333.1?report=genbank

With the FASTA of the reference genome downloaded run the following command to index it:
$bwa index reference.fasta

###################################
# Instruction for running scripts #
###################################

## Stage 1 - Examine FASTQ file quality and define trimming parameters

Move to directory containing forward and reverse FASTQ files
Run fastqc using:
$fastqc *.fastq.gz

Examine each of the html files corresponding to each of the FASTQ files
Decide upon appropriate trimming parameters and enter them into a text file 
called PrinseqSettings_DATE.txt (see notes at bottom of page)

## Stage 2 - Align the FASTQ files against reference and create VCF files
Move to directory containing the forward and reverse FASTQ files and PrinseqSettings_DATE.txt file
Run ProcessRawReads_DATE.sh using:
$bash ProcessRawReads_DATE.sh [FASTQFileEnding] [pathToTrimGalore|false] [pathToPrinseq|false] [pathToReferenceSequence] [pathToPickRandomReads_DATE.pl] [pathToExamineBlastOutput_DATE.pl]

## Stage 3 - Merge the VCF files together
Move to directory containing the VCF files resulting from the previous step - "vcfFiles"
Run MergeVCFFiles_DATE.jar using:
$java -jar MergeVCFFiles_DATA.jar . [pathToAnnotationsOfRefenceGenome]

## Stage 4 - Examine genome coverage of isolates
Run ExamineGenomeCoverage_DATE.pl using:
$perl ExamineGenomeCoverage_DATE.pl [ReadDepthThreshold] genomeCoverage_DATE.txt

## Stage 5 - Filter the merged VCF file
Run FilterVariants_DATE.pl using:
$perl FilterVariants_DATE.pl 1 [ReadDepth] [HighQualityBaseDepth] [MappingQuality] [ProportionReadsSupportingAllele] [SiteCoverageAcrossIsolates] [GeneralQuality] [FQ] merged_DATE.txt

## Stage 6 - Rescuing site information for poorer quality isolates
Run RescueVariantPositionInfo_DATE.pl using:
$perl RescueVariantPositionInfo_DATE.pl [nIsolatesAlleleMustBePresentIn] [HighQualityBaseDepth] [ProportionReadsSupportingAllele] filtered_DATE.txt

## Stage 7 - Calculate isolate coverage in filtered OR filtered-rescued file
Run CalculateIsolateCoverageFromFiltered_DATE.pl using:
$perl CalculateIsolateCoverageFromFiltered_DATE.pl filtered_DATE.txt

## Stage 8 - Create FASTA file from filtered
Run CreateFastaFromFiltered_DATE.pl using:
$perl CreateFastaFromFiltered_DATE.pl 1 [MinNumberSitesBetweenVariantPositions] [pathToReferenceSequence] filtered-rescued_DATE.txt

#########
# NOTES #
#########

If you have any problems please drop me an email at: crispelljoseph@gmail.com

For additional help (including installation instructions of required tools) with particular
scripts use "-help" directly after script name

Full installation instructions for the require tools for ProcessRawReads_DATE.sh can be found
in the "checkToolsAreInstalled" function at the top of the script

PrinseqSettings_DATE.txt should look like this (change values where necessary):
LENGTH=50			# The minimum length of read to be accepted
MEANQUAL=20			# Filter sequence if mean quality score below x
TRIML=20			# Trim sequence at the 5' end by x positions
TRIMR=5				# Trim sequence at the 3' end by x positions
TRIMQUALL=20		# Trim sequence by quality score from the 5' end with this threshold
TRIMQUALR=20		# Trim sequence by quality score from the 3' end with this threshold score
TRIMTYPE="mean"		# Type of quality score calculation to use [min, mean, max, sum]
WINDSIZE=10			# The sliding window size used to calculate quality score by type
TRIMLTAIL=5			# Trim poly A/T > X length at 5' end
TRIMRTAIL=5			# Trim poly A/T > X length at 3' end
