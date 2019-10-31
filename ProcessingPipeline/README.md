# Processing pipeline for Whole Genome Sequenced Mycobacterium bovis FASTQ files
## Created by Joseph Crispell
## Email: crispelljoseph@gmail.com
## Created: 12-10-2017


# Requirements

**FOR FULL INSTALLATION INSTRUCTIONS SEE installation.txt**

* Running on Ubuntu v18
* fastqc - tool to examine quality of FASTQ files
* cutadapt - tool to trim reads and remove adapters
* bwa - alignment algorithm
* samtools - package of tools to handle sequencing files
* blastn - tool to compare a sequence to a database of sequences
* bcftools - package of tools to handle sequencing files
* java jdk - programming language necessary to run some of the scripts
* perl - programming language necessary to run some of the scripts


# Downloading and Indexing Reference Sequence

The recently updated M. bovis reference sequence is available here: https://www.ncbi.nlm.nih.gov/nuccore/LT708304
NOTE - the updated annotation file doesn't contain the repeat_regions or PE/PPE feature annotations
I have an updated version of this annotation file and am working towards replacing these annotations
The old sequence with working annotation file can be found here: https://www.ncbi.nlm.nih.gov/nuccore/BX248333.1?report=genbank

With the FASTA of the reference genome downloaded run the following command to index it:
```
bwa index reference.fasta
```


# Instruction for running scripts

## Stage 1 - Examine FASTQ file quality and define trimming parameters

Move to directory containing forward and reverse FASTQ files
Run `fastqc` using:
```
fastqc --threads 4 *.fastq.gz
```

Examine each of the html files corresponding to each of the FASTQ files
Decide upon appropriate trimming parameters and enter them into a text file 
called PrinseqSettings_DATE.txt (see notes at bottom of page)

## Stage 2 - Align the FASTQ files against reference and create VCF files
Move to directory containing the forward and reverse FASTQ files and PrinseqSettings_DATE.txt file
Run `ProcessRawReads_DATE.sh` using:
```
bash ProcessRawReads_DATE.sh [FASTQFileEnding] [pathToCutadapt|false] [pathToReferenceSequence] [pathToPickRandomReads_DATE.pl] [pathToExamineBlastOutput_DATE.pl]
```
*NOTE:* Set the trimming parameters for cutadapt by creating a file called `Cutadapt_[DATE].txt` with the following structure:
```
UNIVERSALADAPTER1="AGATCGGAAGAG" # Illumina universal adapter sequence (source: https://github.com/s-andrews/FastQC/blob/master/Configuration/adapter_list.txt)
UNIVERSALADAPTER2="CTGTCTCTTATA" # Nextera Transposase sequence (source:https://github.com/s-andrews/FastQC/blob/master/Configuration/adapter_list.txt)
LENGTH=75                        # The minimum length of read to be accepted
QUAL=25                          # Trim low-quality bases from 5' and 3' ends
PROPN=0.5                        # Discard reads with more than this proportion of Ns
TRIML=10                         # Trim sequence at the 5' end by x positions
TRIMR=0                          # Trim sequence at the 3' end by x positions
```

## Stage 3 - Merge the VCF files together
Move to directory containing the VCF files resulting from the previous step - "vcfFiles"
Run `MergeVCFFiles_DATE.jar` using:
```
java -jar MergeVCFFiles_DATE.jar [pathToDirectoryToRunToolIn] [pathToAnnotationsOfRefenceGenome]
```

## Stage 4 - Examine genome coverage of isolates
Run `ExamineGenomeCoverage_DATE.pl` using:
```
perl ExamineGenomeCoverage_DATE.pl [ReadDepthThreshold] genomeCoverage_DATE.txt
```

## Stage 5 - Filter the merged VCF file
Run `FilterVariants_DATE.pl` using:
```
perl FilterVariants_DATE.pl 1 [ReadDepth] [HighQualityBaseDepth] [MappingQuality] [ProportionReadsSupportingAllele] [SiteCoverageAcrossIsolates] [GeneralQuality] [FQ] merged_DATE.txt
```

## Stage 6 - Rescuing site information for poorer quality isolates
Run `RescueVariantPositionInfo_DATE.pl` using:
```
perl RescueVariantPositionInfo_DATE.pl [nIsolatesAlleleMustBePresentIn] [HighQualityBaseDepth] [ProportionReadsSupportingAllele] filtered_DATE.txt
```

## Stage 7 - Calculate isolate coverage in filtered OR filtered-rescued file
Run `CalculateIsolateCoverageFromFiltered_DATE.pl` using:
```
perl CalculateIsolateCoverageFromFiltered_DATE.pl filtered_DATE.txt
```

## Stage 8 - Create FASTA file from filtered
Run `CreateFastaFromFiltered_DATE.pl` using:
```
perl CreateFastaFromFiltered_DATE.pl 1 [MinNumberSitesBetweenVariantPositions] [pathToReferenceSequence] filtered-rescued_DATE.txt
```


# NOTES 

If you have any problems please drop me an email at: crispelljoseph@gmail.com

For additional help (including installation instructions of required tools) with particular
scripts use "-help" directly after script name

Full installation instructions for the require tools for `ProcessRawReads_DATE.sh` can be found
in the "checkToolsAreInstalled" function at the top of the script
