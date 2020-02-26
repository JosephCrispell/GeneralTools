#### Preparation ####

# Load the required libraries
library(ape) # Reading in FASTA
library(lubridate) # Converting to decimal dates

# Set the path
path <- file.path("~", "Desktop", "BuildingBASTAXML")

# Set the date
date <- format(Sys.Date(), "%d-%m-%y")

#### Read in FASTA, metadata and constant site counts ####

# Read in the FASTA file
fastaFile <- file.path(path, "example.fasta")
sequences <- read.dna(fastaFile, format="fasta", as.character=TRUE)

# Read in the metadata
metadataFile <- file.path(path, "metadata.csv")
metadata <- read.table(metadataFile, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Convert sampling dates to decimal dates
metadata$DecimalDate <- decimal_date(as.Date(metadata$Date))

# Read in the constant site counts                          A,C,G,T proportions: 0.175, 0.325, 0.325, 0.175
# - I've estimated these: (genomeSize - nSites) split between A, C, G, T with 65% bias to G and C
constantSiteCounts <- generateConstantSiteCounts(nSites=ncol(sequences))

#### Build BASTA xml(s) ####

# Set the options to be used in the BASTA analyses
popSizeEstimation <- "equal" # "equal" or "varying" ?
relaxedOrStrict <- "relaxed" # "relaxed" or "strict" ?
chainLength <- 300000000
sampleFromPrior <- FALSE

# Build the BASTA xml


#### FUNCTIONS - preparation ####

generateConstantSiteCounts <- function(nSites, genomeSize=4345492, nucleotideProbs=c(0.175, 0.325, 0.325, 0.175)){
  
  # Calculate the number if sites not considered in the FASTA file
  nGenomeSitesNotConsidered <- genomeSize - nSites
  
  # Spread these not considered sites into nucleotide counts accordinag to therii respective probabilities
  nucleotideCounts <- ceiling(nucleotideProbs * nGenomeSitesNotConsidered)
  
  return(nucleotideCounts)
}

generateRandomMetadata <- function(nSequences, path){
  
  # Generate random data for sequences
  metadata <- data.frame("ID"=rownames(sequences), 
                         "Date"=sample(seq(from=as.Date("01-01-2005", format="%d-%m-%Y"),
                                           to=as.Date("01-12-2015", format="%d-%m-%Y"), by=1),
                                       size=nSequences, replace=TRUE),
                         "Species"=sample(c("badger", "cow"), size=nSequences, replace=TRUE))
  
  # Write the data to file
  metadataFile <- file.path(path, "metadata.csv")
  write.table(metadata, file=metadataFile, quote=FALSE, sep=",", row.names=FALSE)
}

#### FUNCTIONS - build BASTA xml ####

addDemeAssignmentBlock <- function(fileLines, metadata){
  
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Deme Assignment -->"
  output <- "\t<typeTraitSet id=\"typeTraitSet\" spec=\"TraitSet\" traitname=\"type\" value=\""
  for(row in 1:nrow(metadata)){
    
    name <- paste0(metadata[row, "ID"], "_", metadata[row, "Species"])
    
    output <- paste0(output, name, "=", metadata[row, "Species"])
    
    if(row < nrow(metadata)){
      output <- paste(output, ",", sep="")
    }
  }
  output <- paste(output, "\">", sep="")
  fileLines[length(fileLines) + 1] <- output
  fileLines[length(fileLines) + 1] <- "\t\t<taxa spec='TaxonSet' alignment='@alignment'/>"
  fileLines[length(fileLines) + 1] <- "\t</typeTraitSet>"
  
  return(fileLines)
}

addTipDateBlock <- function(fileLines, metadata){
  
  # Print out tip dates
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Tip Dates -->"
  fileLines[length(fileLines) + 1] <- "\t<timeTraitSet spec='beast.evolution.tree.TraitSet' id='timeTraitSet' traitname=\"date-forward\""
  output <- "\t\tvalue=\""
  for(row in 1:nrow(metadata)){
    name <- paste0(metadata[row, "ID"], "_", metadata[row, "Species"])
    
    output <- paste0(output, name, "=", metadata[row, "DecimalDate"])
    
    if(row < nrow(metadata)){
      output <- paste(output, ",", sep="")
    }
  }
  output <- paste(output, "\">", sep="")
  fileLines[length(fileLines) + 1] <- output
  fileLines[length(fileLines) + 1] <- "\t\t<taxa spec='TaxonSet' alignment='@alignment'/>"
  fileLines[length(fileLines) + 1] <- "\t</timeTraitSet>"
  
  return(fileLines)
}

addConstantSiteCountsBlock <- function(fileLines, counts){
  
  ## Print out block
  counts <- paste(counts, collapse=" ")
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Constant Site Counts -->"
  fileLines[length(fileLines) + 1] <- paste("\t<data id='alignment' spec='FilteredAlignment' filter='-' data='@alignmentVar' constantSiteWeights='",
                                            counts, "'/>", sep="")
  
  return(fileLines)
}

addDistributionBlock <- function(fileLines){
  fileLines[length(fileLines) + 1] <- ""
  fileLines[length(fileLines) + 1] <- "\t<!-- Initialise distributions to refer to -->"
  fileLines[length(fileLines) + 1] <- "\t<map name=\"Exponential\" >beast.math.distributions.Exponential</map>"
  
  return(fileLines)
}

addSequenceBlock <- function(fileLines, sequences, metadata, sampleFromPrior){
  
  # Start the sequence block
  fileLines[length(fileLines) + 1] <- "\t<!-- Sequence Alignment -->"
  if(sampleFromPrior == FALSE){
    fileLines[length(fileLines) + 1] <- "\t<data id=\"alignmentVar\" dataType=\"nucleotide\">"
  }else{
    fileLines[length(fileLines) + 1] <- "\t<data id=\"alignment\" dataType=\"nucleotide\">"
  }
  
  # Add each of the sequences into the block
  for(row in 1:nrow(metadata)){
    
    # Build the sequence name
    name <- paste0(metadata[row, "ID"], "_", metadata[row, "Species"])
    
    # Build the sequence
    sequence <- paste(sequences[metadata[row, "ID"], 1:ncol(sequences)], collapse="")
    
    # Build the file lines using the name and sequence
    fileLines[length(fileLines) + 1] <- paste0("\t\t<sequence taxon=\"", name, "\">")
    fileLines[length(fileLines) + 1] <- paste0("\t\t\t", ifelse(sampleFromPrior, "N", sequence))
    fileLines[length(fileLines) + 1] <- "\t\t</sequence>"
  }
  
  # Close the sequence block
  fileLines[length(fileLines) + 1] <- "\t</data>"
  
  return(fileLines)
}

startBuildingOutputFileLines <- function(){
  
  # Create a vector to store the fileLines
  fileLines <- c()
  
  # Insert the first two lines
  fileLines[1] <- "<beast version='2.0' namespace='beast.evolution.alignment:beast.core:beast.core.parameter:beast.evolution.tree:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood:beast.evolution.tree:beast.math.distributions:multitypetreeVolz.distributions:multitypetreeVolz.operators:multitypetreeVolz.util'>"
  fileLines[2] <- ""
  
  return(fileLines)
}

buildXMLFile <- function(sequences, metadata, equalOrVaryingPopSizes, path, date, 
                         constantSiteCounts, chainLength,
                         relaxedOrStrict, sampleFromPrior, runNumber=1){
  
  # Define the deme structure:
  # Demes:
  #   0 badger
  #   1 cow
  #
  #        | badger | cow | 
  # _______|________|_____|
  # badger |   0    |  0  |
  # _______|________|_____|
  # cow    |   0    |  0  |
  # _______|________|_____|
  migrationRateMatrix <- matrix(c(NA, 1,
                                  1, NA),
                                byrow=TRUE, nrow=2)
  
  # Note the output XML name
  outputFileName <- paste("BASTA_", equalOrVaryingPopSizes, "_", relaxedOrStrict,
                          "_", date, sep="")
  if(sampleFromPrior == TRUE){
    outputFileName <- paste("BASTA_", equalOrVaryingPopSizes, "_", relaxedOrStrict,
                            "_PRIOR_", date, sep="")
  }
  
  # Create a directory for the output file
  dir.create(paste(path, outputFileName, "/", sep=""))
  
  # Create an array of file lines and add in each necessary block
  fileLines <- startBuildingOutputFileLines()
  
  # Add Sequenced block
  fileLines <- addSequenceBlock(sequences, metadata, sampleFromPrior)
  
  # Add distribution block - if relaxed
  if(relaxedOrStrict == "relaxed"){
    fileLines <- addDistributionBlock(fileLines)
  }
  
  # Add constant site counts block
  if(sampleFromPrior == FALSE){
    fileLines <- addConstantSiteCountsBlock(fileLines, constantSiteCounts)
  }
  
  # Add sampling dates block
  fileLines <- addTipDateBlock(fileLines, metadata)
  
  # Add deme assignment block - I AM HERE !?!?!?!?!?!?!?!?!?!?!?
  fileLines <- addDemeAssignmentBlock(fileLines, metadata)
  
  # Add substitution model block
  fileLines <- addSubstitutionModelBlock(fileLines, tipInfo, relaxedOrStrict,
                                         sampleFromPrior)
  
  # Add branch rates model - if relaxed
  if(relaxedOrStrict == "relaxed"){
    fileLines <- addBranchRateModelBlock(fileLines)
  }
  
  # Add migration model block
  fileLines <- addMigrationModelBlock(fileLines, demeStructure, initialValue=0.1, 
                                      migrationRateMatrix, nDemes=ncol(migrationRateMatrix))
  
  # Add Prior distributions block
  fileLines <- addPriorsBlock(fileLines, nDemes=ncol(migrationRateMatrix), relaxedOrStrict)
  
  # Add Tree likelihood block
  fileLines <- addTreeLikelihoodBlock(fileLines, relaxedOrStrict)
  
  # Add Migration rate likelihood block
  fileLines <- addMigrationRateLikelihoodBlock(fileLines)
  
  # Add Beast settings block
  fileLines <- addBeastSettingsBlock(fileLines, demeStructure,
                                     equalOrVaryingPopSizes == "varying", outputFileName,
                                     initialValue = 0.1, chainLength, relaxedOrStrict,
                                     migrationRateMatrix, nDemes=ncol(migrationRateMatrix))
  
  # Add the end to the XML file
  fileLines <- addEnd(fileLines)
  
  # Write the file out
  writeXMLFile(fileLines, paste(path, outputFileName, "/", outputFileName, ".xml", sep=""))
}

