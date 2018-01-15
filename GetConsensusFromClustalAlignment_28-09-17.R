############
# Packages #
############

library(ape)

### NOTES
# Needs clustalw2.exe to be in path
# Install from here: http://www.clustal.org/download/current/
# Add to path on Windows machine by:
# - Find on computer after installation - Program Files
# - Copy path
# - Start -> Control Panel -> System -> Advanced system settings
#   -> Environment Variables -> Select "Path" from System variables
# - Go to end, type ";", paste path

######################
# Read in FASTA file #
######################

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/"

# Read in the sequences
fastaFile <- paste(path, "test.fasta", sep="")
sequences <- read.FASTA(fastaFile)

###########################################
# Run Clustal multiple sequence alignment #
###########################################

# Run the alignment
alignment <- clustal(sequences)

# Print out the file
output <- paste(path, "alignment.fasta", sep="")
write.dna(alignment, file = output, format = 'fasta', colsep="")

################################
# Build the consensus sequence #
################################

# Read in the alignment
sequences <- readFASTA(output)

# Get the consensus sequence
consensusSequence <- getConsensus(sequences)

# Print out the consensus sequence
outputFile <- paste(path, "Consensus.fasta", sep="")
printSequence(consensusSequence, outputFile)

#############
# FUNCTIONS #
#############

printSequence <- function(consensus, fileName){
  
  # Open an output file
  connection <- file(fileName, "w")
  
  # Write a sequence header
  writeLines(">Consensus", connection)
  
  # Write out the sequence
  writeLines(paste(consensus, collapse=""), connection)
  
  # Close the output file
  close(connection)
}

getConsensus <- function(sequences){
  
  # Get the sequence names
  sequenceNames <- names(sequences)
  
  # Initialise vector to store the consensus
  consensus <- c()
  
  # Examine each position in the alignment
  for(position in 1:length(sequences[[1]])){
    
    # Initialise a vector to store the alleles at the current position
    alleles <- c()
    
    # Examine the current position in each sequence - store each allele
    for(i in 1:length(sequenceNames)){
      alleles[i] <- sequences[[sequenceNames[i]]][position]
    }
    
    # Get the number of each allele present
    alleleCounts <- table(alleles)
      
    # Note the consensus
    maxAllele <- names(which(alleleCounts == max(alleleCounts)))
    if(length(maxAllele) == 1){
      consensus[position] <- maxAllele
    }else{
      consensus[position] <- "N"
    }
  }
  
  return(consensus)
}

readFASTA <- function(fastaFile){
  
  # Open the file and store lines
  connection <- file(fastaFile, "r")
  fileLines <- readLines(connection)
  close(connection)
  
  # Initialise a list to store the sequences
  sequences <- list()
  
  # Examine each of the file lines one by one
  for(i in 1:length(fileLines)){
    
    # Check if line starts with ">"
    if(grepl(fileLines[i], pattern="^>") == TRUE){
      
      # Store the previous sequence
      if(i != 1){
        sequences[[name]] <- strsplit(sequence, split="")[[1]]
      }
      
      # Initialise variables to store sequence and id
      name <- substr(fileLines[i], 2, stop=nchar(fileLines[i]))
      sequence <- ""
    }else{
      
      sequence <- paste(sequence, fileLines[i], sep="")
    }
  }
  
  # Store the last sequence
  sequences[[name]] <- strsplit(sequence, split="")[[1]]
  
  return(sequences)
}
