############
# Packages #
############

library(ape)

### NOTES - WINDOWS
# Needs clustalw2.exe to be in path
# Install from here: http://www.clustal.org/download/current/
# Add to path on Windows machine by:
# - Find on computer after installation - Program Files
# - Copy path
# - Start -> Control Panel -> System -> Advanced system settings
#   -> Environment Variables -> Select "Path" from System variables
# - Go to end, type ";", paste path

### NOTES - UBUNTU
# sudo apt upgrade
# sudo apt install clustalw

######################
# Read in FASTA file #
######################

# Set the path
path <- "/home/josephcrispell/Desktop/"

# Read in the sequences
fastaFile <- paste(path, "test.fasta", sep="")
sequences <- read.FASTA(fastaFile)

###########################################
# Run Clustal multiple sequence alignment #
###########################################

# Run the alignment
alignment <- clustal(sequences)

################################
# Build the consensus sequence #
################################

# Add the consensus sequence into the alignment
alignment <- getConsensus(alignment)

# Print out sequences with consensus
outputFile <- paste(path, "sequencesWithConsensus.fasta", sep="")
write.FASTA(alignment, file=outputFile)

#############
# FUNCTIONS #
#############

getConsensus <- function(sequences){
  
  # Convert between ways of storing sequences
  sequences <- as.alignment(sequences)
  
  # Initialise vector to store the consensus
  consensus <- c()
  
  # Get the sequences and split them into their nucleotides
  nucleotides <- list()
  for(i in 1:sequences$nb){
    nucleotides[[sequences$nam[i]]] <- strsplit(sequences$seq[i], split="")[[1]]
  }
  
  # Examine each position in the alignment
  for(position in 1:length(nucleotides[[1]])){
    
    # Initialise a vector to store the alleles at the current position
    alleles <- c()
    
    # Examine the current position in each sequence - store each allele
    for(i in 1:sequences$nb){
      alleles[i] <- nucleotides[[i]][position]
    }
    
    # Get the number of each allele present
    alleleCounts <- table(alleles)
      
    # Note the consensus
    maxAllele <- names(which(alleleCounts == max(alleleCounts)))
    if(length(maxAllele) == 1){
      consensus[position] <- maxAllele
    }else{
      consensus[position] <- "n"
    }
  }
  
  # Attach the consensus to the alignment
  sequences$seq[sequences$nb + 1] <- paste(consensus, collapse="")
  sequences$nb <- sequences$nb + 1
  sequences$nam <- c(sequences$nam, "consensus")
  
  # Convert between ways of storing sequences
  sequences <- as.DNAbin.alignment(sequences)
  
  return(sequences)
}
