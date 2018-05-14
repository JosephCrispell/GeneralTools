##################
# Load libraries #
##################

library(devtools)
#install_github("JosephCrispell/homoplasyFinder")
library(homoplasyFinder)
library(ape)

################################################
# Run homoplasyFinder on all examples datasets #
################################################

# Set the path
path <- "/home/josephcrispell/Desktop/Research/Homoplasy/TimeTrial/"

# Set the date when datasets were created
date <- "11-05-18"

# Note the range in the number of sequences in each dataset and number of replicates
nSequences <- seq(50, 500, 50)
nReplicates <- 10

# Initialise a table to record the time taken for homoplasyFinder to analyses each dataset
timeTaken <- data.frame("NSequences"=NA, "NReplicates"=NA, "User"=NA, "System"=NA, "Elapsed"=NA, stringsAsFactors=FALSE)
row <- 0

# Build each file name and run in homoplasyFinder
for(n in nSequences){
  
  for(replicate in 1:nReplicates){
    
    # Progress information
    cat(paste("\rAnalysing data with", n, "sequences. Replicate:", replicate, "      "))
    
    # Increment row
    row <- row + 1
    
    # Build file name prefix
    prefix <- paste(path, "Example_", (n*2), "-", n, "_", replicate, "_", date, sep="")
    fastaFile <- paste(prefix, "fasta", sep=".")
    treeFile <- paste(prefix, "tree", sep=".")

    # Run homoplasy finder
    time <- system.time(run(fastaFile, treeFile))
    
    # Store the time taken
    timeTaken[row, ] <- c(n, replicate, time[[1]], time[[2]], time[[3]])
  }
}

# Write the table to file
file <- paste(path, "timeTaken-R_", date, ".csv", sep="")
write.table(timeTaken, file, row.names=FALSE, quote=FALSE, sep=",")

#############
# FUNCTIONS #
#############

run <- function(fastaFile, treeFile){
  
  sequences <- read.dna(fastaFile, skip=1, format="fasta")
  tree <- read.tree(treeFile)
  
  results <- homoplasyFinder(tree, sequences, verbose=FALSE)
}
