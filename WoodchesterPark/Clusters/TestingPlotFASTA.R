# Load the basicPlotteR library
devtools::install_github("JosephCrispell/basicPlotteR")
library(basicPlotteR)

# Set the number of sequences
nSequences <- 25

# Set the number of sites in the sequences
nSites <- 200

# Create a random nucleotide alignment
sequences <- createNucleotideSequenceAlignment(nSequences, nSites, meanNSNPs=10)

# Plot the FASTA file
plotFASTA(sequences)


#### FUNCTIONS ####

createNucleotideSequenceAlignment <- function(nSequences, nSites, meanNSNPs){
  
  # Build a random reference genome
  reference <- sample(c("A", "C", "G", "T"), replace=TRUE, size=nSites)
  
  # Create a matrix to store the sequences
  sequences <- matrix(reference, byrow=TRUE, nrow=nSequences, ncol=nSites)
  sequences[1, ] <- reference
  rownames(sequences)<- paste0("Seq-", 1:nrow(sequences))
  
  # Create mutations in the sequences
  sequences <- mutateSequences(sequences, meanNSNPs)
  
  return(sequences)
}

mutateSequences <- function(sequences, meanNSNPs){
  
  for(row in 1:nrow(sequences)){
    
    # Decide how many mutations will be in the current sequence
    nMutations <- rpois(1, lambda=meanNSNPs)
    
    # Choose the position for the mutations
    positions <- sample(seq_len(ncol(sequences)), size=nMutations)
    
    # Create the mutations
    sequences[row, positions] <- sample(c("A", "C", "G", "T"), replace=TRUE, 
                                        size=nMutations)
  }
  
  return(sequences)
}