#### Libraries ####
library(parallel)

#### Initialise a cluster ####

# Get the number of threads in the current machine
nThreads <- parallel::detectCores()

# Initialise the cluster of threads
clusterOfThreads <- parallel::makeCluster(nThreads)
    
# Register the cluster of threads
doParallel::registerDoParallel(clusterOfThreads, cores=nThreads)

#### Start multi-threading ####

sequences <- createRandomNucleotideAlignment(300, 100000)

start <- Sys.time()
nucleotideFrequences <- clusterApply(cl=clusterOfThreads,
                                     x=sequences,
                                     fun=calculateNucleotideFrequencies)
print(Sys.time() - start)

start <- Sys.time()
for(i in 1:300){
  nucleotideFrequencies <- calculateNucleotideFrequencies(sequences[[i]])
}
print(Sys.time() - start)

# Close the cluster of threads
stopCluster(clusterOfThreads)


#### FUNCTIONS ####

calculateNucleotideFrequencies <- function(sequence){
  
  # Initialise a list to store the nucleotide counts
  frequencies <- list('A'=0, 'C'=0, 'G'=0, 'T'=0)
  
  # Split the sequence into its nucleotides
  nucleotides <- strsplit(sequence, split="")[[1]]
  
  # Count the number of times each nucleotide appears in the given sequence
  for(nucleotide in nucleotides){
    frequencies[[nucleotide]] <- frequencies[[nucleotide]] + 1
  }
  
  return(frequencies)
}

createRandomNucleotideAlignment <- function(n, length){
  
  # Initialise a dataframe to store the sequences
  sequences <- list()
  
  # Create each sequence
  for(i in 1:n){
    sequences[[i]] <- paste(sample(c('A', 'C', 'G', 'T'), size=length, replace=TRUE), collapse="")
  }
  
  return(sequences)
}
