# Load packages
library(ape)
library(phangorn)

# Read in the VNTR and spoligotype information for your isolates
table <- read.table("fileName", header=TRUE, sep=",", stringsAsFactors=FALSE)

# Initialise a genetic distane matrix
distanceMatrix <- matrix(nrow=nrow(table), ncol=nrow(table))

# Compare every isolate to one another
for(i in 1:nrow(table)){
  
  for(j in 1:nrow(table)){
    
    # Skip the diagonal and making the smae comparison twice
    if(i >= j){
      next
    }
    
    # Compare the current pair of isolates
    distance <- 0
    for(column in 1:ncol(table)){
      if(table[i, column] != table[j, column]){
        distance <- distance + 1
      }
    }
    
    # Store that distance in the distance matrix
    distanceMatrix[i, j] <- distance
    distanceMatrix[j, i] <- distance
  }
}

# Build a phylogenetic tree
tree <- nj(distanceMatrix)

# Plot the tree
plot.phylo(tree)
