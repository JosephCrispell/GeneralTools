###############
# Preparation #
###############

# Packages
library(phangorn)
library(gplots)

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Homoplasmy/"

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

######################
# Generate sequences #
######################

# Generate the sequences
sequences <- simulateSequences(n=500, genomeSize=1000000)

####################
# Insert homoplasy #
####################

output <- insertHomoplasy(sequences, ancestralSite="A")
sequences <- output[["sequences"]]
ij <- c(output[["source"]], output[["sink"]])

###############
# Build FASTA #
###############

writeFasta(sequences, paste(path, "example_", date, ".fasta", sep=""))

###################
# Build phylogeny #
###################

buildPhylogeny(sequences, paste(path, "example_", date, ".tree", sep=""), ij)

#############
# FUNCTIONS #
#############

removeUninformativeSites <- function(sequences){
  
  # Initialise a vector to record which sites to remove
  sitesToKeep <- rep(FALSE, length(sequences[[1]]))
  
  # Get the list of isolates
  keys <- names(sequences)
  
  # Examine each site
  for(position in 1:length(sequences[[1]])){
    
    # Compare the first isolate to all others - search for difference = informative
    for(index in 2:length(keys)){
      
      # Does the current isolate differ from the first at the current position?
      if(sequences[[keys[1]]][position] != sequences[[keys[index]]][position]){
        sitesToKeep[position] <- TRUE
        break
      }
    }
    
    cat(paste("\rChecking for uninformative sites... (", position,
              " of ", length(sequences[[1]]), ")", sep=""))
  }
  
  # Remove the uninformative sites from each isolate
  for(key in keys){
    sequences[[key]] <- sequences[[key]][sitesToKeep]
  }
  
  return(sequences)
}

buildPhylogeny <- function(sequences, file, ij){
  
  par(mar=c(0,0,0,0))
  
  # Build genetic distance matrix
  geneticDistances <- buildGeneticDistanceMatrix(sequences)
  
  # Build neighbour joining tree
  tree <- NJ(geneticDistances)
  
  # Plot the tree
  plot.phylo(tree, show.tip.label=TRUE, type="phylogram",
             edge.color="grey", edge.width=3,
             show.node.label=TRUE, label.offset=0.15,
             tip.color=ifelse(tree$tip.label %in% ij, "red", "black"))
  
  # Write the tree to file
  write.tree(tree, file=file, append=FALSE,
             digits=20, tree.names=FALSE)
}

insertHomoplasy <- function(sequences, ancestralSite){
  i <- sample(size=1, x=names(sequences))
  j <- sample(size=1, x=names(sequences)[names(sequences) != i])
  
  positionsToChooseFrom <- which(sequences[[i]] != sequences[[j]] &
                                   sequences[[i]] != ancestralSite)
  
  position <- positionsToChooseFrom[
    sample(size=1, x=1:length(positionsToChooseFrom))]
  
  sequences[[j]][position] <- sequences[[i]][position]
  
  cat(paste("Introduced homoplasy in isolate ", j, " at position ",
            position, " from isolate ", i, " with allele ", 
            sequences[[i]][position], "\n", sep=""))
  
  output <- list(
    "sequences" = sequences,
    "source" = i,
    "sink" = j,
    "position" = position,
    "allele" = sequences[[i]][position]
  )
  
  return(output)
}

simulateSequences <- function(n, genomeSize){
  # Note the nucleotides
  nucleotides <- c("A", "T", "G", "C")
  
  # Set population size
  susceptibles <- 1:n
  
  # Choose random seed
  sources <- c(sample(size=1, susceptibles))
  susceptibles <- susceptibles[susceptibles != sources[1]]
  
  # Store each individual's sequence
  sequences <- c()
  sequences[[as.character(sources[1])]] <- rep("A", genomeSize)
  
  # Run infection
  while(length(susceptibles) > 0){
    
    # Choose a source
    source <- sources[sample(size=1, x=1:length(sources))]
    
    # Choose a sink
    sink <- susceptibles[sample(size=1, x=1:length(susceptibles))]
    sources[length(sources) + 1] <- sink
    susceptibles <- susceptibles[susceptibles != sink]
    
    # Get sources sequence
    sequence <- sequences[[as.character(source)]]
    
    # Mutate sequence
    position <- sample(size=1, x=1:length(sequence))
    sequence[position] <- sample(size=1, x=nucleotides[nucleotides != sequence[position]])
    
    # Transfer sequence
    sequences[[as.character(sink)]] <- sequence
    
    cat(paste(source, " -> ", sink, "\n"))
  }
  
  # Remove uninformative sites
  sequences <- removeUninformativeSites(sequences)
  
  return(sequences)
}

writeFasta <- function(sequences, fileName){
  
  # Open file
  fileConnection <- file(fileName, open="w")
  
  # Write header
  writeLines(paste(length(sequences),
                   length(sequences[[1]])), fileConnection)
  
  # Write the sequences
  for(key in names(sequences)){
    writeLines(paste(">", key, "\n", 
                     paste(sequences[[key]], collapse=""), sep=""),
               fileConnection)
  }
  
  # Close the file
  close(fileConnection)
}

buildGeneticDistanceMatrix <- function(sequences){
  
  matrix <- matrix(nrow=length(sequences), ncol=length(sequences))
  
  
  keys <- names(sequences)
  rownames(matrix) <- keys
  colnames(matrix) <- keys
  
  for(i in 1:length(sequences)){
    
    for(j in 1:length(sequences)){
      
      if(i >= j){
        next
      }
      
      distance <- geneticDistance(sequences[[keys[i]]], sequences[[keys[j]]])
      matrix[i, j] <- distance
      matrix[j, i] <- distance
    }
  }
  
  return(matrix)
}

geneticDistance <- function(a, b){
  
  distance <- 0
  
  for(i in 1:length(a)){
    
    if(a[i] != "N" && b[i] != "N" && a[i] != b[i]){
      
      distance <- distance + 1
    }
  }
  
  return(distance)
}