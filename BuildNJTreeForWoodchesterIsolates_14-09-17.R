# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_13-07-17/"
############
# Packages #
############

library(phangorn)

###################
# Build phylogeny #
###################

# Read in the FASTA file
file <- paste(path, "vcfFiles/sequences_Prox-10_01-08-2017.fasta", sep="")
sequences <- readFasta(file)

# Build the genetic distance matrix
geneticDistances <- buildGeneticDistanceMatrix(sequences)

# Build a neighbour joining tree
tree <- NJ(geneticDistances)

# Plot the tree
plot.phylo(tree, show.tip.label=FALSE, type="fan",
           edge.color="dimgrey", edge.width=3,
           show.node.label=TRUE)











#############
# FUNCTIONS #
#############

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

readFasta <- function(fileName){
  
  # Store all file lines
  connection <- file(fileName, open="r")
  fileLines <- readLines(connection)
  close(connection)
  
  # Initialise a list to store the fasta sequences
  sequences <- list()
  
  # Examine each line - skip first line
  for(i in 2:length(fileLines)){
    
    # Check if sequence header
    if(startsWith(fileLines[i], prefix=">") == TRUE){
      
      # Store previous sequence
      if(i != 2){
        sequences[[name]] <- strsplit(sequence, split="")[[1]]
      }
      
      # Get sequence name
      name <- substr(fileLines[i], start=2, stop=nchar(fileLines[i]))
      name <- strsplit(name, split="_")[[1]][1]
      
      # Reset sequence
      sequence <- ""
    }else{
      sequence <- paste(sequence, fileLines[i], sep="")
    }
  }
  
  # Store last sequence
  sequences[[name]] <- strsplit(sequence, split="")[[1]]
  
  return(sequences)
}
