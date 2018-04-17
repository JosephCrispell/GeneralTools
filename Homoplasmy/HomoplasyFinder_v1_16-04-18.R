#-------------------#
#### Preparation ####
#-------------------#

# Packages
library(geiger) # For the tips function
library(ape) # read.dna()

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/vcfFiles/"

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

#------------------------------#
#### Read in tree and FASTA ####
#------------------------------#

# Read in the tree
treeFile <- paste(path, "mlTree_27-03-18.tree", sep="")
tree <- read.tree(treeFile)

# Note the isolates at the tips of each node
nodes <- getNodes(tree)

# Read in the FASTA file
fastaFile<- paste(path,"sequences_Prox-10_24-03-2018.fasta", sep="")
sequences <- as.alignment(read.dna(fastaFile, format="fasta", skip=1))

# Record the alleles present in FASTA sequences
alleles <- recordAllelesInPopulation(sequences, TRUE)

#--------------------------------------------#
#### Assign alleles to nodes in phylogeny ####
#--------------------------------------------#

# Assign alelles to nodes in the phylogeny - where possible
notAssigned <- assignAlellesToNodes(nodes, alleles, sequences$nam, TRUE)

#-----------------------------------------#
#### Report the homopalsies identified ####
#-----------------------------------------#

results <- reportHomoplasiesIdentified(notAssigned, alleles)

#-----------------#
#### FUNCTIONS ####
#-----------------#

reportHomoplasiesIdentified <- function(notAssigned, alleles){
  
  # Get the positions where homoplasies were identified
  positions <- getPositions(notAssigned)
  
  # Initialise a data frame to store the homoplasy information
  output <- data.frame("Position"=names(positions), "Alleles"=NA, "IsolatesForAlleles"=NA, stringsAsFactors=FALSE)
  
  # Report each homoplasy
  for(row in 1:length(positions)){
    
    # Store the position
    output[row, "Position"] <- names(positions)[row]
    
    # Store the nucleotides associated with the current homoplasy's position
    nucleotides <- positions[[names(positions)[row]]]
    output[row, "Alleles"] <- paste(toupper(nucleotides), collapse=",")
    
    # Get the isolates associated with each nucleotide observed at the current site
    isolates <- paste(alleles[[paste(names(positions)[row], nucleotides[1], sep=":")]], collapse="-")
    for(i in 2:length(nucleotides)){
      
      isolates <- paste(isolates, paste(alleles[[paste(names(positions)[row], nucleotides[1], sep=":")]], collapse="-"), sep=",")
    }
    output[row, "IsolatesForAlleles"] <- isolates
  }
  
  return(output)
}

getPositions <- function(alleles){
  
  positions <- list()
  for(allele in alleles){
    
    parts <- strsplit(allele, split=":")[[1]]
    position <- parts[1]
    if(is.null(positions[[position]]) == FALSE){
      positions[[position]] <- c(positions[[position]], parts[2])
    }else{
      positions[[position]] <- c(parts[2])
    }
  }
  
  return(positions)
}

assignAlellesToNodes <- function(nodes, alleles, isolates, verbose){
  
  # Initialise a hashtable to record which alleles were assigned
  nodeAlleles <- list()
  
  # Initialise a vector to store the assigned alleles
  assigned <- c()

  # Examine each node
  for(i in 1:length((nodes))){

    # Progress
    if(verbose){
      cat(paste("\rAssigning alleles to node", i, "of", length(nodes)))
    }

    # Get the isolates at the tips of the current node
    tips <- nodes[[names(nodes)[i]]]
    
    # Get the rest of the isolates
    isolatesBelow <- isolates[isolates %in% tips == FALSE]

    # Check and see if these isolates are associated with an allele
    for(allele in names(alleles)){
      
      # Skip alleles that have already been assigned
      if(allele %in% assigned){
        next
      }
      
      # Get the isolates with an "N" at the current allele's position
      isolatesWithN <- getIsolatesWithNsAtPosition(allele, alleles)
      
      # Remove the isolates with an "N" from those associated with the current node
      tipsWithoutNs <- tips[tips %in% isolatesWithN == FALSE]
      isolatesBelowWithoutNs <- isolatesBelow[isolatesBelow %in% isolatesWithN == FALSE]
      
      # Get the isolates with the current allele
      isolatesWithAllele <- alleles[[allele]]
      
      # Compare the two sets of alleles - if they match exactly then current allele can be assigned to node
      if(areSetsOfIsolatesTheSame(tipsWithoutNs, isolatesWithAllele) == TRUE){

        nodeAlleles[[names(nodes)[i]]] <- c(nodeAlleles[[names(nodes)[i]]], allele)
        assigned[length(assigned) + 1] <- allele
      
      }else if(areSetsOfIsolatesTheSame(isolatesBelowWithoutNs, isolatesWithAllele) == TRUE){
        nodeAlleles[[names(nodes)[i]]] <- c(nodeAlleles[[names(nodes)[i]]], allele)
        assigned[length(assigned) + 1] <- allele
      }
    }
  }
  
  if(verbose){
    cat("\rFinished assigning alleles to nodes.\t\t\t\t\n")
  }
  
  # Note the alleles that weren't assigned
  notAssigned <- names(alleles)[names(alleles) %in% assigned == FALSE]
  
  return(notAssigned)
}

getIsolatesWithNsAtPosition <- function(allele, alleles){
  
  # Build the allele with an N
  key <- paste(strsplit(allele, split=":")[[1]][1], "n", sep=":")
  
  # Check whether any isolates with an N were noted
  isolates <- c()
  if(is.null(alleles[[key]]) == FALSE){
    isolates <- alleles[[key]]
  }
  
  return(isolates)
}

areSetsOfIsolatesTheSame <- function(a, b){
  result <- FALSE
  if(length(a) == length(b) && length(intersect(a, b)) == length(a)){
    result <- TRUE
  }
  return(result)
}

getNodes <- function(tree){
  nodes <- list()
  
  # Nodes number 1:nTips and then onwards: two methods to calculate total number:
  # - Number of edges + 1
  # - Number of tips plus number internal nodes: length(tree$tip.label) + tree$Nnode
  for(node in 1:(length(tree$tip.label) + tree$Nnode)){
    nodes[[as.character(node)]] <- tips(tree, node)
  }
  
  return(nodes)
}

recordAllelesInPopulation <- function(sequences, verbose){
  
  # Initialise a list to store the alleles found and the sequences they are found in
  alleles <- list()
  
  # Examine each of the sequences
  for(i in 1:length(sequences$nam)){
    
    # Progress
    if(verbose){
      cat(paste("\rReading alleles in sequence", i, "of", length(sequences$nam)))
    }
    
    # Split the current sequence into nucleotides
    nucleotides <- strsplit(sequences$seq[i], split="")[[1]]
    
    # Examine each nucleotide in the current sequence
    for(pos in 1:length(nucleotides)){
      
      # Build an allele ID
      id <- paste(pos, nucleotides[pos], sep=":")
      
      # Check if encountered this allele before
      if(is.null(alleles[[id]]) == FALSE){
        
        alleles[[id]] <- c(alleles[[id]], sequences$nam[i])
      }else{
        alleles[[id]] <- c(sequences$nam[i])
      }
    }
  }
  
  if(verbose){
    cat("\rFinished reading alleles from sequences.\t\t\t\t\n")
  }
  
  return(alleles)
}
