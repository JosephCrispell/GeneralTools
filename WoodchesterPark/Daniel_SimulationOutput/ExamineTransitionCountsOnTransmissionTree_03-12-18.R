#### Load the data ####

# Set the path
path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/BASTA/Simulations/"

# Get a list of sample IDs from the input XML file
xmlFile <- paste0(path, "sim_119-02_all_uniform/2Deme_varying_strict_02-12-18/", "2Deme_varying_strict_02-12-18.xml")
sampled <- getSampleIDs(xmlFile)

# Read in the full transmission tree
treeFile <- paste0(path, "sim_119-02/", "transmissionNetwork.0.edgeList")
edgeList <- read.table(treeFile, header=FALSE, sep="\t", stringsAsFactors=FALSE)

#### Build sampled transmission tree ####

# Convert edge list to adjacency matrix
adjacencyMatrix <- buildAdjacencyMatrix(edgeList)

# Build the sampled adjacency matrix
sampledAdjacencyMatrix <- buildSampledTransmissionTree(adjacencyMatrix, sampled)

#### Count transition events on full and sampled transmission trees ####

# Open a pdf
pdf(paste(path, "sim_119-02_uniform_02-12-18.pdf"), height=7, width=14)
par(mfrow=c(1,2))

# On the full transmission tree
countsOnFull <- countTransitions(adjacencyMatrix, plotAsMatrix=TRUE)

# On the sampled transmission tree
countsOnSampled <- countTransitions(sampledAdjacencyMatrix, plotAsMatrix=TRUE)

dev.off()

#### FUNCTIONS ####

getEdgeList <- function(adjacencyMatrix){
  
  # Intialise a dataframe to store the edges
  edges <- data.frame("Source"=NA, "Sink"=NA, stringsAsFactors=FALSE)
  count <- 0
  
  # Examine each row of the adjacency matrix
  for(row in seq_len(nrow(adjacencyMatrix))){

    # Examine each column of the adjacency matrix
    for(col in seq_len(ncol(adjacencyMatrix))){
      
      # Check if edge present
      if(adjacencyMatrix[row, col] == 1){
        count <- count + 1
        edges[count, "Source"] <- rownames(adjacencyMatrix)[row]
        edges[count, "Sink"] <- colnames(adjacencyMatrix)[col] 
      }
    }
  }
  
  return(edges)
}

plotCountsAsBars <- function(counts){
  
  # Get the set the margins
  currentMar <- par("mar")
  par(mar=c(7, 4.1, 1, 1))
  
  # Create a bar plot
  barplot(c(counts[1, 1], counts[1, 2], counts[1, 3],
            counts[2, 1], counts[2, 2], counts[2, 3],
            counts[3, 1], counts[3, 2], counts[3, 3]),
          space=c(0, 1, 1, 2, 1, 1, 2, 1, 1), las=1,
          ylab="Number of transitions")
  
  # Add labels
  axis(side=1, at=c(0.5, 2.5, 4.5, 7.5, 9.5, 11.5, 14.5, 16.5, 18.5),
       labels=c("root-root", "root-badger", "root-cow", 
                "badger-root", "badger-badger", "badger-cow",
                "cow-root", "cow-badger", "cow-cow"), 
       tick=FALSE, las=2)
  
  # Reset the margins
  par(mar=currentMar)
}

plotCountsAsMatrix <- function(counts){
  
  # Get the set the margins
  currentMar <- par("mar")
  par(mar=c(0, 4.1, 4.1, 0))
  
  # Create an empty plot
  plot(x=NULL, y=NULL, xlim=c(0,3), ylim=c(0,3), bty="n", axes=FALSE, xlab="", ylab="")
  
  # Add box labels
  axis(side=2, at=c(2.5, 1.5, 0.5), labels=c("root", "badger", "cow"), las=1, tick=FALSE)
  axis(side=3, at=c(0.5, 1.5, 2.5), labels=c("root", "badger", "cow"), las=1, tick=FALSE)
  
  # Note the maximum value
  maxCount <- max(counts)
  
  # Add polygons
  polygon(x=c(0, 0, 1, 1), y=c(2, 3, 3, 2), col=rgb(1,0,0, counts[1, 1] / maxCount))
  polygon(x=c(1, 1, 2, 2), y=c(2, 3, 3, 2), col=rgb(1,0,0, counts[1, 2] / maxCount))
  polygon(x=c(2, 2, 3, 3), y=c(2, 3, 3, 2), col=rgb(1,0,0, counts[1, 3] / maxCount))
  
  polygon(x=c(0, 0, 1, 1), y=c(1, 2, 2, 1), col=rgb(1,0,0, counts[2, 1] / maxCount))
  polygon(x=c(1, 1, 2, 2), y=c(1, 2, 2, 1), col=rgb(1,0,0, counts[2, 2] / maxCount))
  polygon(x=c(2, 2, 3, 3), y=c(1, 2, 2, 1), col=rgb(1,0,0, counts[2, 3] / maxCount))
  
  polygon(x=c(0, 0, 1, 1), y=c(0, 1, 1, 0), col=rgb(1,0,0, counts[3, 1] / maxCount))
  polygon(x=c(1, 1, 2, 2), y=c(0, 1, 1, 0), col=rgb(1,0,0, counts[3, 2] / maxCount))
  polygon(x=c(2, 2, 3, 3), y=c(0, 1, 1, 0), col=rgb(1,0,0, counts[3, 3] / maxCount))
  
  # Add text
  text(x=c(0.5, 1.5, 2.5, 0.5, 1.5, 2.5, 0.5, 1.5, 2.5), 
       y=c(2.5, 2.5, 2.5, 1.5, 1.5, 1.5, 0.5, 0.5, 0.5),
       labels=c(counts[1, 1], counts[1, 2], counts[1, 3],
                counts[2, 1], counts[2, 2], counts[2, 3],
                counts[3, 1], counts[3, 2], counts[3, 3]), cex=3)
  
  # Reset the margins
  par(mar=currentMar)
}

countTransitions <- function(adjacencyMatrix, plotAsMatrix=TRUE){
  
  # Initialise a matrix to store the counts
  counts <- matrix(0, nrow=3, ncol=3)
  colnames(counts) <- c("root", "badger", "cow")
  rownames(counts) <- c("root", "badger", "cow")
  
  # Note the row and column names of the input adjacency matrix
  rowNames <- rownames(adjacencyMatrix)
  colNames <- colnames(adjacencyMatrix)
  
  # Examine each row of the adjacency matrix
  for(row in seq_len(nrow(adjacencyMatrix))){
    
    # Get the row's species index
    rowSpecies <- returnSpeciesIndex(rowNames[row])
    
    # Examine each column of the adjacency matrix
    for(col in seq_len(ncol(adjacencyMatrix))){
      
      # Check if edge present
      if(adjacencyMatrix[row, col] == 1){

        # Get the columns's species index
        colSpecies <- returnSpeciesIndex(colNames[col])

        # Record transition
        counts[rowSpecies, colSpecies] <- counts[rowSpecies, colSpecies] + 1
      }
    }
  }
  
  # Plot the counts if requested
  if(plotAsMatrix){
    plotCountsAsMatrix(counts)
  }else{
    plotCountsAsBars(counts)
  }
  
  return(counts)
}

returnSpeciesIndex <- function(name){
  
  # Assume no index
  index <- -1
  
  # If if root
  if(grepl(name, pattern="ROOT")){
    index <- 1
  # Check if cow
  }else if(grepl(name, pattern="Badger")){
    index <- 2
  }else if(grepl(name, pattern="Cow")){
    index <- 3
  }
  
  return(index)
}

buildSampledTransmissionTree <- function(adjacencyMatrix, sampled){
  
  # Step 1: Recursivelyremove unsampled leaves (out_degree = 0)
  sampledAdjacencyMatrix <- removeUnsampledLeaves(adjacencyMatrix, sampled)
  
  # Step 2: Remove unsampled individuals on path to sampled individuals (in_degree = 1 & out_degree = 1)
  sampledAdjacencyMatrix <- removeUnsampledIndividualsOnPathToSampled(sampledAdjacencyMatrix, sampled)
  
  # Step 3: Recursively remove unnecessary unsampled roots (in_degree = 0 & out_degree = 1)
  sampledAdjacencyMatrix <- removeRootIfUnnecessary(sampledAdjacencyMatrix, sampled)
  
  return(sampledAdjacencyMatrix)
}

removeRootIfUnnecessary <- function(sampledAdjacencyMatrix, sampled, round=1){
  
  # Get the row names of the adjacency matrix
  rowNames <- rownames(sampledAdjacencyMatrix)
  
  # Initialise a counter to record if seed removed
  removed <- FALSE
  
  # Examine each row in the adjacency matrix
  for(row in seq_len(nrow(sampledAdjacencyMatrix))){
    
    # Skip sampled individuals
    if(rowNames[row] %in% sampled){
      next
    }
    
    # Calculate the index of the current individual's source
    source <- which(sampledAdjacencyMatrix[, row] == 1)
    
    # Skip individuals with a source
    if(length(source) != 0){
      next
    }
    
    # Identify the indices of the individuals the current individual infected
    sinks <- which(sampledAdjacencyMatrix[row, ] == 1)
    
    # If current individual has no source (it is the root) and it only infects one individual - remove it
    if(length(sinks) == 1){
      removed <- TRUE
      sampledAdjacencyMatrix[row, sinks] <- 0
    }
  }
  
  # Check if removed any unsampled leaves - if so check for unnecessary root again
  if(removed){
    
    cat(paste0("Removed unnecessary root (round = ", round, ")\n"))
    round <- round + 1
    
    sampledAdjacencyMatrix <- removeRootIfUnnecessary(sampledAdjacencyMatrix, sampled, round)
  }else{
    cat(paste0("Removed no root (round = ", round, ")\n"))
    round <- round + 1
  }
  
  return(sampledAdjacencyMatrix)
}

removeUnsampledIndividualsOnPathToSampled <- function(sampledAdjacencyMatrix, sampled){
  
  # Get the row names of the adjacency matrix
  rowNames <- rownames(sampledAdjacencyMatrix)
  
  # Initialise a counter to record how many individuals removed
  nRemoved <- 0
  
  # Examine each row in the adjacency matrix
  for(row in seq_len(nrow(sampledAdjacencyMatrix))){
    
    # Skip sampled individuals
    if(rowNames[row] %in% sampled){
      next
    }
    
    # Identify the index of the current individual's source
    source <- which(sampledAdjacencyMatrix[, row] == 1)
    
    # Skip individuals already removed - no in degree - no source
    if(length(source) == 0){
      next
    }
    
    # Identify the indices of the individuals the current individual infected
    sinks <- which(sampledAdjacencyMatrix[row, ] == 1)
    
    # Skip individuals that infected more than 1 individual
    if(length(sinks) > 1){
      next
    }
    
    # Remove the edge to the current individual
    nRemoved <- nRemoved + 1
    sampledAdjacencyMatrix[source, row] <- 0
    
    # Remove the edge from the current individual
    sampledAdjacencyMatrix[row, sinks] <- 0
    
    # Join the current individual's source and sink
    sampledAdjacencyMatrix[source, sinks] <- 1
  }
  
  # Note how many leaves were removed
  cat(paste0("Removed ", nRemoved, " unsampled individuals on the path to sampled individuals\n"))

  return(sampledAdjacencyMatrix)
}

removeUnsampledLeaves <- function(sampledAdjacencyMatrix, sampled, round=1){
  
  # Get the row names of the adjacency matrix
  rowNames <- rownames(sampledAdjacencyMatrix)
  
  # Initialise a counter to record how many unsampled leaves found
  nRemoved <- 0
  
  # Examine each row in the adjacency matrix
  for(row in seq_len(nrow(sampledAdjacencyMatrix))){
    
    # Skip sampled individuals
    if(rowNames[row] %in% sampled){
      next
    }
    
    # Calculate the index of the current individual's source
    source <- which(sampledAdjacencyMatrix[, row] == 1)
    
    # Skip individuals already removed - no in degree
    if(length(source) == 0){
      next
    }
    
    # Identify the indices of the individuals the current individual infected
    sinks <- which(sampledAdjacencyMatrix[row, ] == 1)
    
    # Skip individuals that infected 1 or more individuals
    if(length(sinks) > 0){
      next
    }

    # Remove the edge to the current individual - an unsampled leaf
    nRemoved <- nRemoved + 1
    sampledAdjacencyMatrix[source, row] <- 0
  }
  
  # Note how many leaves were removed
  cat(paste0("Removed ", nRemoved, " unsampled leaves (round = ", round, ")\n"))
  round <- round + 1
  
  # Check if removed any unsampled leaves - if so check again
  if(nRemoved > 0){
    sampledAdjacencyMatrix <- removeUnsampledLeaves(sampledAdjacencyMatrix, sampled, round)
  }
  
  return(sampledAdjacencyMatrix)
}

buildAdjacencyMatrix <- function(edgeList){
  
  # Get a list of all the unique
  individuals <- list()
  index <- 0
  ids <- c()
  
  # Initialise the adjacency matrix
  adjacencyMatrix <- matrix(0, nrow=nrow(edgeList) + 1, ncol=nrow(edgeList) + 1)
  
  # Fill the adjacency matrix
  for(row in seq_len(nrow(edgeList))){
    
    # Check the source individual exists
    if(is.null(individuals[[edgeList[row, 1]]]) == TRUE){
      index <- index + 1
      individuals[[edgeList[row, 1]]] <- index
      ids[length(ids) + 1] <- edgeList[row, 1]
    }

    # Check the sink individual exists
    if(is.null(individuals[[edgeList[row, 2]]]) == TRUE){
      index <- index + 1
      individuals[[edgeList[row, 2]]] <- index
      ids[length(ids) + 1] <- edgeList[row, 2]
    }
    
    # Fill in the adjacency matrix: row -> column
    adjacencyMatrix[individuals[[edgeList[row, 1]]], individuals[[edgeList[row, 2]]]] <- 1
  }
  
  # Note the IDs associated with each row and column
  colnames(adjacencyMatrix) <- ids
  rownames(adjacencyMatrix) <- ids
  
  return(adjacencyMatrix)
}

indexIndividuals <- function(edgeList){
  
  # Initialise a list to store each unique ID in edge list
  individuals <- list()
  index <- 0
  
  # Examine the edge list and note each unique ID
  for(row in seq_len(nrow(edgeList))){
    
    # Check the first ID
    if(is.null(individuals[[edgeList[row, 1]]]) == TRUE){
      index <- index + 1
      individuals[[edgeList[row, 1]]] <- index
    }
    
    # Check the second ID
    if(is.null(individuals[[edgeList[row, 2]]]) == TRUE){
      index <- index + 1
      individuals[[edgeList[row, 2]]] <- index
    }
  }
  
  return(individuals)
}

getSampleIDs <- function(xmlFile){
  
  # Open a connection to a file to read (open="r")
  connection <- file(xmlFile, open="r")
  
  # Get all lines from file and store in vector
  fileLines <- readLines(connection)
  
  # Close file connection
  close(connection)
  
  # Initialise a vector to store the IDs
  ids <- c()
  
  # Loop through each of the lines in file
  for(line in fileLines){
    
    # Skip lines that aren't sequence information
    if(grepl(line, pattern="<sequence taxon=\"") == FALSE){
      next
    }
    
    # Get the sequence ID from the current line
    ids[length(ids) + 1] = strsplit(strsplit(line, split="\"")[[1]][2], split="_")[[1]][1]
  }
  
  # Substitute "-" with "_" to match the simulation output
  ids <- gsub("-", "_", ids)
  
  return(ids)
}
