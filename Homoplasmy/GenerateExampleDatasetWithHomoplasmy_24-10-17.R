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

# Set the seed
#seed <- sample(x=1:10000000, size=1)
#print(seed)
#set.seed(seed)
#set.seed(5962222)

# Generate the sequences
sequences <- simulateSequences(n=25, genomeSize=10000, mutationRate=1)

# Note the source individual
source <- sequences[["source"]]
sequences[["source"]] <- NULL
cat(paste("Source = ", source, "\n", sep=""))

####################
# Insert homoplasy #
####################

output <- insertHomoplasy(sequences, source)
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

##################
# Plot a heatmap #
##################

plotHeatmap(buildGeneticDistanceMatrix(sequences), order=TRUE)

#############
# FUNCTIONS #
#############

plotHeatmap <- function(geneticDistances, order){
  
  par(mar=c(2, 0, 0, 2))
  
  # Create vectors to store the Row and Column numbers
  colNumbers <- seq(from=1, to=ncol(geneticDistances))
  rowNumbers <- seq(from=1, to=nrow(geneticDistances))
  
  # Define the colour breaks
  colBreaks <- seq(0, max(geneticDistances, na.rm=TRUE), by=0.1)
  
  # Plot the heatmap
  heatmap <- heatmap.2(geneticDistances, # matrix is the input data
                       
                       # Add cell labels
                       cellnote=geneticDistances,
                       notecol="black",
                       
                       # Create the colour scale 
                       col=colorpanel(n=length(colBreaks)-1, low="yellow", high="red"),
                       breaks=colBreaks,
                       
                       # Turn off a density plot
                       density.info="none", 
                       
                       # Turn off the trace
                       trace="none",
                       
                       # Column Labels
                       cexCol=2, # Change the size of the column labels
                       srtCol=0, # Set the angle of the column labels (degrees from horizontal)
                       offsetCol=1, # Set size of space between column labels and heatmap

                       # Cell seperating lines
                       colsep=colNumbers, # What columns to place seperator lines between (All)
                       rowsep=rowNumbers, # What rows to place seperator lines between (All)
                       sepwidth=c(0.01, 0.01), # Set size of the seperator lines (colwidth, rowWidth)
                       sepcolor='white', # Set the colour of the seperator line
                       
                       # Change the size of the margins around the plot: c(column space, row space)
                       margins = c(5, 5), 
                       
                       # Row labels
                       cexRow=2, # Change the size of the Row labels
                       srtRow=0,
                       labRow=rownames(geneticDistances), # Create row labels (blank ones in this case)

                       # Make sure the order of the rows and columns is changed
                       Rowv=order, Colv=order,
                       
                       # Make sure all NA values are kept in        
                       na.rm=TRUE,
                       
                       # Don't plot any dendogram
                       dendrogram="none",
                       
                       # Set up the Key
                       key=FALSE, # Turn the key on
                       
                       # Say where to plot each part of the heatmap
                       #     4     3
                       #     2     1
                       # 1. Heatmap
                       # 2. Row Dendrogram
                       # 3. Column Dendrogram
                       # 4. Key
                       lmat=rbind(4:3, 2:1),
                       
                       # Set the size of the spaces for output plots:
                       #     Colour Key      |   Column Dendrogram
                       #     -------------------------------------
                       #     Row Dendrogram  |   Heatmap
                       #
                       # Note that these will be affected by the width and height you set for the PDF
                       lhei=c(1, 10), # c(row1Width, row2Width)
                       lwid=c(1, 10), # c(column1Width, column2Width)
                       
                       # Note that the input matrix is not symmetric
                       symm = FALSE
  )
}

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
             tip.color=ifelse(tree$tip.label %in% ij, "red", "black"),
             cex=1.5)
  
  # Add scale
  lines(x=c(2,3), y=c(1.5, 1.5), lwd=4)
  text(x=2.5, y=1.1, labels="1 SNP")
  
  # Write the tree to file
  write.tree(tree, file=file, append=FALSE,
             digits=20, tree.names=FALSE)
}

insertHomoplasy <- function(sequences, source){
  
  # Get the isolates IDs
  isolates <- names(sequences)

  # Initialise an array to store the positions to choose from
  positionsToChooseFrom <- c()
  
  # Initialise variables to store the source and sink individuals
  i <- NULL
  j <- NULL
  
  # Search for suitable source and sink individuals
  while(length(positionsToChooseFrom) == 0){
    # Select a source and sink
    i <- sample(size=1, x=isolates[isolates != source])
    j <- sample(size=1, x=isolates[isolates != i & isolates != source])
    
    # Note the positions to choose from - not the same or reference
    positionsToChooseFrom <- 
      which(sequences[[i]] != sequences[[j]] &
            sequences[[i]] != sequences[[source]])
  }
    
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

simulateSequences <- function(n, genomeSize, mutationRate){
  # Note the nucleotides
  nucleotides <- c("A", "T", "G", "C")
  
  # Set population size
  susceptibles <- 1:n
  
  # Choose random seed
  sources <- c(sample(size=1, susceptibles))
  susceptibles <- susceptibles[susceptibles != sources[1]]
  
  # Store each individual's sequence
  sequences <- c()
  sequences[[as.character(sources[1])]] <- sample(x=nucleotides, size=genomeSize,
                                                  replace=TRUE)
  
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
    nMutations <- rpois(n=1, lambda=mutationRate)
    positions <- sample(size=nMutations, x=1:length(sequence))
    for(position in positions){
      sequence[position] <- sample(size=1, x=nucleotides[nucleotides != sequence[position]])
    }

    # Transfer sequence
    sequences[[as.character(sink)]] <- sequence
    
    cat(paste(source, " -> ", sink, "\n"))
  }
  
  # Remove uninformative sites
  sequences <- removeUninformativeSites(sequences)
  
  # Record the source
  sequences[["source"]] <- sources[1]
  
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