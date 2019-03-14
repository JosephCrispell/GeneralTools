#### Packages ####

library(phangorn)
library(gplots)

#### Read in sequence data ####

# Set the path
path <- "/home/josephcrispell/Desktop/Research/Cumbria/"

# Read in the FASTA file
file <- paste(path, "vcfFiles/sequences_Prox-10_12-03-2019.fasta", sep="")
sequences <- readFasta(file)

# Note the number of sites in the FASTA
nSitesInFasta <- length(sequences[[1]])

# Calculate the proportion Ns for each isolate
propNs <- calculatePropNsForSequences(sequences)

#### Read in the isolate data ####

# Mapping
file <- paste(path, "isolateMappingSummary_12-03-19.txt", sep="")
mapping <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep="\t")
mapping$Prop <- mapping$NumberMappedReads / (mapping$NumberMappedReads + mapping$NumberUnmappedReads)
mapping$Isolate <- parseIds(mapping$Isolate)
isolateMapping <- createList("Isolate", "Prop", mapping)

# Genome coverage
file <- paste(path, "vcfFiles/","IsolateVariantPositionCoverage_RESCUED_12-03-2019.txt", sep="")
vpCoverage <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep="\t")
vpCoverage$Isolate <- parseIds(vpCoverage$Isolate)
isolateVPCoverage <- createList("Isolate", "Coverage", vpCoverage)

# Sampling information
file <- paste(path, "17z_metadata_040319.csv", sep="")
samplingInfo <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep=",")
idLabels <- designTipLabels(samplingInfo)

# Note the IDs of the badgers
badgers <- samplingInfo[grepl(samplingInfo$Notes, pattern="Badger"), "Isolate"]

#### Read in the RAxML tree ####

# Read in the tree
treeFile <- paste0(path, "vcfFiles/mlTree_12-03-19.tree")
tree <- read.tree(treeFile)

# Re-root the tree
tree <- root(tree, 
             outgroup=">Ref-1997", # Set the isolate to root the tree on - can also use node
             resolve.root=TRUE) # Ensures new root will be bifurcating

# Drop reference tip
tree <- drop.tip(tree, ">Ref-1997")

# Change tip labels
tree$tip.label <- parseTipLabels(tree$tip.label)
tips <- tree$tip.label
tree$tip.label <- getValues(idLabels, tree$tip.label)

# Scale the branch lengths to represent SNPs
tree$edge.length <- tree$edge.length * nSitesInFasta

#### Ancestral State Estimation ####

# Get the tip states
tipStates <- ifelse(tips %in% badgers, "badger", "cow")

# Estimate ancestral states under a ER model
ancestralStateFitting <- ace(tipStates, tree, 
                             model="ARD", # "ER" Equal rates; "ARD" All rates different
                             type="discrete") # Discrete states

# Count the estimated transitions
transitionCounts <- countTransitions(tree=tree, tipStates=tipStates,
                                     ancestralStateProbs=ancestralStateFitting$lik.anc,
                                     probThreshold=0.5)

#### Prepare tree for tempest ####

# Get the tip dates in standard format in a tip label
tipLabelsWithDates <- createTipLabelsWithStandardDateFormat(tree$tip.label)

# Create a tree with the dated tip labels
treeWithDates <- tree
treeWithDates$tip.label <- tipLabelsWithDates

# Remove tips without dates
for(tip in treeWithDates$tip.label){
  
  # Check if no date available
  if(grepl(tip, pattern="_NA")){
    
    # Drop the tip
    treeWithDates <- drop.tip(treeWithDates, tip)
  }
}

# Print the tree to file



#### Open an output file for plots ####

# Use the treeFile name as the basis to name the plot file
outputFile <- paste0(substr(treeFile, 1, nchar(treeFile) - 4), "pdf")
pdf(outputFile)

#### Plot the phylogeny ####

# Note the tip size scaling factor
factor <- 1

# Plot the phylogeny
plot.phylo(tree, show.tip.label=TRUE, type="phylogram",
           edge.color="dimgrey", edge.width=3,
           show.node.label=FALSE, label.offset=0.15,
           underscore=TRUE, cex=0.5)

# Add node labels
nodelabels(node=1:length(tree$tip.label), 
           cex=(1 - getValues(propNs, tips)) * factor,
           pch=ifelse(tips %in% badgers, 21, 24),
           bg=ifelse(tips %in% badgers, "red", "blue"), 
           col="dimgrey")

# Get the plotting region dimensions = x1, x2, y1, y2 
# (coordinates of bottom left and top right corners)
dimensions <- par("usr")
xLength <- dimensions[2] - dimensions[1]
yLength <- dimensions[4] - dimensions[3]

# Add quality legend
legend(x=dimensions[2] - (0.1 * xLength), 
       y=dimensions[3] + (0.3 * yLength),
       legend=seq(0.5, 1, 0.1), pch=24, bty='n',
       pt.cex=seq(0.5, 1, 0.1) * factor, xpd=TRUE,
       title="Coverage")

# Add species legend
legend(x=dimensions[1] + (0.4 * xLength),
       y=dimensions[3], 
       legend=c("Cow", "Badger"),
       pch=c(17, 16), cex=1, col=c("blue", "red"), 
       text.col=c("blue", "red"), bty='n', xpd=TRUE)

# Add Scale bar
points(x=c(0, 1), 
       y=c(dimensions[3] - (0.05 * yLength), dimensions[3] - (0.05 * yLength)),
       type="l", lwd=3, xpd=TRUE)
text(x=0.5, y=dimensions[3] - (0.075 * yLength), labels="~ 1 SNP", cex=0.75, xpd=TRUE)

# Add the estimated transition counts
plotTransitionCounts(transitionCounts, dimensions)

#### Close the output file ####

dev.off()



#############
# FUNCTIONS #
#############

createTipLabelsWithStandardDateFormat <- function(tipLabels){
  
  # Initialise an array to store the new tip labels
  output <- c()
  
  # Examine each of the tip labels
  for(i in seq_along(tipLabels)){
    
    # Split the tip label into its parts
    parts <- strsplit(tipLabels[i], split="_")[[1]]
    
    # Skip individuals with no date
    if(parts[3] == "NA"){
      output[i] <- paste0(parts[1], "_NA")
      next
    }
    
    # Get the current tips date
    date <- ""
    if(grepl(parts[3], pattern="/")){
      date <- as.character(as.Date(parts[3], format="%d/%m/%Y"))
    }else{
      date <- as.character(as.Date(paste0("15-", parts[3]), format="%d-%b-%y"))
    }
    
    # Create the new tip label
    output[i] <- paste0(parts[1], "_", date)
  }
  
  return(output)
}

plotTransitionCounts <- function(transitionCounts, dimensions){
  
  # Define the boundaries of the plotting area
  width <- dimensions[2] - dimensions[1]
  height <- dimensions[4] - dimensions[3]
  xCoords <- c(dimensions[1], dimensions[2] - 0.7*width)
  yCoords <- c(dimensions[3] + (0.1*height), dimensions[4] - 0.6*height)
  width <- xCoords[2] - xCoords[1]
  height <- yCoords[2] - yCoords[1]
  
  # Bottom left box
  rect(xleft=xCoords[1], ybottom=yCoords[1], xright=xCoords[1] + (0.5*width), ytop=yCoords[2] - (0.5*height), xpd=TRUE)
  
  # Bottom right box
  rect(xleft=xCoords[1] + (0.5*width), ybottom=yCoords[1], xright=xCoords[1] + width, ytop=yCoords[2] - (0.5*height), xpd=TRUE)
  
  # Top left box
  rect(xleft=xCoords[1], ybottom=yCoords[1] + (0.5*height), xright=xCoords[1] + (0.5*width), ytop=yCoords[2], xpd=TRUE)
  
  # Top right 
  rect(xleft=xCoords[1] + (0.5*width), ybottom=yCoords[1] + (0.5*height), xright=xCoords[1] + width, ytop=yCoords[2], xpd=TRUE)
  
  # Add labels
  text(x=c(xCoords[1] + (0.05*width), xCoords[1] + (0.6*width), xCoords[1] - (0.41*width), xCoords[1] - (0.41*width)), 
       y=c(yCoords[2] + 0.07*height, yCoords[2] + 0.07*height, yCoords[1] + (0.25*height), yCoords[1] + (0.75*height)), 
       labels=c("Badger", "Cow", "Cow", "Badger"), col=c("red", "blue", "blue", "red"), xpd=TRUE, pos=4)
  
  # Add counts
  text(x=c(xCoords[1] + (0.25*width), xCoords[1] + (0.25*width), xCoords[1] + (0.75*width), xCoords[1] + (0.75*width)), 
       y=c(yCoords[1] + (0.75*height), yCoords[1] + (0.25*height), yCoords[1] + (0.75*height), yCoords[1] + (0.25*height)),
       labels=as.vector(transitionCounts))
}

parseTipLabels <- function(tipLabels){
  
  # Initialise a vector to store the output
  output <- c()
  
  # Examine each of the tips
  for(i in seq_along(tipLabels)){
    
    # Remove the ">" from the start
    output[i] <- substr(tipLabels[i], 2, nchar(tipLabels[i]))
    
    # Get the first part
    output[i] <- strsplit(output[i], split="_")[[1]][1]
  }
  
  return(output)
}

plotHeatmap <- function(geneticDistances, order, propNs, badgers){
  
  # Note column and row colours
  colours <- rep("black", nrow(geneticDistances))
  colours[rownames(geneticDistances) %in% badgers] <- "red"
  
  # Get coverage values
  coverage <- getValues(propNs, rownames(geneticDistances))
  
  # Set column and row names
  colnames(geneticDistances) <- getValues(idLabels, colnames(geneticDistances))
  rownames(geneticDistances) <- getValues(idLabels, rownames(geneticDistances))
  
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
                       cexCol=1, # Change the size of the column labels
                       srtCol=90, # Set the angle of the column labels (degrees from horizontal)
                       offsetCol=1, # Set size of space between column labels and heatmap
                       colCol=colours, # Set column label colours
                       
                       # Cell seperating lines
                       colsep=colNumbers, # What columns to place seperator lines between (All)
                       rowsep=rowNumbers, # What rows to place seperator lines between (All)
                       sepwidth=c(0.01, 0.01), # Set size of the seperator lines (colwidth, rowWidth)
                       sepcolor='white', # Set the colour of the seperator line
                       
                       # Change the size of the margins around the plot: c(column space, row space)
                       margins = c(12, 12), 
                       
                       # Row labels
                       cexRow=1, # Change the size of the Row labels
                       srtRow=0,
                       labRow=rownames(geneticDistances), # Create row labels (blank ones in this case)
                       colRow=colours,
                       
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
  
  # Get the plotting window dimensions
  dimensions <- par("usr")
  xLength <- dimensions[2] - dimensions[1]
  yLength <- dimensions[4] - dimensions[3]
  
  # Add title
  text(x=dimensions[1] + (0.4 * xLength), 
       y=dimensions[4] + (0.1 * yLength), "Minimum Genetic distance", cex=1, xpd=TRUE)
  
  # Add coverage values down side
  rbPal <- colorRampPalette(c("blue","red"))
  colours <- rbPal(10)[as.numeric(cut(1-coverage,breaks = 10))]
  indices <- rev(heatmap$rowInd)
  text(x=rep(dimensions[1] - (0.06 * xLength), nrow(geneticDistances)),
       y=seq(from=dimensions[4] + (0.01 * yLength), 
             to=dimensions[3] + (0.21 * yLength),
             length.out=nrow(geneticDistances)),
       labels=round(1-coverage, digits=2)[indices], xpd=TRUE, cex=0.75,
       col=colours[indices])
  text(x=dimensions[1] - (0.06 * xLength), y=dimensions[4] + (0.05 * yLength), "Coverage", cex=0.75, xpd=TRUE)
}

calculateRowSums <- function(table){
  
  rowSums <- c()
  for(row in 1:nrow(table)){
    rowSums[row] <- sum(table[row, ], na.rm=TRUE)
  }
  
  return(rowSums)
}

printOutInformativeSequences <- function(informativeSequences, fileName, newLabels){
  fileConnection <- file(fileName, open="w")
  writeLines(paste(length(informativeSequences), "\t", length(informativeSequences[[1]])), fileConnection)
  keys <- names(informativeSequences)
  for(i in 1:length(keys)){
    writeLines(paste(newLabels[i], "\t", paste(informativeSequences[[keys[i]]],
                                               collapse="")), fileConnection)
  }
  close(fileConnection)
}

keepSites <- function(sequences, informativeSites){
  
  output <- list()
  
  for(key in names(sequences)){
    output[[key]] <- sequences[[key]][informativeSites]
  }
  
  return(output)
}

noteInformativeSites <- function(sequences){
  keep <- c()
  ids <- names(sequences)
  
  for(pos in 1:length(sequences[[1]])){
    
    for(i in 1:length(ids)){
      
      for(j in 1:length(ids)){
        
        if(i >= j){
          next
        }
        
        if(ids[i] != "Ref-1997" && ids[j] != "Ref-1997" &&
           sequences[[ids[i]]][pos] != "N" && sequences[[ids[j]]][pos] != "N" && 
           sequences[[ids[i]]][pos] != sequences[[ids[j]]][pos]){
          keep[length(keep) + 1] <- pos
        }
      }
    }
  }
  
  return(unique(keep))
}

designTipLabels <- function(samplingInfo){
  
  tipLabels <- list()
  for(row in 1:nrow(samplingInfo)){
    
    id <- strsplit(samplingInfo[row, "Isolate"], split="-")[[1]][3]
    
    tipLabels[[samplingInfo[row, "Isolate"]]] <- paste(id, samplingInfo[row, "Location"], 
                                                       samplingInfo[row, "Date"], sep="_")
  }
  
  return(tipLabels)
}

getValues <- function(list, keys){
  
  values <- c()
  for(i in 1:length(keys)){
    if(is.null(list[[keys[i]]]) == FALSE){
      values[i] <- list[[keys[i]]]
    }else{
      values[i] <- keys[i]
    }
  }
  
  return(values)
}

createList <- function(keysColumn, valuesColumn, table){
  
  output <- list()
  for(row in 1:nrow(table)){
    output[[table[row, keysColumn]]] <- table[row, valuesColumn]
  }
  
  return(output)
}

parseIds <- function(strings){
  output <- c()
  for(i in 1:length(strings)){
    
    output[i] <- strsplit(strings[i], split="_")[[1]][1]
  }
  
  return(output)
}

removeBlankIsolates <- function(geneticDistances){
  
  output <- geneticDistances
  
  ignore <- c()
  for(i in 1:nrow(output)){
    
    if(sum(output[i, ], na.rm=TRUE) == 0){
      ignore[length(ignore) + 1] <- i
    }
  }
  
  if(length(ignore) > 0){
    output <- output[-ignore, -ignore]
  }
  
  return(output)
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

calculatePropNsForSequences <- function(sequences){
  
  propNs <- list()
  for(key in names(sequences)){
    
    propNs[[key]] <- calculatePropNs(sequences[[key]])
  }
  
  return(propNs)
}

calculatePropNs <- function(sequence){
  
  propNs <- 0
  for(i in 1:length(sequence)){
    
    if(sequence[i] == "N"){
      propNs <- propNs + 1
    }
  }
  
  return(propNs / length(sequence))
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

getNSitesInFASTA <- function(fastaFile){
  
  # Open a connection to a file to read (open="r")
  connection <- file(fastaFile, open="r")
  
  # Get first line of file
  firstLine <- readLines(connection, n=1)
  
  # Close file connection
  close(connection)
  
  # Get the number of sites used in the FASTA file from first line
  nSites <- as.numeric(strsplit(firstLine, " ")[[1]][2])
  
  return(nSites)
}

countTransitions <- function(tree, tipStates, ancestralStateProbs, probThreshold=0.5){
  
  # A more conservative transition count method - recommended by Nicola De Maio
  # 	- Assumes the ancestor represents one of the tips in the past
  # 
  # For a two state problem: A & B
  #  ----A				----A				----A				  ----B	
  # |		  1 AA	 |		0 AA	 |		  1 AA	 |		 0 AA
  # A		  0 AB	 A		1 AB	 A----B	1 AB	 A		 2 AB
  # |				     |				   |				     |
  # ----A			   ----B			 ----A				 ----B
  
  # Nodes in phylogenetic tree are numbered:
  # 1:Ntip(tree) (1:Ntip(tree) + tree$Nnode)
  # 
  # Edges on the phylogenetic tree are under tree$edge
  # A matrix of two columns: FROM, TO - based upon the node indices defined above
  #
  # Ancestral state probs - probability of each state for each INTERNAL node
  #
  # Probability threshold is the threshold to define the state assigned to each internal node
  
  # Note the state probabilities for the tips
  tipStateProbs <- matrix(0, nrow=length(tipStates), ncol=2)
  tipStateProbs[tipStates == "cow", 2] <- 1
  tipStateProbs[tipStates == "badger", 1] <- 1
  
  # Add the state probabilities for the tips onto the ancestral state probs for ease
  stateProbs <- rbind(tipStateProbs, ancestralStateProbs)
  
  # Initialise a list to store the states associated with each internal nodes daughter nodes
  nodes <- list()

  # Examine each edge in the phylogeny
  for(row in 1:nrow(tree$edge)){
    
    # Note the node index of from node
    fromNode <- as.character(tree$edge[row, 1])
    
    # Get the state for the from node
    fromState <- which(stateProbs[as.numeric(fromNode), ] > probThreshold)
    if(length(fromState) == 0){
      fromState <- -1
    }
    
    # Get the state for the to node
    toState <- which(stateProbs[tree$edge[row, 2], ] > probThreshold)
    if(length(toState) == 0){
      toState <- -1
    }
    
    # Check if we have encountered this from node before
    if(is.null(nodes[[fromNode]]) == FALSE){
      
      nodes[[fromNode]]$Daughters <- c(nodes[[fromNode]]$Daughters, toState)
      
      # Create a record for the current node if we haven't encountered it before
    }else{
      nodes[[fromNode]] <- list("Parent"=fromState, "Daughters"=c(toState))
    }
  }
  
  # Initialise a matrix to count the transitions
  transitionCounts <- matrix(0, nrow=2, ncol=2)
  colnames(transitionCounts) <- colnames(stateProbs)
  rownames(transitionCounts) <- colnames(stateProbs)
  
  # Examine each of the node
  for(node in names(nodes)){
    
    # Initialise a variable to count the number of branches beginning and ending with the same state
    nSameStateBranches <- 0
    
    # Skip node if no state available
    if(nodes[[node]][1] == -1){
      next
    }
    
    # Get the state of the parent
    parent <- nodes[[node]]$Parent
    
    # Get the states of the daughters for the current node
    daughters <- nodes[[node]]$Daughters
    
    # Examine each of the daughter states
    for(i in seq_along(daughters)){
      
      # Check if state of daughter is the same as that of parent
      if(daughters[i] != -1 && daughters[i] == parent){
        nSameStateBranches <- nSameStateBranches + 1
        
        # Found branch beginning and ending with different states - count the inter species transmission event
      }else if(daughters[i] != -1){
        transitionCounts[parent, daughters[i]] <- transitionCounts[parent, daughters[i]] + 1
      }
    }
    
    # Count the number of within species transmission events - conservatively assumes parent node is same animal as one of daughters
    if(nSameStateBranches > 1){
      for(i in seq_len(nSameStateBranches - 1)){
        transitionCounts[parent, parent] <- transitionCounts[parent, parent] + 1
      }
    }
  }
  
  return(transitionCounts)
}
