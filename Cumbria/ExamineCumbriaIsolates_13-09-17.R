############
# Packages #
############

library(phangorn)
library(gplots)

#########################
# Read in sequence data #
#########################

# Set the path
path <- "/home/josephcrispell/Desktop/Research/Cumbria/"

# Read in the FASTA file
file <- paste(path, "vcfFiles/sequences_Prox-10_07-03-2019.fasta", sep="")
sequences <- readFasta(file)

# Note the number of sites in the FASTA
nSitesInFasta <- length(sequences[[1]])

############################
# Read in the isolate data #
############################

# Mapping
file <- paste(path, "isolateMappingSummary_07-03-19.txt", sep="")
mapping <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep="\t")
mapping$Prop <- mapping$NumberMappedReads / (mapping$NumberMappedReads + mapping$NumberUnmappedReads)
mapping$Isolate <- parseIds(mapping$Isolate)
isolateMapping <- createList("Isolate", "Prop", mapping)

# Genome coverage
file <- paste(path, "vcfFiles/","IsolateVariantPositionCoverage_RESCUED_07-03-2019.txt", sep="")
vpCoverage <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep="\t")
vpCoverage$Isolate <- parseIds(vpCoverage$Isolate)
isolateVPCoverage <- createList("Isolate", "Coverage", vpCoverage)

# Sampling information
file <- paste(path, "17z_metadata_040319.csv", sep="")
samplingInfo <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep=",")
idLabels <- designTipLabels(samplingInfo)

###################################
# Calculate the genetic distances #
###################################

# Set the margins
par(mfrow=c(1,1))
par(mar=c(5,1,1,2)) # Bottom, Left, Top, Right

# Calculate the proportion Ns for each isolate
propNs <- calculatePropNsForSequences(sequences)

# Build genetic distance matrix
geneticDistances <- buildGeneticDistanceMatrix(sequences)

# Drop referencce and isolates that are only different from the reference
geneticDistances <- removeBlankIsolates(geneticDistances)

##########################
# Read in the RAxML tree #
##########################

# Read in the tree
file <- paste0(path, "vcfFiles/mlTree_07-03-19.tree")
tree <- read.tree(file)

# Change tip labels
tree$tip.label <- parseTipLabels(tree$tip.label)
tips <- tree$tip.label
tree$tip.label <- getValues(idLabels, tree$tip.label)

# Scale the branch lengths to represent SNPs
tree$edge.length <- tree$edge.length * nSitesInFasta

# Plot tree with and without reference
file <- paste(path, "vcfFiles/CumbrianIsolates_SummaryPlots_07-03-19.pdf", sep="")
pdf(file)

#### WITH REFERENCE

badgers <- samplingInfo[grepl(samplingInfo$Notes, pattern="Badger"), "Isolate"]
factor <- 1
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


##### WITHOUT REFERENCE

# Drop reference tip
tree <- drop.tip(tree, "Ref-1997")

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


#############
# SNP table #
#############

# # Build SNP table for sites where there are two isolates with coverage that are different
# informativeSites <- noteInformativeSites(sequences)
# 
# # Build smaller sequences using only the informative sites
# informativeSequences <- keepSites(sequences, informativeSites)
# 
# # Print out SNP table
# newLabels <- getValues(idLabels, names(informativeSequences))
# file <- paste(path, "vcfFiles/InformativeFasta_07-03-19.fasta", sep="")
# printOutInformativeSequences(informativeSequences, file, newLabels)

#################################
# Min Genetic Distance Heatmap? #
#################################
# 
# # Build genetic distance matrix
# geneticDistances <- buildGeneticDistanceMatrix(sequences)
# 
# # Plot ordered heatmap with Reference
# plotHeatmap(geneticDistances, order=TRUE, propNs, badgers)
# 
# # Plot ordered heatmap without Reference
# refRow <- which(rownames(geneticDistances) == "Ref-1997")
# plotHeatmap(geneticDistances[-refRow, -refRow], order=TRUE, propNs, badgers)

# Close the PDF file (must be closed when code is running as well)
dev.off()



#############
# FUNCTIONS #
#############

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
