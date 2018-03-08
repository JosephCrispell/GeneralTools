###############
# Preparation #
###############

# Packages
library(phangorn)
library(gplots)
library(geiger) # For the tips function

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Homoplasmy/"

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

###########
# Testing #
###########

###################
# Run simulations #
###################

n <- 1000
results <- data.frame(PropFound=rep(NA, n), NIncorrect=rep(NA, n),
                      Seed=rep(NA, n))

for(i in 1:n){

  cat(paste("\rRunning simulation", i, "of", n))

  results[i, 3] <- sample(1:100000000, 1)
  set.seed(results[i, 3])

  # Generate the sequences
  simulationOutput <- runSimulation(n=300, mutationRate=0.5, infectiousness=0.001,
                                    samplingProb=0.05, nToSample=150, FALSE)

  # Build the sequences based upon the mutation events
  sequences <- buildSequences(simulationOutput)

  # Insert homoplasies
  homoplasyInsertionInfo <- insertHomoplasies(sequences, n=10, FALSE)
  sequences <- homoplasyInsertionInfo[["sequences"]]

  # Build FASTA
  sequences[["REF"]] <- NULL
  writeFasta(sequences, paste(path, "example_", date, ".fasta", sep=""))

  # Build phylogeny
  buildPhylogeny(sequences, paste(path, "example_", date, ".tree", sep=""),
                 homoplasyInsertionInfo, FALSE)

  # Run HomoplasyFinder
  setwd(path)
  fastaFile <- paste("example_", date, ".fasta", sep="")
  treeFile <- paste("example_", date, ".tree", sep="")
  system(paste("java -jar HomoplasyFinder_06-03-18.jar", 0, fastaFile, treeFile, sep=" "),
         ignore.stdout=FALSE)

  # Get and Check HomoplasyFinder output
  file <- paste(path, "homoplasyReport_", date, ".txt", sep="")
  homoplasyFinderOutput <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                                      check.names=FALSE)
  results[i, c(1,2)] <- reportHowManyHomoplasiesWereFound(homoplasyFinderOutput, homoplasyInsertionInfo)

  propFound <- results[i, 1]
}

file <- paste(path, "TestingHomoplasyFinder_", n, "-", mutationRate,
              "-", infectiousness, "-", samplingProb, "-", nToSample, "_",
              nHomoplasies, "_", nrow(results), "_", date, ".pdf", sep="")
pdf(file)

plot(x=results$NIncorrect, y=results$PropFound, ylim=c(0,1),
     col=rgb(0,0,0, 0.01), bg=rgb(0,0,0, 0.01), cex=4, las=1, pch=21, bty="n",
     ylab=paste("Proportion of homoplasies (n=", 10, ") found", sep=""),
     xlab="Number of homoplasies incorrectly identified", xpd=TRUE,
     main="Testing HomoplasyFinder")
legend("topright", legend=c(
  "Mean values:",
  paste("Proportion found =", round(mean(results$PropFound), digits=1)),
  paste("Number incorrect =", round(mean(results$NIncorrect), digits=1))), bty="n")

# par(mar=c(0.1, 0.1, 4, 0.1))
# plot(results$Seed, xaxt="n", yaxt="n", bty="n", xlab="", ylab="", main="Seeds used in simulations",
#      pch=19, col=rgb(0,0,0, 0.5))
# par(mar=c(5.1, 4.1, 4.1, 2.1))

dev.off()

#############
# FUNCTIONS #
#############

reportHowManyHomoplasiesWereFound <- function(homoplasyFinderOutput, homoplasyInsertionInfo){
  
  # Initialise a variable to record how many of the inserted homoplasies were found
  nFound <- 0
  
  # Examine each of the homoplasies inserted
  for(key in names(homoplasyInsertionInfo)[names(homoplasyInsertionInfo) != "sequences"]){
    
    # Search for current homoplasy in HomoplasyFinder output
    rows <- which(homoplasyFinderOutput$Position == homoplasyInsertionInfo[[key]]$position)
    if(length(row) == 0){
      next
    }
    
    for(row in rows){
      
      # Find the allele index
      alleleIndex <- which(strsplit(homoplasyFinderOutput[row, "Alleles"], split=",")[[1]] == 
                             homoplasyInsertionInfo[[key]]$allele)
      if(length(alleleIndex) == 0){
        next
      }
      
      # Check all the isolates associated with the source an dsink nodes are present
      present <- TRUE
      sourceSinkIsolates <- c(homoplasyInsertionInfo[[key]]$sourceIsolates,
                              homoplasyInsertionInfo[[key]]$sinkIsolates)
      isolates = strsplit(
        strsplit(homoplasyFinderOutput[row, "IsolatesForAlleles"], split=",")[[1]][alleleIndex], 
        split="-")[[1]]
      for(isolate in sourceSinkIsolates){
        
        if(isolate %in% isolates == FALSE){
          present <- FALSE
          break
        }
      }
      if(present == TRUE){
        nFound <- nFound + 1
      }
    }
  }
  
  # Calculate proportion found and number of false positives (sequences also present in insertion info)
  propFound <- nFound / (length(homoplasyInsertionInfo) - 1)
  
  # Calculate proportion found that weren't inserted
  nWronglyFound <- nrow(homoplasyFinderOutput) - nFound
  if(nWronglyFound < 0){
    nWronglyFound <- 0 # If a homoplasy is inserted at the same site then this number could be negative
  }
  
  return(c(propFound, nWronglyFound))
}

buildPhylogeny <- function(sequences, file, homoplasyInsertionInfo, verbose){
  
  # Build genetic distance matrix
  geneticDistances <- buildGeneticDistanceMatrix(sequences)
  
  # Build neighbour joining tree
  tree <- NJ(geneticDistances)

  # Plot tree if wanting verbose information
  if(verbose){
    
    par(mar=c(0,0,0,0))
    
    # Define the colour of the nodes based on the inserted homoplasies
    colours <- c("red", "blue", "green", "cyan", "orange", "darkorchid4",
                 "deeppink", "black", "brown", "darkolivegreen")
    nodeColours <- rep(rgb(0,0,0, 0), length(tree$edge.length) + 1)
    for(homoplasyId in names(homoplasyInsertionInfo)[names(homoplasyInsertionInfo) != "sequences"]){
      
      # Get the colour for the current homoplasy
      colour <- colours[as.numeric(homoplasyId) %% length(colours)]
      if(length(colour) == 0){
        colour <- colours[length(colours)]
      }
      
      # Assign that colour to the source and sink nodes for the current homoplasy
      nodeColours[as.numeric(homoplasyInsertionInfo[[homoplasyId]]$sourceNode)] <- colour
      nodeColours[as.numeric(homoplasyInsertionInfo[[homoplasyId]]$sinkNode)] <- colour
    }
    
    # Plot the tree
    plot.phylo(tree, show.tip.label=TRUE, type="phylogram",
               edge.color="grey", edge.width=2,
               show.node.label=TRUE, label.offset=0.15,
               tip.color="black", cex=1)
    
    # Add node circles to highlight homoplasies
    nodelabels(node=1:(length(tree$edge.length) + 1), pch=19, col=nodeColours)
    
    # Get the axis limits
    axisLimits <- par("usr")
    
    # Add scale
    xPos <- axisLimits[2] - (0.1 * (axisLimits[2] - axisLimits[1]))
    yPos <- axisLimits[3] + (0.075 * (axisLimits[4] - axisLimits[3]))
    lines(x=c(xPos,xPos+1), y=c(yPos, yPos), lwd=4)
    text(x=xPos+0.5, y=axisLimits[3], labels="1 SNP", pos=3)
    
    par(mar=c(5.1, 4.1, 4.1, 2.1))
  }
  
  # Write the tree to file
  write.tree(tree, file=file, append=FALSE,
             digits=20, tree.names=FALSE)
}

insertHomoplasies <- function(n, sequences, verbose){
  
  # Build genetic distance matrix
  geneticDistances <- buildGeneticDistanceMatrix(sequences)
  
  # Build neighbour joining tree
  tree <- NJ(geneticDistances)
  if(verbose == TRUE){
    par(mar=c(0,0,0,0))
    plot.phylo(tree, show.tip.label=TRUE, type="phylogram",
               edge.color="grey", edge.width=2,
               show.node.label=TRUE, label.offset=0.15,
               tip.color="black", cex=1)
    par(mar=c(5.1, 4.1, 4.1, 2.1))
  }
  
  # Get the isolates associated with each node in the phylogeny
  nodes <- getNodes(tree)
  
  # Initialise a list to store the homoplasy information
  output <- list()
  
  # Insert the homoplasies
  for(i in 1:n){
    
    source <- NULL
    sink <- NULL
    sitesUniqueToIsolates <- c()
    used <- c()
    while(is.null(source) == TRUE || is.null(sink) == TRUE || length(sitesUniqueToIsolates) == 0){
      # Randomly select a node to act as source
      source <- sample(names(nodes), size=1)
      
      # Randomly select a node to act as a sink - as long as it isn't on the path to the root
      sink <- randomlySelectSinkNode(nodes, source)
      
      # Check that haven't already used this pair
      if(paste(source, sink, sep=":") %in% used == TRUE){
        next
      }else{
        used[length(used) + 1] <- paste(source, sink, sep=":")
      }
      
      # Find the sites that are common to all the isolates under the current node
      sitesUniqueToIsolates <- findConservedSitesThatAreUniqueToClade(nodes[[source]], sequences)
    }
    
    # Randomly choose a site
    site <- sample(sitesUniqueToIsolates, size=1)
    allele <- sequences[[nodes[[source]][1]]][site]
    
    # Assign the randomly chosen site and allele to those isolates
    sequences <- assignSiteAndAllele(nodes[[sink]], sequences, site, allele)
    
    # Record homoplasy
    if(verbose){
      cat(paste("---------------------------------------------------------------------------\n",
                "Introduced homoplasy: From node: ", source, " to ", sink, 
                "\tPosition: ", site, "\tAllele: ", allele,
                "\nSource Isolates: ", paste(nodes[[source]], collapse=", "),
                "\nSink Isolates: ", paste(nodes[[sink]], collapse=", "), "\n", sep=""))
    }
    output[[as.character(i)]] <- list(
      "sourceNode" = source,
      "sinkNode" = sink,
      "position" = site,
      "allele" = allele,
      "sourceIsolates"=nodes[[source]],
      "sinkIsolates"=nodes[[sink]])
  }
  
  # Add the sequences to the output
  output[["sequences"]] <- sequences
  
  return(output)
}

randomlySelectSinkNode <- function(nodes, source){
  
  # Note all the nodes that don't contain any of the isolates from the source node
  nodesToSelectFrom <- c()
  for(node in names(nodes)[names(nodes) != source]){
    
    use <- TRUE
    for(isolate in nodes[[node]]){
      if(isolate %in% nodes[[source]]){
        use <- FALSE
        break
      }
    }
    
    if(use == TRUE){
      nodesToSelectFrom[length(nodesToSelectFrom) + 1] <- node
    }
  }
  
  selectedNode <- NULL
  if(length(nodesToSelectFrom) > 0){
    selectedNode <- sample(nodesToSelectFrom, size=1)
  }
  
  return(selectedNode)
}

assignSiteAndAllele <- function(isolates, sequences, site, allele){
  for(isolate in names(sequences)){
    
    # Skip isolates that weren't under node
    if(isolate %in% isolates == FALSE){
      next
    }
    
    # Assign the site and and allele to the current isolate
    sequence <- sequences[[isolate]]
    sequence[site] <- allele
    sequences[[isolate]] <- sequence
  }
  
  return(sequences)
}

findConservedSitesThatAreUniqueToClade <- function(isolates, sequences){
  
  # Note the conserved sites in the current isolates
  conserved <- rep(TRUE, length(sequences[[1]]))
  for(pos in 1:length(sequences[[1]])){
    
    # Get the allele from the first isolate
    allele <- sequences[[isolates[1]]][pos]
    
    # Compare that allele to all others - finish if find a different one
    for(isolate in isolates[-1]){
      
      if(allele != sequences[[isolate]][pos]){
        conserved[pos] <- FALSE
        break
      }
    }
  }
  
  # Find conserved sites that are unique to isolates
  uniqueToIsolates <- conserved
  for(pos in 1:length(sequences[[1]])){
    
    # Skip if not conserved in isolates
    if(conserved[pos] == FALSE){
      next
    }
    
    # Check if you can find the isolates allele in any of the other isolates
    allele <- sequences[[isolates[1]]][pos]
    for(id in names(sequences)){
      if(id %in% isolates){
        next
      }
      
      if(allele == sequences[[id]][pos]){
        uniqueToIsolates[pos] <- FALSE
        break
      }
    }
  }
  
  return(c(1:length(sequences[[1]]))[uniqueToIsolates])
}

getNodes <- function(tree){
  nodes <- list()
  for(node in 1:(length(tree$edge.length) + 1)){
    nodes[[as.character(node)]] <- tips(tree, node)
  }
  
  return(nodes)
}

reportHowManyHomoplasiesWereFoundOld <- function(homoplasyFinderOutput, homoplasyInsertionInfo){
  
  # Initialise a variable to record how many of the inserted homoplasies were found
  nFound <- 0
  
  # Examine each of the homoplasies inserted
  for(key in names(homoplasyInsertionInfo)){
    
    # Ignore the sequences
    if(key == "sequences"){
      next
    }
    
    # Search for current homoplasy in HomoplasyFinder output
    row <- which(homoplasyFinderOutput$Position == homoplasyInsertionInfo[[key]]$position)
    if(length(row) == 0){
      next
    }
    
    
    # Find the allele index
    alleleIndex <- which(strsplit(homoplasyFinderOutput[row, "Alleles"], split=",")[[1]] == 
                           homoplasyInsertionInfo[[key]]$allele)
    if(length(alleleIndex) == 0){
      next
    }
    
    # Check if source and sink are present amongst isolates
    isolates = strsplit(
      strsplit(homoplasyFinderOutput[row, "IsolatesForAlleles"], split=",")[[1]][alleleIndex], 
      split="-")[[1]]
    if(homoplasyInsertionInfo[[key]]$source %in% isolates && 
       homoplasyInsertionInfo[[key]]$sink %in% isolates){
      nFound <- nFound + 1
    }
  }
  
  # Calculate proportion found and number of false positives
  propFound <- nFound / (length(homoplasyInsertionInfo) - 1)
  
  # Calculate proportion found that weren't inserted
  nWronglyFound <- nrow(homoplasyFinderOutput) - nFound
  if(nWronglyFound < 0){
    nWronglyFound <- 0 # If a homoplasy is inserted at the same site then this number could be negative
  }
  
  return(c(propFound, nWronglyFound))
}

run <- function(path, date, n, mutationRate, infectiousness,
                samplingProb, nToSample, nHomoplasies, verbose){
  
  ######################
  # Generate sequences #
  ######################
  
  # Generate the sequences
  simulationOutput <- runSimulation(n, mutationRate, infectiousness,
                                    samplingProb, nToSample, verbose)
  
  # Build the sequences based upon the mutation events
  sequences <- buildSequences(simulationOutput)
  
  ####################
  # Insert homoplasy #
  ####################
  
  homoplasyInsertionInfo <- insertHomoplasies(sequences, n=nHomoplasies, verbose)
  sequences <- homoplasyInsertionInfo[["sequences"]]
  
  ###############
  # Build FASTA #
  ###############
  sequences[["REF"]] <- NULL
  writeFasta(sequences, paste(path, "example_", date, ".fasta", sep=""))
  
  ###################
  # Build phylogeny #
  ###################
  
  buildPhylogeny(sequences, paste(path, "example_", date, ".tree", sep=""),
                 homoplasyInsertionInfo, rootOnREF=FALSE, verbose)
  
  ################################
  # Run the HomoplasyFinder tool #
  ################################
  
  setwd(path)
  fastaFile <- paste("example_", date, ".fasta", sep="")
  treeFile <- paste("example_", date, ".tree", sep="")
  system(paste("java -jar HomoplasyFinder_06-03-18.jar", 0, fastaFile, treeFile, sep=" "),
         ignore.stdout=FALSE)
  
  ###########################################
  # Read in the HomoplasyFinder output tool #
  ###########################################
  
  file <- paste(path, "homoplasyReport_", date, ".txt", sep="")
  homoplasyFinderOutput <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                                      check.names=FALSE)
  
  #########################################
  # Calculate sensitivity and specificity #
  #########################################
  
  results <- reportHowManyHomoplasiesWereFound(homoplasyFinderOutput, homoplasyInsertionInfo)
  
  return(results)
}

checkWhichHomoplasiesWereFound <- function(homoplasyFinderOutput, homoplasyInsertionInfo, verbose){
  
  # Get the isolates associated with the homoplasis identified
  isolatesOfHomoplasiesIdentified <- getIsolatesFromHomoplasyInfo(homoplasyFinderOutput)
  
  # Check  if each of the inserted homoplasies were found
  inserted <- names(homoplasyInsertionInfo)
  inserted <- inserted[inserted != "sequences"]
  nFound <- 0
  for(key in inserted){
    
    source <- homoplasyInsertionInfo[[key]]$source
    sink <- homoplasyInsertionInfo[[key]]$sink
    position <- homoplasyInsertionInfo[[key]]$position
    
    for(i in 1:length(isolatesOfHomoplasiesIdentified)){
      
      if(position == isolatesOfHomoplasiesIdentified[[i]]$position &&
         (source %in% isolatesOfHomoplasiesIdentified[[i]]$foundIn == TRUE ||
          source %in% isolatesOfHomoplasiesIdentified[[i]]$from == TRUE) && 
         (sink %in% isolatesOfHomoplasiesIdentified[[i]]$foundIn == TRUE ||
          sink %in% isolatesOfHomoplasiesIdentified[[i]]$from == TRUE)){
        nFound <- nFound + 1
        break
      }
    }
  }
  
  # Calculate proportion homoplasies found
  propFound <- nFound / length(inserted)
  
  # Calculate proportion found that weren't inserted
  nWronglyFound <- length(isolatesOfHomoplasiesIdentified) - nFound
  if(nWronglyFound < 0){
    nWronglyFound <- 0
  }
  
  if(verbose){
    cat(paste("Proportion homoplasies found:", propFound, "\n"))
    cat(paste("Number homoplasies incorrectly identified:", nWronglyFound, "\n"))
  }
  
  return(c(propFound, nWronglyFound))
}

getIsolatesFromHomoplasyInfo <- function(homoplasyFinderOutput){
  
  # Ensure the isolate information are in character format
  homoplasyFinderOutput$`Node/IsolateID` <- 
    as.character(homoplasyFinderOutput$`Node/IsolateID`)
  homoplasyFinderOutput$IsolatesAlleleFoundIn <- 
    as.character(homoplasyFinderOutput$IsolatesAlleleFoundIn)
  
  # Initialise a list to store the isolates associated with each homoplasy
  output <- list()
  
  # Examine each homoplasy
  for(row in 1:nrow(homoplasyFinderOutput)){
    
    output[[row]] <- list(
      "foundIn"=strsplit(homoplasyFinderOutput[row, "Node/IsolateID"], split="-")[[1]],
      "from"=strsplit(homoplasyFinderOutput[row, "IsolatesAlleleFoundIn"], split="-")[[1]],
      "position"=homoplasyFinderOutput[row, "Position"]
    )
  }
  
  return(output)
}

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

buildPhylogenyOld <- function(sequences, file, output, rootOnREF, verbose){
  
  par(mar=c(0,0,0,0))
  
  # Build genetic distance matrix
  geneticDistances <- buildGeneticDistanceMatrix(sequences)
  
  # Build neighbour joining tree
  tree <- NJ(geneticDistances)

  # Root the tree against the reference sequence
  if(rootOnREF == TRUE){
    tree <- root(tree, outgroup="REF")
  }
  
  # Plot tree if wanting verbose information
  if(verbose){
    # Get the source-sink pairs from the homoplasy insertion output
    homoplasyIDs <- names(output)
    homoplasyIDs <- homoplasyIDs[homoplasyIDs != "sequences"]
    homoplasyPairs <- c()
    for(ID in homoplasyIDs){
      homoplasyPairs[length(homoplasyPairs) + 1] <- output[[ID]][["source"]]
      homoplasyPairs[length(homoplasyPairs) + 1] <- output[[ID]][["sink"]]
    }
    
    # Define the tip colors
    # Note will get warning when an isolate is involved in multiple homoplasy events
    tipColours <- rep("black", length(tree$tip.label))
    colours <- c("red", "blue", "green", "cyan", "orange", "darkorchid4",
                 "deeppink", "black", "brown", "darkolivegreen")
    for(i in 1:length(tree$tip.label)){
      
      if(tree$tip.label[i] %in% homoplasyPairs == FALSE){
        next
      }
      tipColours[i] <- colours[ceiling(which(homoplasyPairs == tree$tip.label[i])/2)]
    }
    
    # Plot the tree
    plot.phylo(tree, show.tip.label=TRUE, type="phylogram",
               edge.color="grey", edge.width=2,
               show.node.label=TRUE, label.offset=0.15,
               tip.color=tipColours, cex=1)
    
    # Get the axis limits
    axisLimits <- par("usr")
    
    # Add scale
    xPos <- axisLimits[2] - (0.1 * (axisLimits[2] - axisLimits[1]))
    yPos <- axisLimits[3] + (0.075 * (axisLimits[4] - axisLimits[3]))
    lines(x=c(xPos,xPos+1), y=c(yPos, yPos), lwd=4)
    text(x=xPos+0.5, y=axisLimits[3], labels="1 SNP", pos=3)
  }
  
  # Write the tree to file
  write.tree(tree, file=file, append=FALSE,
             digits=20, tree.names=FALSE)
  
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

insertHomoplasiesOld <- function(sequences, n, verbose){
  
  homoplasyEventInfo <- list()
  
  for(i in 1:n){
    
    # Insert a homoplasy
    output <- insertHomoplasy(sequences, verbose)
    
    # Get the sequences
    sequences <- output[["sequences"]]
    
    # Store the information about the homoplasy inserted
    homoplasyEventInfo[[as.character(i)]] <- list(
      "source" = output[["source"]],
      "sink" = output[["sink"]],
      "position" = output[["position"]],
      "allele" = output[["allele"]]
    )
  }
  
  homoplasyEventInfo[["sequences"]] <- sequences
  
  return(homoplasyEventInfo)
}

insertHomoplasyOld <- function(sequences, verbose){
  
  # Get the isolates IDs
  isolates <- names(sequences)
  
  # Remove reference sequence
  isolates <- isolates[isolates != "REF"]
  
  # Initialise an array to store the positions to choose from
  positionsToChooseFrom <- c()
  
  # Initialise variables to store the source and sink individuals
  i <- NULL
  j <- NULL
  
  # Search for suitable source and sink individuals
  while(length(positionsToChooseFrom) == 0){
    # Select a source and sink
    i <- sample(size=1, x=isolates)
    j <- sample(size=1, x=isolates[isolates != i])
    
    # Note the positions to choose from - not the same or reference
    positionsToChooseFrom <- 
      which(sequences[[i]] != sequences[[j]] &
              sequences[[i]] != sequences[["REF"]])
  }
  
  position <- positionsToChooseFrom[
    sample(size=1, x=1:length(positionsToChooseFrom))]
  
  sequences[[j]][position] <- sequences[[i]][position]
  
  if(verbose){
    cat(paste("Introduced homoplasy in isolate: ", j, "\tFrom: ", i, 
              "\tPosition: ", position, "\tAllele: ", 
              sequences[[i]][position], "\n", sep=""))
  }
  
  output <- list(
    "sequences" = sequences,
    "source" = i,
    "sink" = j,
    "position" = position,
    "allele" = sequences[[i]][position]
  )
  
  return(output)
}

buildSequences <- function(simulationOutput){
  
  # Get the mutation event information from the simulation output
  mutationEvents = simulationOutput[["events"]]
  mutationID = simulationOutput[["id"]]
  
  # Create the reference allele for each event
  nucleotides <- c("A", "C", "G", "T")
  eventsRefAlleles <- sample(nucleotides, size=mutationID, replace=TRUE)
  
  # Create the alternate allele for each event
  eventsAltAlleles <- c()
  for(i in 1:mutationID){
    eventsAltAlleles[i] <- sample(nucleotides[nucleotides != eventsRefAlleles[i]], size=1)
  }
  
  # Build each individuals sequence
  sequences <- list()
  for(key in names(mutationEvents)){
    
    # Skip unsampled individuals
    if(key %in% simulationOutput$sampled == FALSE){
      next
    }
    
    sequence <- eventsRefAlleles
    sequence[mutationEvents[[key]]] <- eventsAltAlleles[mutationEvents[[key]]]
    sequences[[key]] <- sequence[-1]
  }
  
  # Store the reference alleles
  sequences[["REF"]] <- eventsRefAlleles[-1]
  
  return(sequences)
}

runSimulation <- function(n, mutationRate, infectiousness, samplingProb,
                          nToSample, verbose){
  
  nSampled <- 0
  simulationOutput <- NULL
  
  while(nSampled < nToSample){
    
    simulationOutput <- simulateSequences3(n, mutationRate, infectiousness, samplingProb,
                                           nToSample, verbose)
    
    nSampled <- length(simulationOutput$events)
  }
  
  return(simulationOutput)
}

simulateSequences3 <- function(n, mutationRate, infectiousness, samplingProb,
                               nToSample, verbose){
  
  # Set population size
  susceptibles <- 1:n
  
  # Infect a random seed and move into infected class
  seed <- sample(size=1, susceptibles)
  infecteds <- c(seed)
  susceptibles <- susceptibles[susceptibles != seed]
  
  # Initialise an array to store the sampled individuals
  sampled <- c()
  
  # Initialise a list to store each sequence - during and after
  mutationEvents <- list()
  mutationID <- 1
  mutationEvents[[as.character(infecteds[1])]] <- c(mutationID)
  
  # Initialise a variable to record the timeStep
  timeStep <- 0
  
  # Initialise a list to record when each sequence was last checked
  timeStepSequencesLastChecked <- list()
  timeStepSequencesLastChecked[[as.character(infecteds[1])]] <- timeStep
  
  # Initialise a table to record the population state
  counts <- data.frame(Timestep=timeStep, 
                       Susceptibles=length(susceptibles), Infecteds=length(infecteds),
                       Sampled=length(sampled))
  
  # Run infection
  while(length(sampled) < nToSample && length(infecteds) != 0){
    
    # Increment the timestep
    timeStep <- timeStep + 1
    
    # Calculate the force of infection on a susceptible individual
    probInfected <- 1 - ((1 - infectiousness)^length(infecteds))
    
    # Get a distribution of random numbers and random sources to select from
    randomNumbers <- runif(length(susceptibles))
    
    # Note the infected susceptibles
    infectedSusceptibles <- susceptibles[which(randomNumbers < probInfected)]
    
    # Examine each susceptible individual that became infected
    if(length(infectedSusceptibles) > 0){
      
      # Note the sources for the infected individuals
      sources <- infecteds[sample(size=length(infectedSusceptibles),
                                  x=1:length(infecteds), replace=TRUE)]
      
      for(i in 1:length(infectedSusceptibles)){
        
        # Get the current susceptible individual
        susceptible <- infectedSusceptibles[i]
        
        # Get the random source for the current individual
        source <- sources[i]
        
        # Calculate the number of timesteps since the source sequence was mutated
        nTimeSteps <- timeStep - timeStepSequencesLastChecked[[as.character(source)]]
        
        # Mutate the sources sequence
        output <- generateMutationEvents(
          nTimeStepsSinceChecked=nTimeSteps,
          mutationRate=mutationRate,
          mutationEvents=mutationEvents[[as.character(source)]],
          mutationID=mutationID
        )
        mutationEvents[[as.character(source)]] <- output[["events"]]
        mutationID <- output[["id"]]
        timeStepSequencesLastChecked[[as.character(source)]] <- timeStep
        
        # Transmit sequence onto current susceptible individual
        mutationEvents[[as.character(susceptible)]] <- 
          mutationEvents[[as.character(source)]]
        timeStepSequencesLastChecked[[as.character(susceptible)]] <- timeStep
      }
      
      # Remove infected susceptible indivuduals
      susceptibles <- susceptibles[susceptibles %in% infectedSusceptibles == FALSE]
    }
    
    # Sample from infecteds
    randomNumbers <- runif(length(infecteds))
    sampledInfecteds <- infecteds[which(randomNumbers < samplingProb)]
    for(infected in sampledInfecteds){
      
      # Calculate the number of timesteps since the source sequence was mutated
      nTimeSteps <- timeStep - timeStepSequencesLastChecked[[as.character(infected)]]
      
      # Mutate the sources sequence
      output <- generateMutationEvents(
        nTimeStepsSinceChecked=nTimeSteps,
        mutationRate=mutationRate,
        mutationEvents=mutationEvents[[as.character(infected)]],
        mutationID=mutationID
      )
      mutationEvents[[as.character(infected)]] <- output[["events"]]
      mutationID <- output[["id"]]
      timeStepSequencesLastChecked[[as.character(infected)]] <- timeStep
      
      # Move the infected individual into the sampled array
      infecteds <- infecteds[infecteds != infected]
      sampled[length(sampled) + 1] <- infected
    }
    
    # Move all infected susceptible individuals into infecteds array
    infecteds <- c(infecteds, infectedSusceptibles)
    
    # Print current population status
    counts[nrow(counts) + 1, ] <- c(timeStep, length(susceptibles), length(infecteds),
                                    length(sampled))
  }
  
  # Plot the population dynamics
  if(verbose){
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    plot(x=counts$Timestep, y=counts$Susceptibles, type="o", col="blue",
         ylab="N. Individuals", xlab="Timestep", main="Simulation Dynamics", bty="n",
         ylim=c(0, n), xlim=c(0, max(counts$Timestep)), las=1, cex.lab=1.5,
         cex.main=1.5)
    points(x=counts$Timestep, y=counts$Infecteds, type="o", col="red")
    points(x=counts$Timestep, y=counts$Sampled, type="o", col="black")
    legend("right", legend=c("Susceptible", "Infected", "Sampled"),
           text.col=c("blue", "red", "black"), bty="n")
  }
  
  # Put the output information in a list
  output <- list(
    "events"=mutationEvents,
    "id"=mutationID,
    "sampled"=sampled
  )
  
  return(output)
}

generateMutationEvents <- function(nTimeStepsSinceChecked, 
                                   mutationRate, mutationEvents, mutationID){
  
  # Calculate how many mutations occurred
  nMutations <- sum(rpois(n=nTimeStepsSinceChecked, lambda=mutationRate))
  
  # Create each mutation event, give it and id and store it
  if(nMutations != 0){
    
    for(i in 1:nMutations){
      mutationID <- mutationID + 1
      mutationEvents[length(mutationEvents) + 1] <- mutationID
    }
  }
  
  # Store the outputs in a list
  output <- list(
    "events"=mutationEvents,
    "id"=mutationID)
  
  return(output)
}

simulateSequences2 <- function(n, genomeSize, mutationRate, infectiousness, samplingProb,
                               nToSample){
  
  # Note the nucleotides
  nucleotides <- c("A", "T", "G", "C")
  
  # Set population size
  susceptibles <- 1:n
  
  # Infect a random seed and move into infected class
  seed <- sample(size=1, susceptibles)
  infecteds <- c(seed)
  susceptibles <- susceptibles[susceptibles != seed]
  
  # Initialise an array to store the sampled individuals
  sampled <- c()
  
  # Initialise a list to store each sequence - during and after
  sequences <- c()
  sequences[["REF"]] <- sample(x=nucleotides, size=genomeSize, replace=TRUE)
  sequences[[as.character(infecteds[1])]] <- sequences[["REF"]]
  
  # Initialise a variable to record the timeStep
  timeStep <- 0
  
  # Initialise a list to record when each sequence was last checked
  timeStepSequencesLastChecked <- list("REF"=NA)
  timeStepSequencesLastChecked[[as.character(infecteds[1])]] <- timeStep
  
  # Initialise a table to record the population state
  counts <- data.frame(Timestep=timeStep, 
                       Susceptibles=length(susceptibles), Infecteds=length(infecteds),
                       Sampled=length(sampled))
  
  # Run infection
  cat("Starting simulation...")
  while(length(sampled) < nToSample && length(infecteds) != 0){
    
    # Increment the timestep
    timeStep <- timeStep + 1
    
    # Calculate the force of infection on a susceptible individual
    probInfected <- 1 - ((1 - infectiousness)^length(infecteds))
    
    # Get a distribution of random numbers and random sources to select from
    randomNumbers <- runif(length(susceptibles))
    
    # Note the infected susceptibles
    infectedSusceptibles <- susceptibles[which(randomNumbers < probInfected)]
    
    # Examine each susceptible individual that became infected
    if(length(infectedSusceptibles) > 0){
      
      # NOte the sources for the infected individuals
      sources <- infecteds[sample(size=length(infectedSusceptibles),
                                  x=1:length(infecteds), replace=TRUE)]
      
      for(i in 1:length(infectedSusceptibles)){
        
        # Get the current susceptible individual
        susceptible <- infectedSusceptibles[i]
        
        # Get the random source for the current individual
        source <- sources[i]
        
        # Calculate the number of timesteps since the source sequence was mutated
        nTimeSteps <- timeStep - timeStepSequencesLastChecked[[as.character(source)]]
        
        # Mutate the sources sequence
        sequences[[as.character(source)]] <- mutateSequence(
          nTimeStepsSinceChecked=nTimeSteps,
          mutationRate=mutationRate,
          sequence=sequences[[as.character(source)]],
          nucleotides=nucleotides
        )
        timeStepSequencesLastChecked[[as.character(source)]] <- timeStep
        
        # Transmit sequence onto current susceptible individual
        sequences[[as.character(susceptible)]] <- sequences[[as.character(source)]]
        timeStepSequencesLastChecked[[as.character(susceptible)]] <- timeStep
      }
      
      # Remove infected susceptible indivuduals
      susceptibles <- susceptibles[susceptibles %in% infectedSusceptibles == FALSE]
    }
    
    # Sample from infecteds
    randomNumbers <- runif(length(infecteds))
    sampledInfecteds <- infecteds[which(randomNumbers < samplingProb)]
    for(infected in sampledInfecteds){
      
      # Calculate the number of timesteps since the source sequence was mutated
      nTimeSteps <- timeStep - timeStepSequencesLastChecked[[as.character(infected)]]
      
      # Mutate the sources sequence
      sequences[[as.character(infected)]] <- mutateSequence(
        nTimeStepsSinceChecked=nTimeSteps,
        mutationRate=mutationRate,
        sequence=sequences[[as.character(infected)]],
        nucleotides=nucleotides
      )
      timeStepSequencesLastChecked[[as.character(infected)]] <- timeStep
      
      # Move the infected individual into the sampled array
      infecteds <- infecteds[infecteds != infected]
      sampled[length(sampled) + 1] <- infected
    }
    
    # Move all infected susceptible individuals into infecteds array
    infecteds <- c(infecteds, infectedSusceptibles)
    
    # Print current population status
    counts[nrow(counts) + 1, ] <- c(timeStep, length(susceptibles), length(infecteds),
                                    length(sampled))
  }
  cat("\rSimulation complete.\n")
  
  # Plot the population dynamics
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(x=counts$Timestep, y=counts$Susceptibles, type="o", col="blue",
       ylab="N. Individuals", xlab="Timestep", main="Simulation Dynamics", bty="n",
       ylim=c(0, n), xlim=c(0, max(counts$Timestep)), las=1)
  points(x=counts$Timestep, y=counts$Infecteds, type="o", col="red")
  points(x=counts$Timestep, y=counts$Sampled, type="o", col="black")
  legend("right", legend=c("Susceptible", "Infected", "Sampled"),
         text.col=c("blue", "red", "black"), bty="n")
  
  return(sequences)
}

mutateSequence <- function(nTimeStepsSinceChecked, mutationRate, sequence, nucleotides){
  
  if(nTimeStepsSinceChecked != 0){
    nMutations <- sum(rpois(n=nTimeStepsSinceChecked, lambda=mutationRate))
    
    positions <- sample(size=nMutations, x=1:length(sequence))
    
    for(position in positions){
      sequence[position] <- sample(size=1, x=nucleotides[nucleotides != sequence[position]])
    }
  }
  
  return(sequence)
}

mutateSequence <- function(nTimeStepsSinceChecked, mutationRate, sequence, nucleotides){
  
  if(nTimeStepsSinceChecked != 0){
    nMutations <- sum(rpois(n=nTimeStepsSinceChecked, lambda=mutationRate))
    
    positions <- sample(size=nMutations, x=1:length(sequence))
    
    for(position in positions){
      sequence[position] <- sample(size=1, x=nucleotides[nucleotides != sequence[position]])
    }
  }
  
  return(sequence)
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
