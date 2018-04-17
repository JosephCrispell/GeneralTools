#-------------------#
#### Preparation ####
#-------------------#

# Packages
library(phangorn) # Maximum likelihood phylogeny
library(gplots)
library(geiger) # For the tips function
library(grid) # Used to plot lines between plot panels
library(ape) # ladderise() - orders nodes in phylogeny

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Homoplasmy/"

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

#-----------------------#
#### Run simulations ####
#-----------------------#

# Simulation settings
popSize <- 300
mutationRate <- 0.5
infectiousness <- 0.001
samplingProb <- 0.05
nToSample <- 150
nHomoplasiesValues <- 0:10
nSimulations <- 1000

# Initialise a table to store the results of HomoplasyFinder on simulated data
results <- data.frame("NFoundTrue"=rep(NA, nSimulations*length(nHomoplasiesValues)), 
                      "NIncorrectTrue"=rep(NA, nSimulations*length(nHomoplasiesValues)),
                      "NFoundAfter"=rep(NA, nSimulations*length(nHomoplasiesValues)), 
                      "NIncorrectAfter"=rep(NA, nSimulations*length(nHomoplasiesValues)),
                      "Seed"=rep(NA, nSimulations*length(nHomoplasiesValues)), 
                      "NHomoplasies"=rep(NA, nSimulations*length(nHomoplasiesValues)))

# Set the working directory
setwd(path)

# Run the simulations
seeds <- sample(1:100000000, size=nSimulations*length(nHomoplasiesValues), replace=FALSE)
row <- 0
for(nHomoplasies in nHomoplasiesValues){
  for(i in 1:nSimulations){
    row <- row + 1
    
    # Progress
    cat(paste("\rRunning simulation", i, "of", nSimulations, "inserting", nHomoplasies, "homoplasy(ies)     "))

    # Get the current seed
    results[row, "Seed"] <- seeds[row]
    set.seed(results[row, "Seed"])

    # Note the number of homoplasies being inserted
    results[row, "NHomoplasies"] <- nHomoplasies
    
    # Generate the sequences
    simulationOutput <- runSimulation(popSize, mutationRate, infectiousness,
                                      samplingProb, nToSample)

    # Build the sequences based upon the mutation events
    sequences <- buildSequences(simulationOutput)
    sequences[["REF"]] <- NULL
    
    # Build phylogeny
    trueTreeFile <- paste("example-TRUE_", date, ".tree", sep="")
    treeBefore <- buildPhylogeny(sequences, trueTreeFile, maximumLikelihood=TRUE)
    
    # Insert homoplasies
    homoplasyInsertionInfo <- insertHomoplasies(sequences, n=nHomoplasies, tree=treeBefore)
    sequences <- homoplasyInsertionInfo[["sequences"]]
    
    # Build FASTA
    fastaFile <- paste("example_", date, ".fasta", sep="")
    writeFasta(sequences, fastaFile)
    
    # Build phylogeny
    treeFile <- paste("example-AFTER_", date, ".tree", sep="")
    buildPhylogeny(sequences, treeFile, maximumLikelihood=TRUE)
    
    #######
    ## Run HomoplasyFinder using the TRUE tree
    system(paste("java -jar HomoplasyFinder_06-03-18.jar", 0, fastaFile, trueTreeFile, sep=" "),
           ignore.stdout=FALSE)
    
    # Get and Check HomoplasyFinder output
    file <- paste(path, "homoplasyReport_", date, ".txt", sep="")
    homoplasyFinderOutput <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                                        check.names=FALSE)
    results[row, c("NFoundTrue","NIncorrectTrue")] <- 
      reportHowManyHomoplasiesWereFound(homoplasyFinderOutput, homoplasyInsertionInfo)
    
    #######
    ## Run HomoplasyFinder using the TRUE tree
    system(paste("java -jar HomoplasyFinder_06-03-18.jar", 0, fastaFile, treeFile, sep=" "),
           ignore.stdout=FALSE)
    
    # Get and Check HomoplasyFinder output
    file <- paste(path, "homoplasyReport_", date, ".txt", sep="")
    homoplasyFinderOutput <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                                        check.names=FALSE)
    results[row, c("NFoundAfter","NIncorrectAfter")] <- 
      reportHowManyHomoplasiesWereFound(homoplasyFinderOutput, homoplasyInsertionInfo)
  }
}

#### Plot results ####

file <- paste(path, "TestingHomoplasyFinder_", popSize, "-", mutationRate,
              "-", infectiousness, "-", samplingProb, "-", nToSample, "_",
              min(nHomoplasiesValues), "-", max(nHomoplasiesValues), "_", nSimulations, "_", date, ".pdf", sep="")
pdf(file)

plot(x=NULL, y=NULL, xlim=range(results$NHomoplasies), ylim=c(0,1), las=1, bty="n",  
     xlab="Number homoplasies inserted", ylab="Proportion present",
     main="Identifying inserted homoplasies using HomoplasyFinder")
for(i in unique(results$NHomoplasies)){
  
  subset <- results[results$NHomoplasies == i, ]
  counts <- table(subset$NFoundAfter / subset$NHomoplasies)
  if(nrow(counts) == 0){
    counts <- table(subset$NFoundAfter)
  }
  
  points(x=rep(i, length(counts)), y=as.numeric(names(counts)), pch=19, cex=6, col=rgb(0,0,0, counts/nrow(subset)), xpd=TRUE)
  text(x=rep(i, length(counts)), y=as.numeric(names(counts)),
       labels=round(counts/nrow(subset), digits=2),
       col=rgb(1,0,0, 1), cex=1)
}

plot(x=NULL, y=NULL, 
     xlim=range(results$NHomoplasies),
     ylim=range(results$NIncorrectAfter), las=1, bty="n",  
     xlab="Number homoplasies inserted", 
     ylab="Number non-inserted homoplasies Present", 
     main="Identifying non-inserted homoplasies using HomoplasyFinder")
for(i in unique(results$NHomoplasies)){
  
  subset <- results[results$NHomoplasies == i, ]
  counts <- table(subset$NIncorrectAfter)
  if(nrow(counts) == 0){
    counts <- table(subset$NIncorrectAfter)
  }
  
  points(x=rep(i, length(counts)), y=as.numeric(names(counts)), pch=19, cex=6, col=rgb(0,0,0, counts/nrow(subset)), xpd=TRUE)
  text(x=rep(i, length(counts)), y=as.numeric(names(counts)),
       labels=round(counts/nrow(subset), digits=2),
       col=rgb(1,0,0, 1), cex=1)
}

plot(x=results$NFoundTrue, y=results$NFoundAfter,
     xlim=range(results$NHomoplasies), ylim=range(results$NHomoplasies), las=1, bty="n",  
     xlab="Number inserted homoplasies present on tree built without homoplasies", 
     ylab="Number inserted homoplasies present on tree built with homoplasies",
     pch=19, cex=6, col=rgb(0,0,0, 0.005), xpd=TRUE,
     main="Identifying inserted homoplasies on tree\n built with and without them")
for(i in 0:max(results$NHomoplasies)){
  
  for(j in 0:max(results$NHomoplasies)){
    
    count <- nrow(results[results$NFoundTrue == i & results$NFoundAfter == j, ])
    if(count != 0){
      text(x=i, y=j,
           labels=round(count/nrow(results), digits=2),
           col=rgb(0,0,1, 1), cex=1)
    }
  }
}

plot(x=results$NIncorrectTrue, y=results$NIncorrectAfter,
     xlim=range(results$NHomoplasies), ylim=range(results$NHomoplasies), las=1, bty="n",  
     xlab="Number non-inserted homoplasies present on tree built without homoplasies", 
     ylab="Number non-inserted homoplasies present on tree built with homoplasies",
     pch=19, cex=6, col=rgb(0,0,0, 0.005), xpd=TRUE,
     main="Identifying non-inserted homoplasies on tree\n built with and without them")
for(i in 0:max(results$NHomoplasies)){
  
  for(j in 0:max(results$NHomoplasies)){
    
    count <- nrow(results[results$NIncorrectTrue == i & results$NIncorrectAfter == j, ])
    if(count != 0){
      text(x=i, y=j,
           labels=round(count/nrow(results), digits=2),
           col=rgb(0,0,1, 1), cex=1)
    }
  }
}

results$Proportion <- results$NFoundAfter / results$NFoundTrue
results$Proportion[is.nan(results$Proportion)] <- 0
plot(x=results$NFoundTrue, y=results$Proportion,
     pch=19, cex=6, col=rgb(0,0,0, 0.01), bty="n", xpd=TRUE, las=1,
     xlab="Number homoplasies present on tree built without homoplasies",
     ylab="Proportion present on tree built with homoplasies",
     main="Number homoplasies present on true tree that were\n identified on tree built using homoplasies")

for(i in unique(results$NFoundTrue)){
  
  for(j in unique(results$Proportion)){
    
    count <- nrow(results[results$NFoundTrue == i & results$Proportion == j, ])
    if(count != 0){
      
      text(x=i, y=j,
           labels=round(count/nrow(results), digits=2),
           col=rgb(0,0,1, 1), cex=1)
    }
  }
}

dev.off()

# Write the results to file
file <- paste(path, "TestingHomoplasyFinder_", popSize, "-", mutationRate,
              "-", infectiousness, "-", samplingProb, "-", nToSample, "_",
              min(nHomoplasiesValues), "-", max(nHomoplasiesValues), "_", nSimulations, "_", date, ".csv", sep="")
write.table(results, file, row.names=FALSE, quote=FALSE, sep=",")

#### Generation example ####

# Simulation settings
nHomoplasies <- 2

for(i in 1:100){
  file <- paste(path, "Examples/ManuscriptExample_", popSize, "-", mutationRate,
                "-", infectiousness, "-", samplingProb, "-", nToSample, "_", nHomoplasies, "_", i, "_", date, ".pdf", sep="")
  pdf(file)
  
  # Generate the sequences
  par(mfrow=c(1,1))
  simulationOutput <- runSimulation(popSize, mutationRate, infectiousness,
                                    samplingProb, nToSample, TRUE)
  
  # Build the sequences based upon the mutation events
  sequences <- buildSequences(simulationOutput)
  sequences[["REF"]] <- NULL
  
  par(mfrow=c(1,2))
  par(mar=c(0,0,0,0))
  
  # Build phylogeny
  trueTreeFile <- paste(path, "Examples/ManuscriptExample_TRUE_", date, ".tree", sep="")
  treeBefore <- buildPhylogeny(sequences, trueTreeFile, maximumLikelihood=TRUE)
  
  # Insert homoplasies
  homoplasyInsertionInfo <- insertHomoplasies(sequences, tree=treeBefore, n=nHomoplasies, verbose=TRUE)
  sequences <- homoplasyInsertionInfo[["sequences"]]
  plotPhylogeny(treeBefore, homoplasyInsertionInfo=homoplasyInsertionInfo, verbose=TRUE, tipCex=0.6, margins=FALSE,
                showScale=TRUE, showTips=FALSE)
  
  # Get tip coordinates
  tipsPriorToHomoplasies <- getTipCoordinates(simulationOutput$sampled)
  
  # Build FASTA
  writeFasta(sequences, paste(path, "Examples/ManuscriptExample_", date, ".fasta", sep=""))
  
  # Build phylogeny
  treeAfter <- buildPhylogeny(sequences, paste(path, "Examples/ManuscriptExample_AFTER_", date, ".tree", sep=""), 
                              maximumLikelihood=TRUE)
  plotPhylogeny(treeAfter, homoplasyInsertionInfo=NULL, verbose=TRUE, tipCex=0.6, margins=FALSE, showScale=FALSE,
                direction="leftwards", showTips=FALSE)
  
  # Get tip coordinates
  tipsPostHomoplasies <- getTipCoordinates(simulationOutput$sampled)
  
  # Plot lines between tips that have changed location
  plotLinesBetweenTips(tipsPriorToHomoplasies, tipsPostHomoplasies, col=rgb(1,0,0, 0.5), lwd=2)
  
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  dev.off()
}

#############
# FUNCTIONS #
#############

plotLinesBetweenTips <- function(tipLocationsA, tipLocationsB, col="black", lwd=1, lty=1){
  
  pushViewport(viewport())
  popViewport()
  
  for(key in names(tipLocationsA)){
    
    pushViewport(viewport())
    if(tipLocationsA[[key]][2] != tipLocationsB[[key]][2]){
      grid.lines(x = c(tipLocationsA[[key]][1], tipLocationsB[[key]][1]), 
                 y = c(tipLocationsA[[key]][2], tipLocationsB[[key]][2]), 
                 gp = gpar(col=col, lty=lty, lwd=lwd))
    }
    popViewport()
  }
}

getTipCoordinates <- function(tipLabels){
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tips <- list()
  for(i in 1:length(tipLabels)){
    tips[[as.character(tipLabels[i])]] <- c(grconvertX(lastPP$xx[i], "user", "ndc"), 
                                            grconvertY(lastPP$yy[i], "user", "ndc"))
  }
  
  return(tips)
}

reportHowManyHomoplasiesWereFound <- function(homoplasyFinderOutput, homoplasyInsertionInfo, nHomoplasiesValues){
  
  # Initialise a variable to record how many of the inserted homoplasies were found
  nFound <- 0
  
  # Examine each of the homoplasies inserted
  homoplasies <- 1:(length(homoplasyInsertionInfo) - 2)
  if(length(homoplasyInsertionInfo) - 2 != 0){
    
    for(key in as.character(homoplasies)){
      
      # Search for current homoplasy in HomoplasyFinder output
      rows <- which(homoplasyFinderOutput$Position == homoplasyInsertionInfo[[key]]$position)
      if(length(rows) == 0){
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
  }
  
  # Calculate proportion found and number of false positives (sequences also present in insertion info)
  #propFound <- nFound / (length(homoplasyInsertionInfo) - 1)
  
  # Calculate proportion found that weren't inserted
  nWronglyFound <- nrow(homoplasyFinderOutput) - nFound
  if(nWronglyFound < 0){
    nWronglyFound <- 0 # If a homoplasy is inserted at the same site then this number could be negative
  }
  
  return(c(nFound, nWronglyFound))
}

findNode <- function(nodes, isolates, verbose=FALSE){
  
  # Initialise a the id of the node to return
  id <- NULL
  
  # Examine each node
  for(key in names(nodes)){
    
    # Check if the current node represents the set of isolates
    if(length(nodes[[key]]) == length(isolates) && length(which(isolates %in% nodes[[key]])) == length(isolates)){
      id <- key
      break
    }
  }
  
  if(verbose == TRUE && is.null(id) == TRUE){
    cat(paste("Unable to find node associated with isolates:", isolates, "\n"))
  }
  
  return(id)
}

plotPhylogeny <- function(tree, homoplasyInsertionInfo=NULL, treeType="phylogram", tipCex=1, showTips=TRUE,
                          showScale=TRUE, direction="rightwards", labelOffset=0.15, verbose=FALSE, margins=FALSE){
  # Set the margins
  if(margins == TRUE){
    par(mar=c(0,0,0,0))
  }
  
  # Order the nodes in the phylogeny by clade
  tree <- ladderize(tree)
  
  # Label nodes associated with homoplasies - if homoplasy information available
  nodeColours <- rep(rgb(0,0,0, 0), length(tree$edge.length) + 1)
  if(is.null(homoplasyInsertionInfo) == FALSE){
    
    # Get the isolates associated with each node in the phylogeny
    nodes <- getNodes(tree)
    
    # Define the colour of the nodes based on the inserted homoplasies
    colours <- c("blue", "green", "cyan", "orange", "darkorchid4",
                 "deeppink", "black", "brown", "darkolivegreen")
    for(homoplasyId in names(homoplasyInsertionInfo)[1:(length(homoplasyInsertionInfo)-2)]){
      
      # Get the colour for the current homoplasy
      colour <- colours[as.numeric(homoplasyId) %% length(colours)]
      if(length(colour) == 0){
        colour <- colours[length(colours)]
      }
      
      # Identify the source and sink nodes on the current phylogeny (note node numbers may have changed)
      sourceNode <- findNode(nodes, homoplasyInsertionInfo[[homoplasyId]]$sourceIsolates, verbose)
      sinkNode <- findNode(nodes, homoplasyInsertionInfo[[homoplasyId]]$sinkIsolates, verbose)
      
      # Assign that colour to the source and sink nodes for the current homoplasy
      nodeColours[as.numeric(sourceNode)] <- colour
      nodeColours[as.numeric(sinkNode)] <- colour
    }
  }
  
  # Plot the tree
  plot.phylo(tree, show.tip.label=showTips, type=treeType,
             edge.color="grey", edge.width=2,
             show.node.label=FALSE, label.offset=labelOffset,
             tip.color="black", cex=tipCex, direction=direction)
  
  # Add node circles to highlight homoplasies
  nodelabels(node=1:(length(tree$edge.length) + 1), pch=19, col=nodeColours, cex=1)
  
  if(showScale == TRUE){
    # Get the axis limits
    axisLimits <- par("usr")
    
    # Add scale
    xPos <- axisLimits[2] - (0.1 * (axisLimits[2] - axisLimits[1]))
    yPos <- axisLimits[3] + (0.075 * (axisLimits[4] - axisLimits[3]))
    lines(x=c(xPos,xPos+1), y=c(yPos, yPos), lwd=4, xpd=TRUE)
    text(x=xPos+0.5, y=axisLimits[3], labels="1 SNP", pos=3, xpd=TRUE)  
  }
  
  # Reset the margins
  if(margins == TRUE){
    par(mar=c(5.1, 4.1, 4.1, 2.1))
  }
}

buildPhylogeny <- function(sequences, file, maximumLikelihood=FALSE){
  
  # Build genetic distance matrix
  geneticDistances <- dist.dna(as.DNAbin.alignment(as.alignment(sequences)), model="JC69")
  
  # Build neighbour joining tree
  njTree <- NJ(geneticDistances)
  
  # Build maximum likelihood tree if requested
  if(maximumLikelihood == TRUE){
    
    # Remove negative branch lengths from NJ tree
    njTree$edge.length[njTree$edge.length < 0] <- 0
    
    # Convert the sequences into a phyDat format
    sequencesPhyDat <- phyDat(sequences, type="DNA", levels=NULL)
    
    # Compute likelihood of tree given sequences
    likelihoodObject <- pml(njTree, sequencesPhyDat)
    
    # Set maximum likelihood controls
    controls <- pml.control(maxit=100000, trace=0)
    
    # Run maximum likelihood
    fittingOutput <- optim.pml(likelihoodObject, 
                               optNni = TRUE,       # Optimise topology
                               optInv = FALSE,       # Optimise proportion of variable sites
                               model = "JC",       # Substitution model
                               rearrangement="NNI", # Nearest Neighbour Interchanges
                               control=controls)
    
    # Get the tree
    tree <- fittingOutput$tree
    
    # Convert edge lengths to SNPs
    tree$edge.length <- tree$edge.length * length(sequences[[1]])
    
  }else{
    tree <- njTree
    tree$edge.length <- tree$edge.length * length(sequences[[1]])
  }
  
  # Write the tree to file
  write.tree(tree, file=file, append=FALSE,
             digits=20, tree.names=FALSE)
  
  return(tree)
}

insertHomoplasies <- function(n, sequences, tree, verbose=FALSE){
  
  # Get the isolates associated with each node in the phylogeny
  nodes <- getNodes(tree)
  
  # Initialise a list to store the homoplasy information
  output <- list()
  
  # Initialise arrays to the store the nodes and sites already used for homoplasies
  usedNodePairs <- c()
  usedSites <- c()
  
  # Insert the homoplasies
  if(n != 0){
    for(i in 1:n){
      
      # Reset the variables associated with the previous homoplasy
      source <- NULL
      sink <- NULL
      sitesUniqueToIsolatesThatArentUsed <- c()
      
      # Randomly pick source and sink nodes and identify potential positions to create a homoplasy with
      while(is.null(source) == TRUE || is.null(sink) == TRUE || length(sitesUniqueToIsolatesThatArentUsed) == 0){
        
        # Randomly select a node to act as source
        source <- sample(names(nodes), size=1)
        
        # Randomly select a node to act as a sink - as long as it isn't on the path to the root
        sink <- randomlySelectSinkNode(nodes, source)
        
        # Check that haven't already used this node pair
        if(paste(source, sink, sep=":") %in% usedNodePairs == TRUE){
          next
        }else{
          usedNodePairs[length(usedNodePairs) + 1] <- paste(source, sink, sep=":")
        }
        
        # Find the sites that are common and unique to all the isolates under the current node
        sitesUniqueToIsolates <- findConservedSitesThatAreUniqueToClade(nodes[[source]], sequences)
        sitesUniqueToIsolatesThatArentUsed <- sitesUniqueToIsolates[which(sitesUniqueToIsolates %in% usedSites == FALSE)]
      }
      
      # Randomly choose a site
      site <- sitesUniqueToIsolatesThatArentUsed[sample(1:length(sitesUniqueToIsolatesThatArentUsed), size=1)]
      allele <- sequences[[nodes[[source]][1]]][site]
      usedSites[length(usedSites) + 1] <- site
      
      # Assign the randomly chosen site and allele to those isolates
      sequences <- assignSiteAndAllele(nodes[[sink]], sequences, site, allele)
      
      # Record homoplasy
      if(verbose){
        cat(paste("---------------------------------------------------------------------------\n",
                  "Introduced homoplasy: From node: ", source, " to ", sink, 
                  "\tPosition: ", site, "\tAllele: ", allele,
                  "\nSource Isolate(s): ", paste(nodes[[source]], collapse=", "),
                  "\nSink Isolate(s): ", paste(nodes[[sink]], collapse=", "), "\n", sep=""))
      }
      output[[as.character(i)]] <- list(
        "sourceNode" = source,
        "sinkNode" = sink,
        "position" = site,
        "allele" = allele,
        "sourceIsolates"=nodes[[source]],
        "sinkIsolates"=nodes[[sink]])
    }
  }
  
  # Add the sequences to the output
  output[["sequences"]] <- sequences
  output[["tree"]] <- tree
  
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
    
    # Assign the site and allele to the current isolate
    sequences[[isolate]][site] <- allele
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
      
      # Skip isolates associated with current clade
      if(id %in% isolates){
        next
      }
      
      # Check if conserved allele in clade found in an isolate outside of clade - skip if found
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
  
  # Nodes number 1:nTips and then onwards: two methods to calculate total number:
  # - Number of edges + 1
  # - Number of tips plus number internal nodes: length(tree$tip.label) + tree$Nnode
  for(node in 1:(length(tree$tip.label) + tree$Nnode)){
    nodes[[as.character(node)]] <- tips(tree, node)
  }
  
  return(nodes)
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
                          nToSample, verbose=FALSE){
  
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
    
    # Check some individuals became infected
    if(length(infectedSusceptibles) > 0){
      
      # Note the sources for the infected individuals
      sources <- infecteds[sample(size=length(infectedSusceptibles),
                                  x=1:length(infecteds), replace=TRUE)]
      
      # Examine each susceptible individual that became infected
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
    legend("left", legend=c("Susceptible", "Infected", "Sampled"),
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
  
  # Create each mutation event, give it an id and store it
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
