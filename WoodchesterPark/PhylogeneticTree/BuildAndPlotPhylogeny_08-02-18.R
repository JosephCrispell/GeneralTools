#### Packages ####
library(ape)
library(phangorn)
library(geiger) # For the tips function
library(plotrix) # For drawing circles
library(ips) # Using RAxML

#~~~~~~~~~~~~~~~~~~~~~#
#### Path and Date ####
#~~~~~~~~~~~~~~~~~~~~~#

path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/"

date <- format(Sys.Date(), "%d-%m-%y")

#~~~~~~~~~~~~~~~~~#
#### Run RaXML ####
#~~~~~~~~~~~~~~~~~#

# Note the fasta file
fastaFile <- paste(path, "vcfFiles/sequences_withoutHomoplasies_27-03-18.fasta", sep="")

# Run RAxML to produce a Maximum Likelihood phylogeny with bootstrap support values
# TAKES AGES!!! - WP data ~5 hours
treeBS <- runRAXML(fastaFile, date, nBootstraps=100, nThreads=2)

# Convert the branch lengths to SNPs
treeBS$edge.length <- treeBS$edge.length * getNSitesInFASTA(fastaFile)

# Parse the tip labels
treeBS$tip.label <- parseIsolateLabels(treeBS$tip.label)

# Make bootstrap values numeric
treeBS$node.label <- as.numeric(treeBS$node.label)
treeBS$node.label[is.na(treeBS$node.label)] <- 0

# Re-root tree
treeBS <- root(treeBS, outgroup="Ref-1997")

# Take an initial look at the phylogeny
viewRAxMLTree(treeBS, path)

#~~~~~~~~~~~~~~~~~~~~~~~~#
#### Clade definition ####
#~~~~~~~~~~~~~~~~~~~~~~~~#

# Open a pdf
file <- paste(path, "vcfFiles/", "mlTree_CladesAndLocations_", date, ".pdf", sep="")
pdf(file, height=10, width=10)

# Define branch colours by clade and plot the tree
nodesDefiningClades <- c(527, 496, 562, 421, 486)
cladeColours <- c("cyan", "magenta", "green", "darkorchid4", "brown")
plotTree(treeBS, plotBSValues=TRUE,
         nodes=nodesDefiningClades,
         colours=cladeColours)
plotTree(treeBS, plotBSValues=FALSE,
         nodes=nodesDefiningClades,
         colours=cladeColours)


# Plot the isolate locations and colour by clade
tipLabelsWithSamplingTimes <- 
  plotIsolateLocations(treeBS,
    badgerIsolateFile=paste(path, "IsolateData/", "BadgerInfo_08-04-15_LatLongs_XY_Centroids.csv",
                            sep=""),
    cattleIsolateFile=paste(path, "IsolateData/",
                            "CattleIsolateInfo_AddedNew_TB1484-TB1490_22-03-18.csv",
                            sep=""),
    colours=cladeColours, nodes=nodesDefiningClades,
    badgerCentre=c(381761.7, 200964.3), expand=75000, thresholdDistance=3500)

tipLabelsWithSamplingTimes <- 
  plotIsolateLocations(treeBS,
                       badgerIsolateFile=paste(path, "IsolateData/", "BadgerInfo_08-04-15_LatLongs_XY_Centroids.csv",
                                               sep=""),
                       cattleIsolateFile=paste(path, "IsolateData/",
                                               "CattleIsolateInfo_AddedNew_TB1484-TB1490_22-03-18.csv",
                                               sep=""),
                       colours=cladeColours, nodes=nodesDefiningClades,
                       badgerCentre=c(381761.7, 200964.3), expand=7000, thresholdDistance=3500)

# Close the pdf
dev.off()

# Write the clusters out to file
noteCladesOfIsolates(treeBS, nodesDefiningClades, 
                     file=paste(path, "vcfFiles/", "clusters_", date, ".csv", sep=""))

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Write tree to file ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Create a file
file <- paste(path, "vcfFiles/", "mlTree_", date, ".tree", sep="")

write.tree(treeBS, file = file, append = FALSE,
           digits = 20, tree.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### Prepare for TempEst ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Add dates to tips
treeBS$tip.label <- tipLabelsWithSamplingTimes

# Print full tree with dated tips
write.tree(treeBS, append = FALSE, digits = 20, tree.names = FALSE,
           file = paste(path, "vcfFiles/",
                        "mlTree_DatedTips_", date, ".tree", sep=""))

# Select the BASTA clade and print tree to file
selectBASTAClade(treeBS, node=433, plot=TRUE)

#-----------------#
#### FUNCTIONS ####
#-----------------#

viewRAxMLTree <- function(treeBS, filePath=path){

  # Set the margins
  par(mfrow=c(1,1))
  par(mar=c(0,0,0,0)) # Bottom, Left, Top, Right
    
  # Plot initial tree to find nodes defining clades
  pdf(paste(filePath, "test.pdf", sep=""), height=40, width=40)
    
  plot.phylo(treeBS, "fan")
  nodelabels()
    
  dev.off()
    
  # Reset margins
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

runRAXML <- function(fastaFile, date, nBootstraps, nThreads){
  
  # Create a directory for the output file
  directory <- paste(path, "vcfFiles/RAxML_", date, sep="")
  dir.create(directory)
  
  # Set the Working directory - this will be where the output files are dumped
  setwd(directory)
  
  # Note where RaxML is located
  raxmlExecutable <- shQuote("C:/Users/Joseph Crisp/Desktop/RAxML_v8.2.10/raxmlHPC-PTHREADS.exe", type="cmd")
  
  # Using system command - convert FASTA file path into windows format (deals with spaces in directory/file names)
  fastaFile <- shQuote(fastaFile, type="cmd")
  
  # Build analysis name
  analysisName <- paste("RaxML-R_", date, sep="")
  
  # Build the command
  model <- "GTRCAT" # No rate heterogenity
  seeds <- sample(1:100000000, size=2, replace=FALSE) # For parsimony tree and boostrapping
  
  command <- paste(raxmlExecutable, 
                   " -f a", # Algorithm: Rapid boostrap inference
                   " -N ", nBootstraps,
                   " -T ", nThreads,
                   " -m ", model, " -V", # -V means no rate heterogenity
                   " -p ", seeds[1], " -x ", seeds[2], # Parsimony and boostrapping seeds
                   " -n ", analysisName,
                   " -s ", fastaFile, sep="")
  system(command, intern=TRUE)
  
  # Get the tree and read it in
  treeBS <- getTreeFileWithSupportValues()
  
  return(treeBS)
}

getTreeFileWithSupportValues <- function(){
  
  # Get files in current working directory
  files <- list.files()
  
  # Select the tree file with BS support values
  treeBSFile <- files[grepl(files, pattern="RAxML_bipartitions[.]") == TRUE]
  
  # Open the file
  treeBS <- read.tree(treeBSFile)
}

defineTipShapesForSpecies <- function(tipLabels, cow, badger){
  
  shapes <- c()
  for(i in 1:length(tipLabels)){
    
    if(grepl(pattern="TB", x=tipLabels[i]) == TRUE){
      shapes[i] <- cow
    }else if(grepl(pattern="WB", x=tipLabels[i]) == TRUE){
      shapes[i] <- badger
    }else{
      shapes[i] <- cow
    }
  }
  
  return(shapes)
}

selectBASTAClade <- function(treeBS, node, plot=FALSE, filePath=path,
                             currentDate=date){

  if(plot == TRUE){
    pdf(file=paste(filePath, "vcfFiles/mlTree_BastaClade_",
                   currentDate, ".pdf", sep=""))
    
    branchColours <- defineBranchColoursOfClade(treeBS, node, "black", "lightgrey")
    
    plot.phylo(treeBS, "fan", edge.color=branchColours, edge.width=3,
               show.tip.label=FALSE)

    dev.off()
  }
  
  # Get basta clade
  bastaClade <- extract.clade(treeBS, node)
  
  # Print out tree
  write.tree(bastaClade, append=FALSE, digits=20, tree.names=FALSE,
             file=paste(filePath, "vcfFiles/",
                        "mlTree_BASTAClade_DatedTips_",
                        currentDate, ".tree", sep=""))
}

defineBranchColoursOfClade <- function(tree, nodeDefiningClade,
                                       colour, defaultColour){
  branchColours <- rep(defaultColour, dim(tree$edge)[1])
  clade <- tips(tree, node=nodeDefiningClade)
  branchesInClades <- which.edge(tree, clade)
  branchColours[branchesInClades] <- colour
  
  return(branchColours)
}

noteCladesOfIsolates <- function(treeBS, nodesDefiningClades, file){
  
  # Initialise two arrays to store the isolate IDs and clades
  isolates <- c()
  clades <- c()
  
  # Examine each clade
  for(i in 1:length(nodesDefiningClades)){
    tipsInClade <- tips(treeBS, nodesDefiningClades[i])
    
    for(tip in tipsInClade){
      isolates[length(isolates) + 1] <- tip
      clades[length(clades) + 1] <- i - 1
    }
  }
  
  # Combine the arrays into table
  isolateClades <- data.frame(ID=isolates, Cluster=clades, stringsAsFactors=FALSE)

  # Print out table
  write.table(isolateClades, file, quote=FALSE, sep=",", row.names=FALSE)
}

plotIsolateLocations <- function(treeBS,
                                 badgerIsolateFile, cattleIsolateFile,
                                 colours, badgerCentre, expand,
                                 thresholdDistance, nodes){
  
  # Read in the badger sampling information
  badgerInfo <- read.table(badgerIsolateFile, header=TRUE, stringsAsFactors=FALSE,
                         sep=",")
  badgerIsolateLocations <- noteBadgerIsolateSamplingLocations(badgerInfo)
  
  # Cattle Isolates
  cattleInfo <- read.table(cattleIsolateFile, header=TRUE, sep=",", 
                           stringsAsFactors=FALSE)
  cattleIsolateLocations <- noteCattleIsolateSamplingLocations(cattleInfo)
  
  # Create the clade colours - apply alpha
  cladeColoursRGB <- getRGBsOfColours(colours, alpha=0.75)

  # Note the isolates in each clade
  isolatesInClades <- findIsolatesInClades(treeBS, nodes)
  
  # Create an empty plot
  par(mar=c(0,0,0,0))
  plot(x=NULL, y=NULL, yaxt="n", xaxt="n", bty="n", ylab="",
       xlim=c(badgerCentre[1] - expand, badgerCentre[1] + expand), 
       ylim=c(badgerCentre[2] - expand, badgerCentre[2] + expand), asp=1,
       xlab="")
  
  # Plot a minimum convex polygon around the 
  # cattle and badger sampling locations for each cluster
  for(i in 1:length(colours)){
    
    # Get the isolates associated with the current clade
    isolates <- isolatesInClades[[as.character(i)]]
    
    # Get the coordinates of each isolate
    isolateCoordinates <- getXandYCoordinatesOfIsolates(isolates,
                                                        cattleIsolateLocations,
                                                        badgerIsolateLocations)
    
    # Remove NA rows - where couldn't find coordinates for isolates
    isolateCoordinates <- isolateCoordinates[is.na(isolateCoordinates$X) == FALSE, ]
    
    # Plot the points
    points(isolateCoordinates, 
           pch=ifelse(isolateCoordinates$Species == "BADGER", 19, 17),
           col=cladeColoursRGB[i], cex=2)
    
    # Add a convex hull around the points
    addPolygon(isolateCoordinates$X, isolateCoordinates$Y, colours[i])
  }
  
  # Add inner circle from BASTA deme assignment diagram
  draw.circle(x=badgerCentre[1], y=badgerCentre[2], radius=thresholdDistance,
              border="black", lty=2)
  text(x=badgerCentre[1], y=badgerCentre[2] - (thresholdDistance + 500),
       labels=paste(round(thresholdDistance/1000, digits=2), "km radius"))
  
  # Add legend
  legend("bottomleft", legend=c("Cow", "Badger"),
         pch=c(17, 16), col="black", pt.cex=2,
         text.col="black", bty='n')
  
  # Add the cluster numbers
  legend("bottomright", legend=addTextToArray("Cluster ", 0:(length(colours)-1), ""),
         text.col=colours, bty="n", cex=2)
  
  # Add sampling times for each isolate  to tip labels - needed for TempEst
  tipLabelsWithSamplingTimes <- addSamplingTimes(
    treeBS$tip.label, badgerInfo, cattleInfo)

  return(tipLabelsWithSamplingTimes)
}

addSamplingTimes <- function(tipLabels, badgerInfo, cattleInfo){
  
  # Initialise an array to store the new tip labels
  newTipLabels <- c()
  
  # Note the format of the dates
  dateFormat <- "%d/%m/%Y"
  
  # Examine each tip
  for(i in 1:length(tipLabels)){
    
    # Check if reference
    if(tipLabels[i] == "Ref-1997"){
      newTipLabels[i] <- paste(tipLabels[i], "_1997-01-01", sep="")
      
      # Check if badger or cow isolate
    }else if(grepl(x=tipLabels[i], pattern="WB") == TRUE){
      
      row <- which(badgerInfo$WB_id == tipLabels[i])
      
      newTipLabels[i] <- paste(tipLabels[i], "_", 
                               as.character(as.Date(badgerInfo[row, "date"], dateFormat)),
                               sep="")
      
    }else{
      
      row <- which(cattleInfo$StrainId == tipLabels[i])
      
      newTipLabels[i] <- paste(tipLabels[i], "_", 
                               as.character(as.Date(cattleInfo[row, "DateCultured"],
                                                    dateFormat)),
                               sep="")
    }
  }
  
  return(newTipLabels)
}

addTextToArray <- function(text, array, sep){
  
  output <- c()
  for(i in 1:length(array)){
    output <- paste(text, array, sep=sep)
  }
  return(output)
}

findIsolatesInClades <- function(tree, nodesDefiningClades){
  isolatesInClades <- list()
  for(i in 1:length(nodesDefiningClades)){
    isolatesInClades[[as.character(i)]] <- tips(tree, nodesDefiningClades[i])
  }
  
  return(isolatesInClades)
}

getRGBsOfColours <- function(colours, alpha){
  
  output <- c()
  for(i in 1:length(colours)){
    output[i] <- convertToRGB(colours[i], alpha)
  }
  
  return(output)
}

convertToRGB <- function(colour, alpha){
  
  rgbInfo <- col2rgb(colour)
  
  output <- rgb(rgbInfo["red", 1], rgbInfo["green", 1], rgbInfo["blue", 1], alpha=alpha*255,
                maxColorValue=255)
  
  return(output)
}

getXandYCoordinatesOfIsolates <- function(isolates, cattleIsolateLocations, badgerIsolateLocations){
  
  # Initialise a dataframe to store the X and Y coordinates of each isolate
  coords <- data.frame(X=rep(NA, length(isolates)), Y=rep(NA, length(isolates)), 
                       Species=rep(NA, length(isolates)), stringsAsFactors=FALSE)
  
  # Examine each isolate
  for(row in 1:length(isolates)){
    
    # Is the current isolate from a badger?
    if(grepl(x=isolates[row], pattern="WB") == TRUE){
      
      if(is.null(badgerIsolateLocations[[isolates[row]]]) == FALSE){
        coords[row, c(1,2)] <- badgerIsolateLocations[[isolates[row]]]
        coords[row, "Species"] <- "BADGER"
      }
      
    }else{
      if(is.null(cattleIsolateLocations[[isolates[row]]]) == FALSE){
        coords[row, c(1,2)] <- cattleIsolateLocations[[isolates[row]]]
        coords[row, "Species"] <- "COW"
      }
    }
  }
  
  return(coords)
}

addPolygon <- function(xValues, yValues, borderColour){
  hullPoints <- chull(xValues, yValues)
  hullPoints <- c(hullPoints, hullPoints[1])
  polygon(xValues[hullPoints], yValues[hullPoints], col = NA, border = borderColour)
}

noteCattleIsolateSamplingLocations <- function(cattleInfo){
  
  isolates <- list()
  
  for(row in 1:nrow(cattleInfo)){
    
    coords <- c()
    
    # Does centroid information exist for the current isolate?
    if(is.na(cattleInfo[row, "Mapx"]) == FALSE){
      
      coords[1] <- cattleInfo[row, "Mapx"]
      coords[2] <- cattleInfo[row, "Mapy"]
      
    }
    
    # Store sampling coordinates if found
    if(length(coords) > 0 && is.na(cattleInfo[row, "StrainId"]) == FALSE){
      isolates[[cattleInfo[row, "StrainId"]]] <- coords
    }
  }
  
  return(isolates)
}

noteBadgerIsolateSamplingLocations <- function(metadata){
  
  isolates <- list()
  
  for(row in 1:nrow(metadata)){
    
    coords <- c()
    
    # Does centroid information exist for the current isolate?
    if(is.na(metadata[row, "GroupCentroidX"]) == FALSE){
      
      coords[1] <- metadata[row, "GroupCentroidX"]
      coords[2] <- metadata[row, "GroupCentroidY"]
      
      # Does X and Y exist for sampled group?
    }else if(is.na(metadata[row, "SampledGrpX"]) == FALSE){
      
      coords[1] <- metadata[row, "SampledGrpX"]
      coords[2] <- metadata[row, "SampledGrpY"]
    }
    
    # Store sampling coordinates if found
    if(length(coords) > 0){
      isolates[[metadata[row, "WB_id"]]] <- coords
    }
  }
  
  return(isolates)
}

plotTree <- function(treeBS, nodes, colours, plotBSValues=FALSE){
  
  branchColours <- defineBranchColoursOfClades(treeBS, nodes,
                                               colours, "lightgrey")
  
  # Define the characteristics of the tips based on species
  tipShapes <- defineTipShapesForSpecies(treeBS$tip.label, 24, 21)
  tipColour <- defineTipColourBySpecies(treeBS, "blue", "red", "lightgrey",
                                        nodes)
  
  # Plot the phylogenetic tree
  plot.phylo(treeBS, show.tip.label=FALSE, "fan",
             edge.color=branchColours, edge.width=3)
  
  # Add bootstrap values
  if(plotBSValues == TRUE){
    nodelabels(pch=20, frame="none",
               col=rgb(0,0,0, treeBS$node.label/max(treeBS$node.label)))
    
    legend("bottomright", title="Boostrap values", legend=c(100, 75, 50), pch=20, 
           col=c(rgb(0,0,0, 1), rgb(0,0,0, 0.75), rgb(0,0,0, 0.5)), bty="n")
  }
  
  # Add node labels
  tiplabels(pch=tipShapes, bg=tipColour, col="dimgrey")
  
  # Add Legends
  legend("bottomleft", legend=c("Cow", "Badger"),
         pch=c(17, 16), cex=1, col=c("blue", "red"), 
         text.col=c("blue", "red"), bty='n')
  text(x=20, y=0, labels="AF2122/97")
  
  # Add Scale bar
  points(x=c(-20, 30), y=c(-130, -130), type="l", lwd=3)
  text(x=5, y=-135, labels="50 SNPs", cex=1, xpd=TRUE)
  
  # Add Clade labels
  text(x=-75, y=95, labels="3", col=colours[4], cex=2)
  text(x=92, y=-82, labels="4", col=colours[5], cex=2)
  text(x=-123, y=-30, labels="1", col=colours[2], cex=2)
  text(x=-92.5, y=-75, labels="0", col=colours[1], cex=2)
  text(x=15, y=-120, labels="2", col=colours[3], cex=2)

  # Reset margins
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

prepareAndViewTree <- function(treeBS, plot=FALSE, filePath=path, fastaFileName){
  
  # Parse the tip labels
  treeBS$tip.label <- parseIsolateLabels(treeBS$tip.label)
  
  # Get number of sites used in FASTA file - use to convert p-distances into n. SNPs
  nSites <- getNSitesInFASTA(fastaFileName)
  treeBS$edge.length <- treeBS$edge.length * nSites
  
  # Re-root tree
  treeBS <- root(treeBS, outgroup="Ref-1997")
  
  if(plot == TRUE){
    # Set the margins
    par(mfrow=c(1,1))
    par(mar=c(0,0,0,0)) # Bottom, Left, Top, Right
    
    # Plot initial tree to find nodes defining clades
    pdf(paste(filePath, "test.pdf", sep=""), height=40, width=40)
      
    plot.phylo(treeBS, "fan")
    nodelabels()
      
    dev.off()
    
    # Reset margins
    par(mar=c(5.1, 4.1, 4.1, 2.1))
  }
  
  return(treeBS)
}

testSubstitutionModels <- function(sequencesPhyDat){

  # Run model testing to select appropriate model
  modelTestResults <- modelTest(sequencesPhyDat, model = c("JC", "HKY", "GTR"))
  
  # Get best model
  cat(paste("Best substitution model:", 
            modelTestResults$Model[which.min(modelTestResults$AIC)], "\n"))
}

buildAndBootStrapMLTree <- function(nBootstraps, substitutionModel){

  #~~~~~~~~~~~~~~~~~~~~~~~~#
  #### Preliminary tree ####
  #~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # Build the distance matrix
  distanceMatrix <- dist.dna(sequencesDNAbin, model="JC69")
  cat("Built distance matrix")
  
  # Build neighbour joining tree
  initialTree <- nj(distanceMatrix)
  cat("\rBuilt initial tree")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #### Maximum Likelihood ####
  #~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # Compute likelihood of tree given sequences
  likelihoodObject <- pml(initialTree, sequencesPhyDat)
  
  # Set maximum likelihood controls
  controls <- pml.control(maxit=100000, trace=0)
  
  # Run maximum likelihood
  fittingOutput <- optim.pml(likelihoodObject, 
                             optNni = TRUE,       # Optimise topology
                             optInv = TRUE,       # Optimise proportion of variable sites
                             model = substitutionModel,       # Substitution model
                             rearrangement="NNI", # Nearest Neighbour Interchanges
                             control=controls)
  cat("\rRan Maximum Likelihood tree estimation")
  
  #~~~~~~~~~~~~~~~~~~~~~#
  #### Bootstrapping ####
  #~~~~~~~~~~~~~~~~~~~~~#
  
  # Bootstrap the result of maximum likelihood
  bootstrapResults <- bootstrap.pml(fittingOutput, bs = nBootstraps, optNni = TRUE,
                                    jumble=TRUE)
  cat("\rRan Boostrapping of Maximum Likelihood tree estimation")
  
  # Get phylogenetic tree with bootstrap values
  treeBS <- plotBS(fittingOutput$tree, bootstrapResults, p = 50, type="phylogram")
  cat("\rBuilt tree. Finished..\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n")
  
  return(treeBS)
}

getTipsInClades <- function(tree, nodesDefiningClades){
  tipsInClades <- list()
  for(node in nodesDefiningClades){
    tipsInClades[[as.character(node)]] <- tips(tree, node)
  }
  return(tipsInClades)
}

defineTipColourBySpecies <- function(tree, cow, badger, defaultColour, nodesDefiningClades){
  
  tipColours <- rep(defaultColour, length(tree$tip.label))
  tipsInClades <- getTipsInClades(tree, nodesDefiningClades)
  for(tipIndex in 1:length(tree$tip.label)){
    
    for(nodeIndex in 1:length(nodesDefiningClades)){
      
      if(tree$tip.label[tipIndex] %in% tipsInClades[[as.character(nodesDefiningClades[nodeIndex])]] == TRUE){
        if(grepl(pattern="TB|HI-|AF-", x=tree$tip.label[tipIndex]) == TRUE){
          tipColours[tipIndex] <- cow
        }else if(grepl(pattern="WB", x=tree$tip.label[tipIndex]) == TRUE){
          tipColours[tipIndex] <- badger
        }else{
          tipColours[tipIndex] <- defaultColour
        }
        break
      }
    }
  }
  
  return(tipColours)
}

defineTipSizeBySequencingQuality <- function(tipLabels, isolateQuality){
  
  tipQuality <- c()
  for(i in 1:length(tipLabels)){
    if(tipLabels[i] != "Ref-1997"){
      
      tipQuality[i] <- isolateQuality[[tipLabels[i]]]
    }else{
      tipQuality[i] <- 1
    }
  }
  
  return(tipQuality)
}

getIsolateQuality <- function(table){
  isolateQuality <- list()
  for(i in 1:nrow(table)){
    isolateQuality[[table[i, "Isolate"]]] <- table[i, "Coverage"]
  }
  
  return(isolateQuality)
}

defineBranchColoursOfClades <- function(tree, nodesDefiningClades,
                                        CladeColours, defaultColour){
  branchColours <- rep(defaultColour, dim(tree$edge)[1])
  for(i in 1:length(nodesDefiningClades)){
    clade <- tips(tree, node=nodesDefiningClades[i])
    branchesInClades <- which.edge(tree, clade)
    branchColours[branchesInClades] <- cladeColours[i]
  }
  return(branchColours)
}

parseIsolateLabels <- function(isolateNames){
  
  output <- c()
  for(i in 1:length(isolateNames)){
    
    if(isolateNames[i] != "Ref-1997"){
      output[i] <- strsplit(isolateNames[i], split="_")[[1]][1]
    }else{
      output[i] <- isolateNames[i]
    }
    
    if(grepl(isolateNames[i], pattern=">") == TRUE){
      output[i] <- substr(output[i], start=2, stop=nchar(output[i]))
    }
  }
  
  return(output)
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

