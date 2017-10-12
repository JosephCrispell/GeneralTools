# Load the ape package
library(ape)
library(geiger) # For the tips function
library(plotrix)

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_13-07-17/"

###################################
# Get the Maximum Likelihood Tree #
###################################

# Read in the newick tree
file <- paste(path, "vcfFiles/",
              "mlTree_29-09-2017.tree", sep="")
tree <- read.tree(file=file)

# Drop badger beside reference - WB129 - THIS WILL MESS UP NODES DEFINING CLADES
# tree <- drop.tip(tree, "WB129")

# Convert Branch lengths to SNPs
fastaLength <- 8893
tree$edge.length <- tree$edge.length * fastaLength

##################################
# Get the Isolate Coverage Table #
##################################

file <- paste(path,"vcfFiles/",
              "IsolateVariantPositionCoverage_RESCUED_29-09-2017.txt", sep="")
table <- read.table(file, header=TRUE, stringsAsFactors=FALSE)

table$Isolate <- getIsolateIDFromFileNames(table$Isolate)

##############################
# Plot the Phylogenetic Tree #
##############################

file <- paste(path, "vcfFiles/", "mlTree_CladesAndLocations_02-10-17.pdf", sep="")
pdf(file, height=10, width=10)

# Set the margins
par(mfrow=c(1,1))
par(mar=c(0,0,0,0)) # Bottom, Left, Top, Right

plotType <- "fan" # "phylogram", "cladogram", "fan", "unrooted", "radial"

# Plot initial tree to find nodes defining clades
#pdf(paste(path, "vcfFiles/", "test.pdf", sep=""), height=40, width=40)
# 
#plot.phylo(tree, "fan")
#nodelabels()
# 
#dev.off()

# Define branch colours by clade
nodesDefiningClades <- c(521, 305, 332, 382) # use nodelabels() to show node numbers
cladeColours <- c("cyan", "pink", "green", "darkorchid4")
branchColours <- defineBranchColoursOfClades(tree, nodesDefiningClades,
                                             cladeColours, "lightgrey")

# Get each isolate's quality
isolateQuality <- getIsolateQuality(table)

# Plot the phylogenetic tree
plot.phylo(tree, show.tip.label=FALSE, plotType,
           edge.color=branchColours, edge.width=3,
           show.node.label=TRUE)

# Add node labels
nodelabels(node=1:length(tree$tip.label), 
           cex=defineTipSizeBySequencingQuality(tree$tip.label, isolateQuality),
           pch=defineTipShapesForSpecies(tree$tip.label, 24, 21),
           bg=defineTipColourBySpecies(tree, "blue", "red", "lightgrey", nodesDefiningClades),
           col="dimgrey")

# Add Legends
#text(x=140, y=-130, labels="Variant Position Coverage:", col="black", cex=1)
#addLegendForQuality("bottomright", 1)
legend("bottomleft", legend=c("Cow", "Badger"),
       pch=c(17, 16), cex=1, col=c("blue", "red"), 
       text.col=c("blue", "red"), bty='n')
text(x=20, y=0, labels="AF2122/97")

# Add Scale bar
points(x=c(-20, 30), y=c(-130, -130), type="l", lwd=3)
text(x=5, y=-135, labels="50 SNPs", cex=1)

# Add Clade labels
text(x=92.5, y=-82, labels="0", col=cladeColours[1], cex=2)
text(x=90, y=80, labels="1", col=cladeColours[2], cex=2)
text(x=19.5, y=128.5, labels="2", col=cladeColours[3], cex=2)
text(x=-82, y=-95, labels="3", col=cladeColours[4], cex=2)

################################
# Get the sampling information #
################################

# Read in the badger sampling information
fileName <- paste(path, "IsolateData/", "BadgerInfo_08-04-15_LatLongs_XY_Centroids.csv",
                  sep="")
metadata <- read.table(fileName, header=TRUE, stringsAsFactors=FALSE, sep=",")

# Get the locations of each of the isolates
badgerIsolateLocations <- noteBadgerIsolateSamplingLocations(metadata)

# Cattle Isolates
file <- paste(path, "IsolateData/",
              "CattleIsolateInfo_LatLongs_plusID_outbreakSize_Coverage_AddedStrainIDs.csv", sep="")
cattleInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Get the locations of each of the isolates
cattleIsolateLocations <- noteCattleIsolateSamplingLocations(cattleInfo)

####################################################
# Plot the spatial locations of isolates in clades #
####################################################

# Create the clade colours - apply alpha
cladeColoursRGB <- getRGBsOfColours(cladeColours, alpha=0.75)
cex=2

# Note the isolates in each clade
isolatesInClades <- findIsolatesInClades(tree, nodesDefiningClades)

# Note the centre of the badger territories
badgerCentre <- c(381761.7, 200964.3)
expand <- 7000

# Create an empty plot
par(mar=c(0,0,0,0))
plot(x=NULL, y=NULL, yaxt="n", xaxt="n", bty="n", ylab="",
     xlim=c(badgerCentre[1] - expand, badgerCentre[1] + expand), 
     ylim=c(badgerCentre[2] - expand, badgerCentre[2] + expand), asp=1,
     xlab="")

# Plot a minimum convex polygon around the 
# cattle and badger sampling locations for each cluster
for(i in 1:length(cladeColours)){
  
  # Get the isolates associated with the current clade
  isolates <- isolatesInClades[[as.character(i)]]
  
  # Get the coordinates of each isolate
  isolateCoordinates <- getXandYCoordinatesOfIsolates(isolates, cattleIsolateLocations,
                                                      badgerIsolateLocations)
  
  # Remove NA rows - where couldn't find coordinates for isolates
  isolateCoordinates <- isolateCoordinates[is.na(isolateCoordinates$X) == FALSE, ]
  
  # Plot the points
  points(isolateCoordinates, 
         pch=ifelse(isolateCoordinates$Species == "BADGER", 19, 17),
         col=cladeColoursRGB[i], cex=cex)
  
  # Add a convex hull around the points
  addPolygon(isolateCoordinates$X, isolateCoordinates$Y, cladeColours[i])
}

# Add inner circle from BASTA deme assignment diagram
thresholdDistance <- 3500
draw.circle(x=badgerCentre[1], y=badgerCentre[2], radius=thresholdDistance,
            border="black", lty=2)
text(x=badgerCentre[1], y=badgerCentre[2] - (thresholdDistance + 500),
     labels=paste(round(thresholdDistance/1000, digits=2), "km radius"))

# Add legend
legend("bottomleft", legend=c("CATTLE", "BADGERS"),
       pch=c(17, 16), col="black", pt.cex=cex,
       text.col="black", bty='n')

# Add the cluster numbers
legend("bottomright", legend=addTextToArray("Cluster ", 0:3, ""),
       text.col=cladeColours, bty="n", cex=2)

dev.off()

##########################################################
# Print file noting which isolates are in which clusters #
##########################################################

# Note the clades of isolates in clades
isolateClades <- noteCladesOfIsolates(tree, nodesDefiningClades)

# Print out table
file <- paste(path, "vcfFiles/", "clusters_02-10-17.csv", sep="")
write.table(isolateClades, file, quote=FALSE, sep=",", row.names=FALSE)

#############
# FUNCTIONS #
#############

noteCladesOfIsolates <- function(tree, nodesDefiningClades){
  
  # Initialise two arrays to store the isolate IDs and clades
  isolates <- c()
  clades <- c()

  # Examine each clade
  for(i in 1:length(nodesDefiningClades)){
    tipsInClade <- tips(tree, nodesDefiningClades[i])
    
    for(tip in tipsInClade){
      isolates[length(isolates) + 1] <- tip
      clades[length(clades) + 1] <- i - 1
    }
  }
  
  # Combine the arrays into table
  output <- data.frame(ID=isolates, Cluster=clades, stringsAsFactors=FALSE)
  
  return(output)
}

addTextToArray <- function(text, array, sep){
  
  output <- c()
  for(i in 1:length(array)){
    output <- paste(text, array, sep=sep)
  }
  return(output)
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

findIsolatesInClades <- function(tree, nodesDefiningClades){
  isolatesInClades <- list()
  for(i in 1:length(nodesDefiningClades)){
    isolatesInClades[[as.character(i)]] <- tips(tree, nodesDefiningClades[i])
  }
  
  return(isolatesInClades)
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
        if(grepl(pattern="TB", x=tree$tip.label[tipIndex]) == TRUE){
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

addLegendForQuality <- function(position, cex){

  sizes <- seq(0.6, 1, 0.05)
  
  legend(position, legend=sizes, col="black", pch=24, bty='n',
         pt.cex=sizes, cex=cex)
}

getIsolateQuality <- function(table){
  isolateQuality <- list()
  for(i in 1:nrow(table)){
    isolateQuality[[table[i, "Isolate"]]] <- table[i, "Coverage"]
  }
  
  return(isolateQuality)
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

getUpperBoundsFromCuts <- function(cuts){
  
  bounds <- c()
  for(i in 1:length(cuts)){
    bounds[i] <- as.numeric(strsplit(strsplit(cuts[i], ",")[[1]][2], "]")[[1]][1])
  }
  
  return(bounds)
}

addLegend <- function(position, colours, nBreaks, cex){
  
  colourPalette <- colorRampPalette(colours)
  
  cuts <- levels(cut(table$PercentageCoverage, breaks=nBreaks))
  
  bounds <- getUpperBoundsFromCuts(cuts)
  
  bounds <- round(bounds, 2)
  
  legend(position, legend=bounds, col=colourPalette(nBreaks), pch=20, bty='n',
         cex=cex)
  
}

getIsolateIDFromFileNames <- function(fileNames){
  isolates <- c()
  for(i in 1:length(fileNames)){
    isolates[i] <- strsplit(fileNames[i], split="_")[[1]][1]
  }
  
  return(isolates)
}

assignIsolatesContinuousColoursByCoverage <- function(table, colours, nBreaks){
  
  colourPalette <- colorRampPalette(colours)
  coloursPerRow <- colourPalette(nBreaks)[
    as.numeric(cut(table$PercentageCoverage, breaks=nBreaks))]
  
  isolateColours <- list()
  for(i in 1:nrow(table)){
    isolateColours[[table[i, 1]]] <- coloursPerRow[i]
  }
  
  return(isolateColours)
}

returnTipColoursForIsolates <- function(tipLabels, assignedColours){
  
  tipColours <- c()
  for(i in 1:length(tipLabels)){
    if(tipLabels[i] != "Ref-1997"){
      
      tipColours[i] <- assignedColours[[tipLabels[i]]]
    }else{
      tipColours[i] <- "black"
    }
  }
  
  return(tipColours)
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
