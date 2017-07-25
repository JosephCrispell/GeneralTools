# Load the ape package
library(ape)
library(geiger) # For the tips function

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/"

###################################
# Get the Maximum Likelihood Tree #
###################################

# Read in the newick tree
file <- paste(path, "allVCFs-IncludingPoor/vcfFiles/",
              "mlTree_Prox-10_plusRef_rmResequenced_SNPCov-0.1_28-10-16.tree", sep="")
tree <- read.tree(file=file)

# Convert Branch lengths to SNPs
fastaLength <- 9464
tree$edge.length <- tree$edge.length * fastaLength

##################################
# Get the Isolate Coverage Table #
##################################

file <- paste(path,"allVCFs-IncludingPoor/vcfFiles/",
              "isolateGenomeCoverageSummary_28-10-16.txt", sep="")
table <- read.table(file, header=TRUE, stringsAsFactors=FALSE)

table$IsolateID <- getIsolateIDFromFileNames(table$IsolateID)

##############################
# Plot the Phylogenetic Tree #
##############################

file <- paste(path, "mlTree_CladesAndLocations_22-06-17.pdf")
pdf(file, height=10, width=10)

# Set the margins
par(mar=c(0,0,0,0)) # Bottom, Left, Top, Right

plotType <- "fan" # "phylogram", "cladogram", "fan", "unrooted", "radial"

# Plot initial tree to find nodes defining clades
#plot.phylo(tree, plotType)
#nodelabels()

# Define branch colours by clade
#nodesDefiningClades <- c(515, 318, 521, 405, 433) # use nodelabels() to show node numbers
nodesDefiningClades <- c(291, 440, 460, 412, 324) # use nodelabels() to show node numbers
cladeColours <- c("cyan", "pink", "green", "orange", "purple")
branchColours <- defineBranchColoursOfClades(tree, nodesDefiningClades, cladeColours)

# Get each isolate's quality
isolateQuality <- getIsolateQuality(table)

# Plot the phylogenetic tree
plot.phylo(tree, show.tip.label=FALSE, plotType,
           edge.color=branchColours, edge.width=3,
           show.node.label=TRUE)

# Add node labels
nodelabels(node=1:length(tree$tip.label), 
           cex=defineTipSizeBySequencingQuality(tree$tip.label, isolateQuality),
           pch=defineTipShapesForSpecies(tree$tip.label, 17, 16),
           col=defineTipColourBySpecies(tree$tip.label, "blue", "red"))

# Add Legends
text(x=132, y=-84, labels="Coverage:", col="black", cex=0.7)
addLegendForQuality("bottomright", 0.8)
text(x=-113, y=-120, labels="Species:", col="black", cex=0.7)
legend("bottomleft", legend=c("CATTLE", "BADGERS"),
       pch=c(17, 16), cex=0.65, col=c("blue", "red"), 
       text.col=c("blue", "red"), bty='n')

# Add Scale bar
points(x=c(-20, 30), y=c(-114, -114), type="l", lwd=3)
text(x=5, y=-118, labels="50 SNPs", cex=0.8)

# Add Clade labels
text(x=52, y=102, labels="0", col=cladeColours[1], cex=2)
text(x=55, y=-100, labels="1", col=cladeColours[2], cex=2)
text(x=95, y=-75, labels="2", col=cladeColours[3], cex=2)
text(x=-45, y=-105, labels="3", col=cladeColours[4], cex=2)
text(x=-102, y=50, labels="4", col=cladeColours[5], cex=2)

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
              "CattleIsolateInfo_LatLongs_plusID_outbreakSize_Coverage_AddedTB1453-TB1456.csv", sep="")
cattleInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Get the locations of each of the isolates
cattleIsolateLocations <- noteCattleIsolateSamplingLocations(cattleInfo)

####################################################
# Plot the spatial locations of isolates in clades #
####################################################

# Create the clade colours - apply alpha
cladeColours <- c("cyan", "pink", "green", "orange", "purple")
cladeColoursRGB <- getRGBsOfColours(cladeColours, alpha=0.75)
cex=2

# Note the isolates in each clade
isolatesInClades <- findIsolatesInClades(tree, nodesDefiningClades)

# Note the centre of the badger territories
badgerCentre <- c(381761.7, 200964.3)
expand <- 6500

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

# Add legend
legend("bottomleft", legend=c("CATTLE", "BADGERS"),
       pch=c(17, 16), col="black", pt.cex=cex,
       text.col="black", bty='n')

# Add the cluster numbers
legend("bottomright", legend=addTextToArray("Cluster ", 0:4, ""),
       text.col=cladeColours, bty="n", cex=2)

# Add Scale
legend("bottom", legend=(paste(round(expand/1000, digits=2), "KM")), bty="n")

######
######
######

# Create the clade colours - apply alpha
cladeColours <- c("cyan", "pink", "green", "orange", "purple")
cladeColoursRGB <- getRGBsOfColours(cladeColours, alpha=0.75)
cex=2

# Note the isolates in each clade
isolatesInClades <- findIsolatesInClades(tree, nodesDefiningClades)

# Note the centre of the badger territories
badgerCentre <- c(381761.7, 200964.3)
expand <- 6500

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
thresholdDistance <- 3000
draw.circle(x=badgerCentre[1], y=badgerCentre[2], radius=thresholdDistance,
            border="black")

# Add legend
legend("bottomleft", legend=c("CATTLE", "BADGERS"),
       pch=c(17, 16), col="black", pt.cex=cex,
       text.col="black", bty='n')

# Add the cluster numbers
legend("bottomright", legend=addTextToArray("Cluster ", 0:4, ""),
       text.col=cladeColours, bty="n", cex=2)

# Add Scale
legend("bottom", legend=(paste(round(expand/1000, digits=2), "KM")), bty="n")

dev.off()

#############
# FUNCTIONS #
#############

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

defineTipColourBySpecies <- function(tipLabels, cow, badger){
  colours <- c()
  for(i in 1:length(tipLabels)){
    
    if(grepl(pattern="TB", x=tipLabels[i]) == TRUE){
      colours[i] <- cow
    }else if(grepl(pattern="WB", x=tipLabels[i]) == TRUE){
      colours[i] <- badger
    }else{
      colours[i] <- "black"
    }
  }
  
  return(colours)
}

addLegendForQuality <- function(position, cex){

  sizes <- seq(0.1, 1, 0.1)
  
  legend(position, legend=sizes, col="black", pch=16, bty='n',
         pt.cex=sizes, cex=cex)
}

getIsolateQuality <- function(table){
  isolateQuality <- list()
  for(i in 1:nrow(table)){
    isolateQuality[[table[i, "IsolateID"]]] <- table[i, "PercentageCoverage"]
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
                                        CladeColours){
  branchColours <- rep("black", dim(tree$edge)[1])
  for(i in 1:length(cladeColours)){
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
