#### Preparation ####

# Libraries
library(ape) # Read and plot phylogeny
library(phytools) # For getDescendants function
library(basicPlotteR) # For setting alphas

# Set the path
path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/"

# Get today's date
date <- format(Sys.Date(), "%d-%m-%y")

#### Read in the isolate data ####

# Read in the cattle data
cattleIsolateFile <- paste(path, "IsolateData/CattleIsolateInfo_AddedNew_TB1484-TB1490_22-03-18.csv", sep="")
cattleInfo <- read.table(cattleIsolateFile, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Read in the badger data
badgerIsolateFile <- paste(path, "IsolateData/BadgerInfo_08-04-15_LatLongs_XY_Centroids.csv", sep="")
badgerInfo <- read.table(badgerIsolateFile, header=TRUE, stringsAsFactors=FALSE, sep=",")

#### Read in BASTA phylogeny and plot ####

# Read in the BASTA phylogeny
treeFile <- paste0(path, "vcfFiles/mlTree_BASTAClade_DatedTips_27-03-18.tree")
tree <- read.tree(treeFile)

# Remove date from tree tips
tree$tip.label <- parseLabels(tree$tip.label)

# Take a look at the tree
#viewRAxMLTree(bastaClade, path)

# Plot the phylogeny
nodesDefiningClades <- c(292, 316, 187)
cladeColours <- c("magenta", "green", "cyan")
file <- paste0(path, "ESM_Figures/BASTAPhylogeny_", date, ".pdf")
plotTree(tree, nodes=nodesDefiningClades, colours=cladeColours, file=file)

#### Plot the tips in space ####

file <- paste0(path, "ESM_Figures/BASTATipLocations_", date, ".pdf")
plotIsolateLocations(tree, badgerInfo=badgerInfo, cattleInfo=cattleInfo, colours=cladeColours, nodes=nodesDefiningClades,
                     badgerCentre=c(381761.7, 200964.3), expand=7000, file=file, alpha=0.5)


#### Plot the lifespans of the animals in each clade ####

# Note the cluster interested in
cluster <- 4

# Read in the sampled animal lifespans
file <- paste(path, "InterSpeciesClusters/sampledAnimalsLifeHistories_22-11-2018.txt", sep="")
table <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep="\t")
table <- table[is.na(table$Isolates) == FALSE & grepl(table$Clusters, pattern=cluster - 1), ]
table <- removeAnimalsNotAssociatedWithBASTAClade(table, tree)

# Plot the lifespan of each animal in cluster
file <- paste0(path, "ESM_Figures/Clade4_Lifespans_", date, ".pdf")
plotAnimalLifespans(table, alpha=0.25, lwd=4, file=file, axisCex=1.5)

#### Plot the network relationships between the animals ####

file <- paste0(path, "ESM_Figures/Clade4_Movements_", date, ".pdf")
plotAnimalMovements(table, badgerCentre=c(381761.7, 200964.3), expand=2000, lwd=4, lineAlpha=0.25, file=file)


#### FUNCTIONS ####

plotAnimalMovements <- function(animalInfo, badgerCentre, expand, lineAlpha=1, file=NULL, ...){
  
  # Open a pdf file if requested
  if(is.null(file) == FALSE){
    pdf(file)
  }
  
  # Get and set the plotting margins
  currentMar <- par()$mar
  par(mar=c(0,0,0,0))
  
  # Create an empty plot
  plot(x=NULL, y=NULL, yaxt="n", xaxt="n", bty="n", ylab="",
       xlim=c(badgerCentre[1] - expand, badgerCentre[1] + expand), 
       ylim=c(badgerCentre[2] - expand, badgerCentre[2] + expand), asp=1,
       xlab="")
  
  # Examine every animal
  for(row in 1:nrow(animalInfo)){
    
    # Get the X coordinates
    xCoords <- as.numeric(strsplit(animalInfo[row, "Xs"], split=",")[[1]])
        
    # Get the Y coordinates
    yCoords <- as.numeric(strsplit(animalInfo[row, "Ys"], split=",")[[1]])
    
    # Remove NA locations
    NAs <- which(is.na(xCoords) | is.na(yCoords))
    if(length(NAs) > 0){
      xCoords <- xCoords[-NAs]
      yCoords <- yCoords[-NAs]
    }
    
    # Remove sequential duplicate locations
    indices <- identifySequentialDuplicateLocations(xCoords, yCoords)
    xCoords <- xCoords[indices]
    yCoords <- yCoords[indices]
    
    # Round locations?
    xCoords <- round(xCoords, digits=0)
    yCoords <- round(yCoords, digits=0)
    
    # Plot each location of the movements
    points(x=xCoords, y=yCoords, 
           pch=ifelse(animalInfo[row, "Species"] == "BADGER", 19, 17),
           col=ifelse(animalInfo[row, "Species"] == "BADGER", "red", "blue"), 
           cex=2)
    
    # Plot movements as arrows
    if(length(xCoords) > 1){

      # Plot each movement
      for(i in 2:length(xCoords)){
        points(x=c(xCoords[i-1], xCoords[i]), y=c(yCoords[i-1], yCoords[i]), type="l",
               col=ifelse(animalInfo[row, "Species"] == "BADGER", rgb(1,0,0, lineAlpha), rgb(0,0,1, lineAlpha)), ...)
      }
    }
  }
  
  # Add scale
  axisLimits <- par()$usr
  xLength <- axisLimits[2] - axisLimits[1]
  yLength <- axisLimits[4] - axisLimits[3]
  xPad <- 0.08*xLength
  points(x=c(axisLimits[2] - xPad - 1000, axisLimits[2] - xPad), y=c(axisLimits[3] + 0.1*yLength, axisLimits[3] + 0.1*yLength),
         type="l", lwd=4)
  text(x=axisLimits[2] - xPad - 500, y=axisLimits[3] + 0.07*yLength, labels="1 KM", cex=2)
  
  # Add Legend
  legend("bottom", legend=c("Cow", "Badger"), text.col=c("blue", "red"), pch=c(17, 19), cex=2, bty="n", col=c("blue", "red"))
  
  # Reset the plotting margins
  par(mar=currentMar)
  
  # Open a pdf file if requested
  if(is.null(file) == FALSE){
    dev.off()
  }
}

identifySequentialDuplicateLocations <- function(xCoords, yCoords){
  
  # Initialise a vector to store the indices of the unique locations
  indices <- c(1)
  
  # Initialise variables to store previous 
  previous <- c(xCoords[1], yCoords[2])
  
  # Check if more than one location
  if(length(xCoords) > 1){
    
    # Examine the locations
    for(i in 2:length(xCoords)){
      
      # Compare the current location to the previous
      if(xCoords[i] != previous[1] || yCoords[i] != previous[2]){
        indices[length(indices) + 1] <- i
        previous[1] = xCoords[i]
        previous[2] = yCoords[i]
      }
    }
  }
  
  return(indices)
}

plotAnimalLifespans <- function(animalInfo, alpha=0.5, axisCex=1, file=NULL, ...){
  
  # Open a pdf file if requested
  if(is.null(file) == FALSE){
    pdf(file)
  }
  
  # Remove animals without movement dates
  animalInfo <- animalInfo[is.na(animalInfo$MovementDates) == FALSE, ]
  
  # Get and set the plotting margins
  currentMar <- par()$mar
  par(mar=c(3, 0, 0, 0))
  
  # Note earliest and latest dates to be considered
  dateRange <- getDateRange(animalInfo)
  
  # Convert the detection date column to a date
  animalInfo$DetectionDate <- as.Date(animalInfo$DetectionDate, format="%d-%m-%Y")
  
  # Order the animals by their detection date
  animalInfo <- animalInfo[order(animalInfo$DetectionDate), ]
  
  # Create an empty plot
  plot(x=NULL, y=NULL, xlim=dateRange, ylim=c(1, nrow(animalInfo)), bty="n", xaxt="n", xlab="")
  
  # Add X axis
  at <- as.Date(c("1995-01-01", "2000-01-01", "2005-01-01", "2010-01-01"), format="%Y-%m-%d")
  axis(side=1, at=at, labels=format(at, "%Y"), cex.axis=axisCex)
  
  # Add a line for every animal
  for(row in 1:nrow(animalInfo)){

    # Get the movement dates for the current animal
    dates <- strsplit(animalInfo[row, "MovementDates"], split=",")[[1]]
    dates <- as.Date(dates, format="%d-%m-%Y")
      
    # Plot the negative period of life
    negativePeriod <- c(dates[1], animalInfo[row, "DetectionDate"])
    points(x=negativePeriod, y=c(row, row), type="l", 
           col=ifelse(animalInfo[row, "Species"] == "BADGER", rgb(1,0,0, alpha), rgb(0,0,1, alpha)), ...)
    
    # Plot the infected period of life
    infectedPeriod <- c(animalInfo[row, "DetectionDate"], dates[length(dates)])
    points(x=infectedPeriod, y=c(row, row), type="l", 
           col=ifelse(animalInfo[row, "Species"] == "BADGER", rgb(1,0,0), rgb(0,0,1)), ...)
  }
  
  # Add Legend
  legend("bottomright", legend=c("Cow", "Badger"), text.col=c("blue", "red"), cex=2, bty="n")
  
  # Open a pdf file if requested
  if(is.null(file) == FALSE){
    dev.off()
  }
  
  # Reset the plotting margins
  par(mar=currentMar)
}

getDateRange <- function(animalInfo){
  
  # Initialise a start and end date
  start <- as.Date("01-01-2020", format="%d-%m-%Y")
  end <- as.Date("01-01-1900", format="%d-%m-%Y")
  
  # Examine the movement dates associated with every animal
  for(row in seq_len(nrow(animalInfo))){
    
    # Skip if no movement dates available
    if(is.na(animalInfo[row, "MovementDates"])){
      next
    }
    
    # Get the movement dates for the current animal
    dates <- strsplit(animalInfo[row, "MovementDates"], split=",")[[1]]
    dates <- as.Date(dates, format="%d-%m-%Y")
    
    # Check if first date before current start
    if(dates[1] < start){
      start <- dates[1]
    }
    
    # Check if last date after current end
    if(dates[length(dates)] > end){
      end <- dates[length(dates)]
    }
  }
  
  return(c(start, end))
}

removeAnimalsNotAssociatedWithBASTAClade <- function(clusterInfo, tree){
  
  remove <- c()
  for(row in 1:nrow(table)){
    
    # Get the isolates associated with the current animal
    isolates <- strsplit(table[row, "Isolates"], split=",")[[1]]
    
    # Count how many isolates were in BASTA analyses
    nFound <- length(which(isolates %in% tree$tip.label))
    
    # Ignore animals if no isolates found
    if(nFound < 1){
      remove[length(remove) + 1] <- row
    }
  }
  table <- table[-remove, ]
  
  return(table)
}

parseLabels <- function(tipLabels){
  
  output <- c()
  for(i in seq_along(tipLabels)){
    output[i] <- strsplit(tipLabels[i], split="_")[[1]][1]
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

plotIsolateLocations <- function(tree, badgerInfo, cattleInfo, colours, badgerCentre, expand, nodes, file=NULL, alpha=0.75){
  
  # Open a pdf file if requested
  if(is.null(file) == FALSE){
    pdf(file)
  }
  
  # Read in the badger sampling information
  badgerIsolateLocations <- noteBadgerIsolateSamplingLocations(badgerInfo)
  
  # Cattle Isolates
  cattleIsolateLocations <- noteCattleIsolateSamplingLocations(cattleInfo)
  
  # Create the clade colours - apply alpha
  cladeColoursRGB <- setAlpha(colours, alpha=alpha)
  
  # Note the isolates in each clade
  isolatesInClades <- findIsolatesInClades(tree, nodes)
  
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
    points(x=isolateCoordinates$X, y=isolateCoordinates$Y, 
           pch=ifelse(isolateCoordinates$Species == "BADGER", 19, 17),
           col=cladeColoursRGB[i], cex=2)
    
    # Add a convex hull around the points
    addPolygon(isolateCoordinates$X, isolateCoordinates$Y, colours[i])
  }
  
  # Add legend
  legend("bottomleft", legend=c("Cow", "Badger"),
         pch=c(17, 16), col="black", pt.cex=2,
         text.col="black", bty='n')
  
  # Open a pdf file if requested
  if(is.null(file) == FALSE){
    dev.off()
  }
}

addPolygon <- function(xValues, yValues, borderColour){
  hullPoints <- chull(xValues, yValues)
  hullPoints <- c(hullPoints, hullPoints[1])
  polygon(xValues[hullPoints], yValues[hullPoints], col = NA, border = borderColour)
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

plotTree <- function(tree, nodes, colours, scaleSize=5, file=NULL){
  
  # Open a pdf file if requested
  if(is.null(file) == FALSE){
    pdf(file)
  }
  
  # Get and set the margins
  currentMar <- par()$mar
  par(mar=c(0, 0, 0, 0))
  
  # Define the branch colours - branches in clades are coloured by clade
  branchColours <- defineBranchColoursOfClades(tree, nodes,
                                               colours, "lightgrey")
  
  # Define the characteristics of the tips based on species
  tipShapes <- defineTipShapesForSpecies(tree$tip.label, 24, 21)
  tipColour <- defineTipColourBySpecies(tree, "blue", "red", "lightgrey",
                                        nodes)
  
  # Plot the phylogenetic tree
  plot.phylo(tree, show.tip.label=FALSE, "fan",
             edge.color=branchColours, edge.width=3)
  
  # Add node labels
  tiplabels(pch=tipShapes, bg=tipColour, col="dimgrey")
  
  # Add Legends
  legend("bottomleft", legend=c("Cow", "Badger"),
         pch=c(17, 16), cex=2, col=c("blue", "red"), 
         text.col=c("blue", "red"), bty='n')

  # Add Scale bar
  axisLimits <- par()$usr
  xLength <- axisLimits[2] - axisLimits[1]
  yLength <- axisLimits[4] - axisLimits[3]
  xPad <- 0.05*xLength
  yPad <- 0.1*yLength
  points(x=c(axisLimits[2] - xPad, axisLimits[2] - xPad - scaleSize), 
         y=c(axisLimits[3] + 0.15*yLength, axisLimits[3] + 0.15*yLength), 
         type="l", lwd=4, xpd=TRUE)
  text(x=axisLimits[2] - xPad - (0.5*scaleSize), y=axisLimits[3] + 0.12*yLength, cex=2, xpd=TRUE,
       labels=ifelse(scaleSize > 1, paste0(scaleSize, " SNPs"), paste0(scaleSize, " SNP")))
  
  # Reset margins
  par(mar=currentMar)
  
  # Close the pdf if requested
  if(is.null(file) == FALSE){
    dev.off()
  }
}
