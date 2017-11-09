#########################
# Read in sequence data #
#########################

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_13-07-17/"

# Read in the FASTA file
file <- paste(path, "vcfFiles/sequences_Prox-10_29-09-2017.fasta", sep="")
sequences <- readFasta(file)

###############################
# Calculate Genetic Distances #
###############################

geneticDistances <- buildGeneticDistanceMatrix(sequences)

####################################
# Get the location of each isolate #
####################################

# Read in the badger isolate information
fileName <- paste(path, "IsolateData/", "BadgerInfo_08-04-15_LatLongs_XY_Centroids.csv",
                  sep="")
badgerInfo <- read.table(fileName, header=TRUE, stringsAsFactors=FALSE, sep=",")

# Read in the cattle isolate information
file <- paste(path, "IsolateData/", 
              "CattleIsolateInfo_LatLongs_plusID_outbreakSize_Coverage_AddedStrainIDs.csv", sep="")
cattleInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Note the X and Y coordinates of each isolate and put in list
isolateLocations <- noteIsolateLocations(badgerInfo, cattleInfo)

#######################################
# Get the isolate cluster assignments #
#######################################

# Read in the cluster assignment file
fileName <- paste(path, "vcfFiles/", "clusters_02-10-17.csv",
                  sep="")
clusters <- read.table(fileName, header=TRUE, stringsAsFactors=FALSE, sep=",")

# Create a list to note each isolates cluster
isolateClusters <- noteIsolateClusters(clusters)


###################################################################################################
# Calculate mean genetic distance to badger isolates and spatial distance to territories centroid #
###################################################################################################

# Note the badger territory centroid
badgerCentre <- c(381761.7, 200964.3) # Taken from PlotIsolateLocationsForBASTADemeAssignment_24-05-17.R

# Create a table containing the mean genetic distances to badgers and 
# spatial distances to the badger territory centroids
distances <- calculateMeanGeneticDistanceAndSpatialDistanceToBadgers(sequences,
                                                                     geneticDistances,
                                                                     isolateLocations)

# Add species column
distances$Species = "COW"
distances$Species[grepl(distances$Isolate, pattern="WB")] <- "BADGER"

# Add colour column based upon clusters
cladeColours <- c("cyan", "pink", "green", "darkorchid4") # Taken from PlotMLTreeIncPoorAndRef_13-10-16.R
distances <- addIsolateClusterColours(distances, isolateClusters, cladeColours)

#################
# Create a plot #
#################

# Get current date
date <- format(Sys.time(), "%d-%m-%y")

# Open a pdf
file <- paste(path, "SpatialAndGeneticDistanceToBadgers_", date, ".pdf", sep="")
pdf(file)

plot(distances$MeanGeneticDistanceToBadgers, 
     distances$SpatialDistanceToCentroid / 1000,
     las=1, 
     xlab="Mean Genetic Distance to Badger Isolates (SNPs)",
     ylab="Spatial Distance (km)",
     main="Comparing Spatial Distance to Woodchester Park\nversus Genetic Distance to Badger Isolates",
     pch=20,
     col=ifelse(distances$Species == "BADGER", rgb(1,0,0, 0.5), rgb(0,0,1, 0.5)))
legend("topleft", legend=c("BADGER", "COW"), text.col=c("red", "blue"), bty="n")

subset <- distances[distances$MeanGeneticDistanceToBadgers < 40, ]

plot(subset$MeanGeneticDistanceToBadgers, 
     subset$SpatialDistanceToCentroid / 1000,
     las=1, 
     xlab="Mean Genetic Distance to Badger Isolates (SNPs)",
     ylab="Spatial Distance (km)",
     main="Comparing Spatial Distance to Woodchester Park\nversus Genetic Distance to Badger Isolates",
     pch=20,
     col=ifelse(subset$Species == "BADGER", rgb(1,0,0, 0.5), rgb(0,0,1, 0.5)))
legend("topleft", legend=c("BADGER", "COW"), text.col=c("red", "blue"), bty="n")

plot(subset$MeanGeneticDistanceToBadgers, 
     subset$SpatialDistanceToCentroid / 1000,
     las=1, 
     xlab="Mean Genetic Distance to Badger Isolates (SNPs)",
     ylab="Spatial Distance (km)",
     main="Comparing Spatial Distance to Woodchester Park\nversus Genetic Distance to Badger Isolates",
     pch=ifelse(subset$Species == "BADGER", 19, 17),
     col=subset$ClusterColour)
legend("topleft", legend=paste("Cluster ", 0:3, ""),
       text.col=cladeColours, bty="n", cex=0.75)
legend("topright", legend=c("BADGER", "COW"), text.col="black", bty="n", pch=c(19, 17),
       cex=0.75)


dev.off()

#############
# FUNCTIONS #
#############

addIsolateClusterColours <- function(distances, isolateClusters, clusterColours){
  distances$ClusterColour <- rgb(0,0,0, 0.5)
  for(row in 1:nrow(distances)){
    
    if(is.null(isolateClusters[[distances[row, "Isolate"]]]) == FALSE){
      
      distances[row, "ClusterColour"] <- setAlpha(
        clusterColours[isolateClusters[[distances[row, "Isolate"]]] + 1],
        alpha=0.5)
    }
  }
  
  return(distances)
}

setAlpha <- function(colour, alpha){
  
  rgbValues <- col2rgb(colour)
  
  # Note that col2rgb returns rgbvlues from 0 to 255
  rgbColour <- rgb(rgbValues["red", 1],
                   rgbValues["green", 1],
                   rgbValues["blue", 1], 
                   alpha=alpha*255, 
                   maxColorValue=255)
  
  return(rgbColour)
}

noteIsolateClusters <- function(){
  isolateClusters <- list()
  for(row in 1:nrow(clusters)){
    isolateClusters[[clusters[row, "ID"]]] <- clusters[row, "Cluster"]
  }
  
  return(isolateClusters)
}

calculateMeanGeneticDistanceAndSpatialDistanceToBadgers <- function(sequences, 
                                                                    geneticDistances,
                                                                    isolateLocations){
  # Get a vector of the isolate IDs
  isolates <- names(sequences)
  isolates <- isolates[isolates != "Ref-1997"] # Remove reference
  
  # Initialise a dataframe to store mean genetic distance and spatial distance
  distances <- data.frame(Isolate=isolates, 
                          MeanGeneticDistanceToBadgers=rep(NA, length(isolates)),
                          SpatialDistanceToCentroid=rep(NA, length(isolates)),
                          stringsAsFactors=FALSE)
  
  # Note the indices of badger isolates
  badgerIndices <- which(grepl(isolates, pattern="WB"))
  
  # Examine each isolate
  for(i in 1:length(isolates)){
    
    # Remove badger Index is current isolate is a badger
    badgers <- badgerIndices[badgerIndices != i]
    
    # Calculate the mean genetic distance of the current isolate to all badgers
    distances[i, "MeanGeneticDistanceToBadgers"] <- 
      mean(geneticDistances[i, badgers], na.rm=TRUE)
    
    # Calculate the spatial distance of current isolate to territory centroid
    if(is.na(isolateLocations[[isolates[i]]][1]) == FALSE){
      distances[i, "SpatialDistanceToCentroid"] <- 
        euclideanDistance(badgerCentre[1], badgerCentre[2], 
                          isolateLocations[[isolates[i]]][1],
                          isolateLocations[[isolates[i]]][2])
    }
  }
  
  return(distances)
}

euclideanDistance <- function(x1, y1, x2, y2){
  return(sqrt(sum((x1 - x2)^2 + (y1 - y2)^2)))
}

noteIsolateLocations <- function(badgerInfo, cattleInfo){
  
  # Initialise a list to store the isolate coordinates
  isolateLocations <- list()
  
  # Add badger isolate locations
  for(row in 1:nrow(badgerInfo)){
    if(is.na(badgerInfo[row, "GroupCentroidX"]) == FALSE){
      isolateLocations[[badgerInfo[row, "WB_id"]]] <- c(
        badgerInfo[row, "GroupCentroidX"],
        badgerInfo[row, "GroupCentroidY"]
      )
    }else{
      isolateLocations[[badgerInfo[row, "WB_id"]]] <- c(
        badgerInfo[row, "SampledGrpX"],
        badgerInfo[row, "SampledGrpY"]
      )
    }
  }
  
  # Add cattle isolate locations
  for(row in 1:nrow(cattleInfo)){
    isolateLocations[[cattleInfo[row, "StrainId"]]] <- c(
      cattleInfo[row, "Mapx"],
      cattleInfo[row, "Mapy"]
    )
  }
  
  return(isolateLocations)
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
    
    if(i %% 10 == 0){
      cat(paste("Finished calculate distances for isolate:", i, "\n"))
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
