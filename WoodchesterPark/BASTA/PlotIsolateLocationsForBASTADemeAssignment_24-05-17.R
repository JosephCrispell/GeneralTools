##################
# Load Libraries #
##################

library(maptools) # Read shape file
library(rgeos) # Polygon centroids
library(ape) # Phylogenetic tree
library(geiger) # Phylogenetic tree tools
library(plotrix) # Draw circle

# Create a path variable
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_13-07-17/"

#######################################################################
# Calculate the Territory Centroids of Each Social Group in Each Year #
#######################################################################

# Note the years that territory shape files are available for
years <- c(2000:2011)
shapeFileNames <- c("territories_2000.shp", 
                    "territories_2001_autumn.shp",
                    "territories_2002_spring.shp",
                    "territories_2003.shp",
                    "territories_2004.shp",
                    "territories_2005.shp",
                    "MCPs2006 corrected.shp",
                    "MCPs_07.shp",
                    "MCP_2008.shp",
                    "MCP_2009.shp",
                    "MCPs 2010.shp",
                    "MCP_2011.shp")

# Initialise list to record each social group's location in each year
groupsCentroidsPerYear <- initialiseTableToStoreTerritoryCentroids(years)
groupsRows <- list()
index <- 0

# Initialise a list to record the territories from each year
territoriesForEachYear <- list()

# Initialise two arrays to store the min and max X and Y values for the territories
min <- c(999999999999, 99999999999)
max <- c(0, 0)

# Examine the territories in each year
for(i in 1:length(years)){
  
  year <- years[i]
  shapeFileName <- shapeFileNames[i]
  
  # Read in the shape file
  file <- paste(path, "BadgerCaptureData/BadgerTerritoryMarkingData/",
                "Baitmarking ", year, "/", shapeFileName, sep="")
  territories <- readShapePoly(file) # Generates SpatialPolygonsDataFrame
  
  # Extract the polygon coordinates
  territoryCoords <- getPolygonCoords(territories)
  min <- updateMin(territoryCoords, min)
  max <- updateMax(territoryCoords, max)
  
  # Get the full social group names
  territoryIDs <- getSocialGroupNames(territories@data, year)
  
  # Assign polygons to their social group names
  socialGroupTerritories <- assignTerritoriesToSocialGroupNames(territoryIDs, territoryCoords)
  territoriesForEachYear[[as.character(year)]] <- socialGroupTerritories
  
  # Calculate the territory centroids - mean X and Y
  territoryCentroids <- calculateTerritoryCentroids(territoryCoords, territoryIDs)
  
    # Give an index to all groups present
  groupsRows <- indexGroups(names(territoryCentroids), groupsRows, index)
  
  # Note the centroids of each group in the current year
  groupsCentroidsPerYear <- addGroupCentroidsFromCurrentYear(groupsCentroidsPerYear, groupsRows,
                                                             territoryCentroids, year)
  index <- nrow(groupsCentroidsPerYear)
}

# Add social group names to territory centroids table
groupsCentroidsPerYear <- addSocialGroupNames(groupsCentroidsPerYear, groupsRows)

###################################################
# Read in the cattle farm land parcel information #
###################################################

# Read in the shape file
file <- paste(path, "LandParcelData/RPA_CLAD_ASI_CURRENT_SP_ST-SO_NE-SE/",
              "RPA_CLAD_ASI_CURRENT_SP_ST-SO_NE-SE_REDUCED.shp", sep="")
landParcels <- readShapePoly(file) # Generates SpatialPolygonsDataFrame

# Get the coordinates for each polygon
landParcelCoords <- getPolygonCoords(landParcels)
min <- landParcelCoords[["min"]]
max <- landParcelCoords[["max"]]

# Note the polygon IDs that are associated with each CPH
cphs <- noteLandParcelsAssociatedWithCPHs(landParcels@data)

######################################
# Note isolates in clade of interest #
######################################

# Read in the newick tree
file <- paste(path, "vcfFiles/",
              "mlTree_29-09-2017.tree", sep="")
tree <- read.tree(file=file)

# Get a list of the isolates in the clade
node <- 301
cladeTips <- tips(tree, node=node)

##########################################
# Get badger isolates sampling locations #
##########################################

# Read in the isolate badgerInfo
fileName <- paste(path, "IsolateData/", "BadgerInfo_08-04-15_LatLongs_XY_Centroids.csv",
                  sep="")
badgerInfo <- read.table(fileName, header=TRUE, stringsAsFactors=FALSE, sep=",")

# Keep only rows for the isolates
badgerInfo <- badgerInfo[badgerInfo$WB_id %in% cladeTips, ]

##########################################
# Get cattle isolates sampling locations #
##########################################

# Cattle Isolates
file <- paste(path, "IsolateData/", 
              "CattleIsolateInfo_LatLongs_plusID_outbreakSize_Coverage_AddedStrainIDs.csv", sep="")
cattleInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Select info only for cattle isolates in clade
cattleInfo <- cattleInfo[cattleInfo$StrainId %in% cladeTips, ]

# Get a list of sampled CPHs
sampledCPHs <- names(cphs)[names(cphs) %in% cattleInfo$CPH]

###########################
# Start creating the plot #
###########################

# Define inner circle radius
thresholdDistance <- 3500
outerDistance <- 8100

# Open a PDF
file <- paste(path, "BASTA/", "DemeAssignmentDiagram_02-10-17.pdf", sep="")
pdf(file)

# Find centre point of badger locations
badgerCentre <- findBadgerGroupsCentre(groupsCentroidsPerYear)
expand <- 8000

##############################################
# Plot the badger territories from all years #
##############################################

# Plot the badger territories from each year onto a single plot
plotTerritories(territoriesForEachYear, badgerCentre, expand, rgb(1,0,0, 0.1))

##############################################
# Plot the land parcels of the sampled herds #
##############################################

# Plot polygons from sampled herds
for(cph in sampledCPHs){
  
  for(key in cphs[[cph]]){
    polygon(landParcelCoords[[key]], border=rgb(0,0,1, 0.2))
  }
}

#################################
# Plot the sampled cattle herds #
#################################
for(row in 1:nrow(cattleInfo)){
  points(x=cattleInfo[row, "Mapx"], y=cattleInfo[row, "Mapy"], pch=17,
         col=rgb(0,0,1, 0.5), cex=1.5)
}

#########################################
# Add Badger sampling locations to plot #
#########################################

# Plot the badger isolate sampling locations
addBadgerIsolatesLocations(badgerInfo=badgerInfo, groupsRows=groupsRows, 
                           groupsCentroidsPerYear=groupsCentroidsPerYear,
                           col=rgb(1,0,0, 0.5), cex=2)

###############################################
# Add circles to define inner and outer demes #
###############################################

# Add deme defining circles
draw.circle(x=badgerCentre[1], y=badgerCentre[2], radius=thresholdDistance,
            border="black")
draw.circle(x=badgerCentre[1], y=badgerCentre[2], radius=outerDistance,
            border="black")
text(x=badgerCentre[1], y=(badgerCentre[2] - thresholdDistance) + 500,
     labels=paste("Inner: ", paste(thresholdDistance / 1000, "km", sep="")))
text(x=badgerCentre[1], y=(badgerCentre[2] - outerDistance) + 500,
     labels=paste("Outer: ", paste(outerDistance / 1000, "km", sep="")))

##############
# Add Legend #
##############

# Add legend
legend("bottomright", legend=c("Badger", "Cow"), text.col=c("red", "blue"),
       bty="n", pch=c(20, 17), pt.cex=c(2, 1.5), col=c("red", "blue"))


#####
#####
#####


###########################################
# Add colours to distinguish cattle herds #
###########################################


##############################################
# Plot the badger territories from all years #
##############################################

# Plot the badger territories from each year onto a single plot
plotTerritories(territoriesForEachYear, badgerCentre, expand, rgb(0,0,0, 0.2))

##############################################
# Plot the land parcels of the sampled herds #
##############################################

# Create colour for each sampled herd
sampledCPHColours <- assignColoursToEachSampledCPH(sampledCPHs)

# Plot polygons from sampled herds
alpha <- 0.8
for(cph in sampledCPHs){
  
  for(key in cphs[[cph]]){
    polygon(landParcelCoords[[key]], border=setAlpha(sampledCPHColours[[cph]], alpha))
  }
}

#################################
# Plot the sampled cattle herds #
#################################
alpha <- 0.9
for(row in 1:nrow(cattleInfo)){
  points(x=cattleInfo[row, "Mapx"], y=cattleInfo[row, "Mapy"], pch=17,
         col=setAlpha(sampledCPHColours[[cattleInfo[row, "CPH"]]], alpha), cex=1.5)
}

#########################################
# Add Badger sampling locations to plot #
#########################################

# Plot the badger isolate sampling locations
addBadgerIsolatesLocations(badgerInfo=badgerInfo, groupsRows=groupsRows, 
                           groupsCentroidsPerYear=groupsCentroidsPerYear,
                           col=rgb(0,0,0, 0.5), cex=2)

###############################################
# Add circles to define inner and outer demes #
###############################################

# Add deme defining circles
draw.circle(x=badgerCentre[1], y=badgerCentre[2], radius=thresholdDistance,
            border="black")
draw.circle(x=badgerCentre[1], y=badgerCentre[2], radius=outerDistance,
            border="black")
text(x=badgerCentre[1], y=(badgerCentre[2] - thresholdDistance) + 500,
     labels=paste("Inner: ", paste(thresholdDistance / 1000, "km", sep="")))
text(x=badgerCentre[1], y=(badgerCentre[2] - outerDistance) + 500,
     labels=paste("Outer: ", paste(outerDistance / 1000, "km", sep="")))

##############
# Add Legend #
##############

# Add legend
legend("bottomright", legend=c("Badger", "Cow"), text.col=c("black", "black"),
       bty="n", pch=c(20, 17), pt.cex=c(2, 1.5))


#####
#####
#####

#####################################
# Add colours to distinguish clades #
#####################################

# Note the clades
nodesDefiningClades <- c(521, 305, 332, 382) # use nodelabels() to show node numbers
cladeColours <- c("cyan", "pink", "green", "darkorchid4")
isolatesInClades <- findIsolatesInClades(tree, nodesDefiningClades)

# Note which CPHs are associated with which clusters
groupOrHerdClusters <- getSampledGroupsOrHerdsAssociatedWithClusters(isolatesInClades,
                                                                     badgerInfo,
                                                                     cattleInfo)

##############################################
# Plot the badger territories from all years #
##############################################

# Plot the badger territories from each year onto a single plot
cladeColoursRGB <- getRGBsOfColours(cladeColours, alpha=0.2)
plotTerritoriesWithClusterColours(territoriesForEachYear, badgerCentre, expand,
                                  rgb(0,0,0, 0.1), cladeColoursRGB, groupOrHerdClusters,
                                  lwd=1.5)

# Not the clade of each isolate
isolateClades <- invertIsolatesInClades(isolatesInClades)

##############################################
# Plot the land parcels of the sampled herds #
##############################################

# Plot polygons from sampled herds
cladeColoursRGB <- getRGBsOfColours(cladeColours, alpha=0.5)
for(cph in sampledCPHs){
  
  for(key in cphs[[cph]]){

    if(is.null(groupOrHerdClusters[[cph]]) == FALSE){
      
      for(cluster in groupOrHerdClusters[[cph]]){
        colour <- cladeColoursRGB[as.numeric(cluster)]
        polygon(landParcelCoords[[key]], border=colour, lty=2, lwd=2)
      }
    }else{
      polygon(landParcelCoords[[key]], border=rgb(0,0,0, 0.5), lty=2, lwd=2)
    }
  }
}

#################################
# Plot the sampled cattle herds #
#################################
cladeColoursRGB <- getRGBsOfColours(cladeColours, alpha=1)
for(row in 1:nrow(cattleInfo)){
  
  if(is.null(isolateClades[[cattleInfo[row, "StrainId"]]]) == FALSE){
    
    points(x=cattleInfo[row, "Mapx"], y=cattleInfo[row, "Mapy"], pch=17,
           col=cladeColoursRGB[as.numeric(isolateClades[[cattleInfo[row, "StrainId"]]])],
           cex=1.5)
    
  }else{
    points(x=cattleInfo[row, "Mapx"], y=cattleInfo[row, "Mapy"], pch=17,
           col=rgb(0,0,0, 1), cex=1.5)
  }
}

#########################################
# Add Badger sampling locations to plot #
#########################################

# Plot the badger isolate sampling locations
cladeColoursRGB <- getRGBsOfColours(cladeColours, alpha=0.75)
addBadgerIsolatesLocationsWithClusterColours(badgerInfo=badgerInfo, 
                                             groupsRows=groupsRows,
                                             groupsCentroidsPerYear=groupsCentroidsPerYear, 
                                             defaultColour=rgb(0,0,0, 0.5),
                                             cex=2,
                                             cladeColoursRGB=cladeColoursRGB,
                                             isolateClades=isolateClades)

###############################################
# Add circles to define inner and outer demes #
###############################################

# Add deme defining circles
draw.circle(x=badgerCentre[1], y=badgerCentre[2], radius=thresholdDistance,
            border="black")
draw.circle(x=badgerCentre[1], y=badgerCentre[2], radius=outerDistance,
            border="black")
text(x=badgerCentre[1], y=(badgerCentre[2] - thresholdDistance) + 500,
     labels=paste("Inner: ", paste(thresholdDistance / 1000, "km", sep="")))
text(x=badgerCentre[1], y=(badgerCentre[2] - outerDistance) + 500,
     labels=paste("Outer: ", paste(outerDistance / 1000, "km", sep="")))

##############
# Add Legend #
##############

# Add legend
legend("bottomright", legend=c("Badger", "Cow"), text.col=c("black", "black"),
       bty="n", pch=c(20, 17), pt.cex=c(2, 1.5),
       col=c("black", "black"))

legend(x=387150, y=196000, legend=c("Cluster-0", "Cluster-1","Cluster-2","Cluster-3"),
       text.col=cladeColours, bty="n", cex=0.75)


dev.off()


#############
# FUNCTIONS #
#############

invertIsolatesInClades <- function(isolatesInClades){
  
  cladesOfIsolates <- list()
  for(clade in names(isolatesInClades)){
    
    for(isolate in isolatesInClades[[clade]]){
      cladesOfIsolates[[isolate]] <- clade
    }
  }
  
  return(cladesOfIsolates)
}

getSampledGroupsOrHerdsAssociatedWithClusters <- function(isolatesInClades, badgerInfo,
                                                          cattleInfo){
  groupOrHerdClusters <- list()
  for(clade in names(isolatesInClades)){
    
    for(isolate in isolatesInClades[[clade]]){
      
      # Get the sampled herd/group for the current isolate
      sampledGroupHerd <- getSampledHerdGroupOfIsolate(isolate, badgerInfo, cattleInfo)
      
      # Store the association
      if(is.null(groupOrHerdClusters[[sampledGroupHerd]]) == TRUE){
        groupOrHerdClusters[[sampledGroupHerd]] <- c(clade)
      }else{
        groupOrHerdClusters[[sampledGroupHerd]] <- c(groupOrHerdClusters[[sampledGroupHerd]],
                                                     clade)
      }
    }
  }
  
  # Make the cluster numbers unique for each herd/group
  for(key in names(groupOrHerdClusters)){
    groupOrHerdClusters[[key]] <- unique(groupOrHerdClusters[[key]])
  }
  
  return(groupOrHerdClusters)
}

getSampledHerdGroupOfIsolate <- function(isolate, badgerInfo, cattleInfo){
  
  # Check if badger or cow isolate
  if(grepl(x=isolate, pattern="WB") == TRUE){
    
    # Get row index
    row <- which(badgerInfo$WB_id == isolate)
    
    # Get sampled group
    sampledGroupHerd <- badgerInfo[row, "Social.Group.Trapped.At"]

    # Remove any spaces
    sampledGroupHerd <- paste(strsplit(sampledGroupHerd, split=" ")[[1]], collapse="")
    
  }else{
    
    # Get row index
    row <- which(cattleInfo$StrainId == isolate)
    
    # Get sampled group
    sampledGroupHerd <- cattleInfo[row, "CPH"]
  }
  
  return(sampledGroupHerd)
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

setAlpha <- function(colour, alpha){
  
  rgbValues <- col2rgb(colour)
  
  # Note that col2rgb() returns rgb values from 0 to 255
  rgbColour <- rgb(rgbValues["red", 1], rgbValues["green", 1], rgbValues["blue", 1],
                   alpha=alpha * 255, maxColorValue=255)
  
  return(rgbColour)
}

assignColoursToEachSampledCPH <- function(sampledCPHs){
  
  # Create a colour palette
  colourPalette <- colorRampPalette(c("red", "orange", "pink", "brown", "green", 
                                      "grey", "cyan", "purple", "violet", "blue"))
  colours = sample(colourPalette(length(sampledCPHs)))
  #colours <- palette(rainbow(nClusters))
  
  # Assign a colour to each CPH
  sampledCPHColours <- list()
  for(i in 1:length(sampledCPHs)){
    sampledCPHColours[[sampledCPHs[i]]] <- colours[i]
  }
  
  return(sampledCPHColours)
}

noteLandParcelsAssociatedWithCPHs <- function(polygonInfo){
  polygonInfo <- landParcels@data
  polygonInfo$PU_CPH <- as.character(polygonInfo$PU_CPH)
  cphs <- list()
  for(row in 1:nrow(polygonInfo)){
    
    if(is.null(cphs[[polygonInfo[row, "PU_CPH"]]]) == FALSE){
      
      cphs[[polygonInfo[row, "PU_CPH"]]] <- c(cphs[[polygonInfo[row, "PU_CPH"]]],
                                              as.character(row - 1))
      
    }else{
      cphs[[polygonInfo[row, "PU_CPH"]]] <- c(as.character(row - 1))
    }
  }
  
  return(cphs)
}

findBadgerGroupsCentre <- function(groupsCentroidsPerYear){
  
  Xs <- c()
  Ys <- c()
  index <- 0
  
  for(row in 1:nrow(groupsCentroidsPerYear)){
    
    for(col in 2:ncol(groupsCentroidsPerYear)){
      
      if(is.na(groupsCentroidsPerYear[row, col]) == FALSE){
        index <- index + 1
        
        values <- as.numeric(
          strsplit(groupsCentroidsPerYear[row, col], split=":")[[1]])
        Xs[index] <- values[1]
        Ys[index] <- values[2]
      }
    }
  }
  
  means <- c(mean(Xs), mean(Ys))
  
  return(means)
}

addBadgerIsolatesLocations <- function(badgerInfo, groupsRows, groupsCentroidsPerYear,
                                       col, cex){
  # Plot the isolate sampling locations (social group centroids)
  for(row in 1:nrow(badgerInfo)){
    
    year <- strsplit(badgerInfo[row, "date"], split="/")[[1]][3]
    group <- paste(strsplit(badgerInfo[row, "Social.Group.Trapped.At"], split=" ")[[1]],
                   collapse="")
    
    territoryCentroid <- strsplit(groupsCentroidsPerYear[groupsRows[[group]], 
                                                         as.character(year)],
                                  split=":")
    
    if(length(territoryCentroid) > 0){
      territoryCentroid <- as.numeric(territoryCentroid[[1]])
      points(x=territoryCentroid[1], y=territoryCentroid[2], pch=20,
             col=col, cex=cex)
    }
  }
}

addBadgerIsolatesLocationsWithClusterColours <- function(badgerInfo, groupsRows,
                                                         groupsCentroidsPerYear, 
                                                         defaultColour, cex,
                                                         cladeColoursRGB,
                                                         isolateClades){
  # Plot the isolate sampling locations (social group centroids)
  for(row in 1:nrow(badgerInfo)){
    
    year <- strsplit(badgerInfo[row, "date"], split="/")[[1]][3]
    group <- paste(strsplit(badgerInfo[row, "Social.Group.Trapped.At"], split=" ")[[1]],
                   collapse="")
    
    territoryCentroid <- strsplit(groupsCentroidsPerYear[groupsRows[[group]], 
                                                         as.character(year)],
                                  split=":")
    
    if(length(territoryCentroid) > 0){
      territoryCentroid <- as.numeric(territoryCentroid[[1]])
      
      if(is.null(isolateClades[[badgerInfo[row, "WB_id"]]]) == FALSE){
        
        points(x=territoryCentroid[1], y=territoryCentroid[2], pch=20,
               col=cladeColoursRGB[
                 as.numeric(isolateClades[[badgerInfo[row, "WB_id"]]])],
               cex=cex)
        
      }else{
        points(x=territoryCentroid[1], y=territoryCentroid[2], pch=20,
               col=defaultColour, cex=cex)
      }
    }
  }
}

getSequenceIDsFromFastaFile <- function(fileName){
  
  # Get the lines from the FASTA file
  connection <- file(fileName, open='r')
  fileLines <- readLines(connection)
  close(connection)
  
  # Initialise a vector to store the IDs
  ids <- c()
  index <- 0
  
  # Read the fasta file line by line
  for(line in fileLines[-1]){
    
    # Skip sequence lines
    if(grepl(x=line, pattern=">") == FALSE){
      next
    }
    
    index <- index + 1
    ids[index] <- strsplit(substr(x=line, start=2, stop=nchar(line)), split="_")[[1]][1]
  }
  
  return(ids)
}

updateMax <- function(territoryCoords, max){
  
  newValues <- territoryCoords[["max"]]
  
  if(newValues[1] > max[1]){
    max[1] <- newValues[1]
  }
  
  if(newValues[2] > max[2]){
    max[2] <- newValues[2]
  }
  
  return(max)
}

updateMin <- function(territoryCoords, min){
  
  newValues <- territoryCoords[["min"]]
  
  if(newValues[1] < min[1]){
    min[1] <- newValues[1]
  }
  
  if(newValues[2] < min[2]){
    min[2] <- newValues[2]
  }
  
  return(min)
}

plotTerritories <- function(territoriesForEachYear, centre, expand,
                            colour){
  
  par(mar=c(0,0,0,0))
  plot(x=NULL, y=NULL, yaxt="n", xaxt="n", bty="n", ylab="",
       xlim=c(centre[1] - expand, centre[1] + expand), 
       ylim=c(centre[2] - expand, centre[2] + expand), asp=1,
       xlab="")
  
  for(year in names(territoriesForEachYear)){
    
    territoryCoords <- territoriesForEachYear[[year]]
    for(id in names(territoryCoords)){
      
      polygon(territoryCoords[[id]], border=colour)
    }
  }
}

plotTerritoriesWithClusterColours <- function(territoriesForEachYear, centre, expand,
                            defaultColour, cladeColours, groupOrHerdClusters, lwd){
  
  par(mar=c(0,0,0,0))
  plot(x=NULL, y=NULL, yaxt="n", xaxt="n", bty="n", ylab="",
       xlim=c(centre[1] - expand, centre[1] + expand), 
       ylim=c(centre[2] - expand, centre[2] + expand), asp=1,
       xlab="")
  
  for(year in names(territoriesForEachYear)){
    
    territoryCoords <- territoriesForEachYear[[year]]
    for(id in names(territoryCoords)){
      
      colour <- defaultColour
      if(is.null(groupOrHerdClusters[[id]]) == FALSE){
        
        for(cluster in groupOrHerdClusters[[id]]){
          
          colour <- cladeColours[as.numeric(cluster)]
          polygon(territoryCoords[[id]], border=colour, lwd=lwd)
        }
      }else{
        polygon(territoryCoords[[id]], border=colour, lwd=lwd)
      }
    }
  }
}

addSocialGroupNames <- function(groupsCentroidsPerYear, groupsRows){
  
  for(group in names(groupsRows)){
    groupsCentroidsPerYear[groupsRows[[group]], "SocialGroup"] <- group
  }
  
  return(groupsCentroidsPerYear)
}

indexGroups <- function(groups, groupsRows, index){
  
  for(i in 1:length(groups)){
    
    if(is.null(groupsRows[[groups[i]]]) == TRUE){
      index <- index + 1
      groupsRows[[groups[i]]] <- index
    }
  }
  
  return(groupsRows)
}

addGroupCentroidsFromCurrentYear <- function(groupsCentroidsPerYear, groupsRows, territoryCentroids, year){
  
  # Update the territory centroids for the current year
  for(group in names(territoryCentroids)){
    
    coords <- territoryCentroids[[group]]
    
    groupsCentroidsPerYear[groupsRows[[group]], as.character(year)] <- paste(coords[1], ":", coords[2], sep="")
  }
  
  return(groupsCentroidsPerYear)
}

initialiseTableToStoreTerritoryCentroids <- function(years){
  groupsCentroidsPerYear <- data.frame(SocialGroup=c(NA), stringsAsFactors=FALSE)
  for(year in years){
    groupsCentroidsPerYear[, as.character(year)] <- NA
  }
  
  return(groupsCentroidsPerYear)
}

calculateTerritoryCentroids <- function(territoryCoords,
                                        territoryIDs){
  
  territoryCentroids <- list()
  polygonIDs <- names(territoryCoords)
  for(id in polygonIDs){
    
    if(id == "min" | id == "max"){
      next
    }
    
    territoryCentroids[[toupper(territoryIDs[[id]])]] <- c(
      mean(territoryCoords[[id]][, 1]),
      mean(territoryCoords[[id]][, 2]))
  }
  
  return(territoryCentroids)
}

assignTerritoriesToSocialGroupNames <- function(territoryIDs, territoryCoords){
  socialGroupTerritories <- list()
  for(id in names(territoryIDs)){
    socialGroupTerritories[[toupper(territoryIDs[[id]])]] <- territoryCoords[[id]]
  }
  
  return(socialGroupTerritories)
}

getSocialGroupNames <- function(polygonInfo, year){
  territoryIDs <- list()
  rowNames <- rownames(polygonInfo)
  
  column <- 1
  if(year == 2007){
    column <- 2
  }
  
  for(row in 1:nrow(polygonInfo)){
    socialGroupName <- as.character(polygonInfo[row, column])
    territoryIDs[[rowNames[row]]] <- removeSep(socialGroupName, sep=" ")
  }
  
  return(territoryIDs)
}

removeSep <- function(string, sep){
  parts <- strsplit(x=string, split=sep)[[1]]
  return(paste(parts, collapse=""))
}

getPolygonCoords <- function(spatialPolygonsDataFrame){
  polygonCoords <- list()
  
  polygonCoords[["min"]] <- c(99999999, 9999999)
  polygonCoords[["max"]] <- c(0, 0)
  
  for(i in 1:length(spatialPolygonsDataFrame@polygons)){
    polygonCoords[[spatialPolygonsDataFrame@polygons[[i]]@ID]] <- 
      spatialPolygonsDataFrame@polygons[[i]]@Polygons[[1]]@coords
    
    rangeX <- range(spatialPolygonsDataFrame@polygons[[i]]@Polygons[[1]]@coords[, 1])
    rangeY <- range(spatialPolygonsDataFrame@polygons[[i]]@Polygons[[1]]@coords[, 2])
    
    if(rangeX[1] < polygonCoords[["min"]][1]){
      polygonCoords[["min"]][1] <- rangeX[1]
    }
    if(rangeX[2] > polygonCoords[["max"]][1]){
      polygonCoords[["max"]][1] <- rangeX[2]
    }
    if(rangeY[1] < polygonCoords[["min"]][2]){
      polygonCoords[["min"]][2] <- rangeY[1]
    }
    if(rangeY[2] > polygonCoords[["max"]][2]){
      polygonCoords[["max"]][2] <- rangeY[2]
    }
  }
  
  return(polygonCoords)
}


