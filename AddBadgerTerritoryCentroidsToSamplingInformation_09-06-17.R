##################
# Load Libraries #
##################

library(maptools) # Read shape file
library(rgeos) # Polygon centroids

# Create a path variable
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/"

###################################################################
# Get the badger territory centroids each year they are available #
###################################################################

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
  file <- paste(path, "BadgerTerritoryMarkingData/",
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

##########################################
# Read badger sampling information table #
##########################################

# Read in the isolate metadata
fileName <- paste(path, "IsolateData/", "BadgerInfo_08-04-15_LatLongs_XY.csv",
                  sep="")
metadata <- read.table(fileName, header=TRUE, stringsAsFactors=FALSE, sep=",")

###################################
# Add badger centroids onto table #
###################################

# Add the badger territory centroids
metadata <- addBadgerIsolatesLocations(metadata, groupsRows, groupsCentroidsPerYear)

# Print the table to file
fileName <- paste(path, "IsolateData/", "BadgerInfo_08-04-15_LatLongs_XY_Centroids.csv",
                  sep="")
write.table(metadata, file=fileName, row.names=FALSE,
            quote=FALSE, sep=",")


#############
# FUNCTIONS #
#############

addBadgerIsolatesLocations <- function(metadata, groupsRows, groupsCentroidsPerYear){
  
  metadata$GroupCentroidX <- rep(NA, nrow(metadata))
  metadata$GroupCentroidY <- rep(NA, nrow(metadata))
  
  for(row in 1:nrow(metadata)){
    
    year <- strsplit(metadata[row, "date"], split="/")[[1]][3]
    group <- paste(strsplit(metadata[row, "Social.Group.Trapped.At"], split=" ")[[1]],
                   collapse="")
    
    # Fix JACKS MIREY - sometimes referred to as JACKS
    if(group == "JACKSMIREY"){
      group <- "JACKS"
    }else if(group == "TOPSETT"){
      group <- "TOP"
    }else if(group == "NA"){
      next
    }
    
    # Get the group's territory - if available
    if(is.null(groupsRows[[group]]) == FALSE){
      territoryCentroid <- strsplit(groupsCentroidsPerYear[groupsRows[[group]], 
                                                           year],
                                    split=":")[[1]]
    }else{
      print(group)
    }
    
    if(length(territoryCentroid) > 0){
      metadata[row, "GroupCentroidX"] <- as.numeric(territoryCentroid[1])
      metadata[row, "GroupCentroidY"] <- as.numeric(territoryCentroid[2])
    }
  }
  
  return(metadata)
}

addGroupCentroidsFromCurrentYear <- function(groupsCentroidsPerYear, groupsRows, territoryCentroids, year){
  
  # Update the territory centroids for the current year
  for(group in names(territoryCentroids)){
    
    coords <- territoryCentroids[[group]]
    
    groupsCentroidsPerYear[groupsRows[[group]], as.character(year)] <- paste(coords[1], ":", coords[2], sep="")
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

removeSep <- function(string, sep){
  parts <- strsplit(x=string, split=sep)[[1]]
  return(paste(parts, collapse=""))
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

initialiseTableToStoreTerritoryCentroids <- function(years){
  groupsCentroidsPerYear <- data.frame(SocialGroup=c(NA), stringsAsFactors=FALSE)
  for(year in years){
    groupsCentroidsPerYear[, as.character(year)] <- NA
  }
  
  return(groupsCentroidsPerYear)
}

addSocialGroupNames <- function(groupsCentroidsPerYear, groupsRows){
  
  for(group in names(groupsRows)){
    groupsCentroidsPerYear[groupsRows[[group]], "SocialGroup"] <- group
  }
  
  return(groupsCentroidsPerYear)
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
