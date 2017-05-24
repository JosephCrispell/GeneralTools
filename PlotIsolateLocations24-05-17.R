##################
# Load Libraries #
##################

library(maptools) # Read shape file
library(rgeos) # Polygon centroids

# Create a path variable
path <- "/Users/josephcrisp1/Desktop/BadgerTerritoryShapeFiles/"

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
groupsCentroidsPerYear <- data.frame(SocialGroup=c(NA), stringsAsFactors=FALSE)
for(year in years){
  groupsCentroidsPerYear[, as.character(year)] <- NA
}
groupRows <- list()

# Examine the territories in each year
for(i in 1:length(years)){
  
  year <- years[i]
  shapeFileName <- shapeFileNames[i]
  
  # Read in the shape file
  file <- paste(path, "Baitmarking ", year, "/", shapeFileName, sep="")
  territories <- readShapePoly(file) # Generates SpatialPolygonsDataFrame
  
  # Extract the polygon coordinates
  territoryCoords <- getPolygonCoords(territories)
  
  # Get the full social group names
  territoryIDs <- getSocialGroupNames(territories@data, year)
  
  # Assign polygons to their social group names
  socialGroupTerritories <- assignTerritoriesToSocialGroupNames(territoryIDs, territoryCoords)
  
  # Calculate the territory centroids - mean X and Y
  territoryCentroids <- calculateTerritoryCentroids(territoryCoords, territoryIDs)
  
  # Note the centroids of each group in the current year
}


#############
# FUNCTIONS #
#############

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


