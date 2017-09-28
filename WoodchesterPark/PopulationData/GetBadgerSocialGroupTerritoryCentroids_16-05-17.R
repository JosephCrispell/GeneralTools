
library(maptools) # Read shape file
library(rgeos) # Polygon centroids

# Create a path variable
path <- "/Users/josephcrisp1/Desktop/BadgerTerritoryShapeFiles/"

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
groupsCentroidsPerYear <- list()

# Open a PDF file
file <- paste(path, "SocialGroupTerritoriesPerYear_16-05-17.pdf", sep="")
pdf(file)

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
  territoryIDs <- getSocialGroupNames(territories@data)
  
  # Calculate the territory centroids - mean X and Y
  territoryCentroids <- calculateTerritoryCentroids(
    territoryCoords, territoryIDs)
  
  # Record the centroids of each group in the current year
  groupsCentroidsPerYear <- recordCentroidsOfGroupInGivenYear(groupsCentroidsPerYear=groupsCentroidsPerYear,
                                                             year=year, years=years, 
                                                             territoryCentroids=territoryCentroids)
}

# Print the Centroid of each social group for each year available
outputTable <- createOutputTable(groupCentroidsPerYear, years)
file <- paste(path, "SocialGroupsCentroidsPerYear_16-05-17.txt", sep="")
write.table(x=outputTable, file=file, sep="\t", row.names=FALSE, quote=FALSE)

dev.off()


############
# FUNCTION #
############

createOutputTable <- function(groupCentroidsPerYear, years){
  
  socialGroups <- names(groupsCentroidsPerYear)
  outputTable <- data.frame(SocialGroup=socialGroups, stringsAsFactors=FALSE)
  for(year in years){
    outputTable[, as.character(year)] <- rep(NA, length(socialGroups))
  }
  
  for(i in 1:length(socialGroups)){
    
    outputTable[i, "SocialGroup"] <- socialGroups[i]
    groupCentroids <- groupsCentroidsPerYear[[socialGroups[i]]]
    
    for(row in 1:nrow(groupCentroids)){
      outputTable[i, row + 1] <- paste(groupCentroids[row, "X"], groupCentroids[row, "Y"], sep=":")
    }  
  }
  
  return(outputTable)
}

recordCentroidsOfGroupInGivenYear <- function(groupsCentroidsPerYear, year, years, territoryCentroids){
  for(group in names(territoryCentroids)){
    
    if(is.null(groupsCentroidsPerYear[[group]]) == FALSE){
      groupsCentroidsPerYear[[group]][(year-years[1])+1, "X"] <- territoryCentroids[[group]][1]
      groupsCentroidsPerYear[[group]][(year-years[1])+1, "Y"] <- territoryCentroids[[group]][2]
    }else{
      groupsCentroidsPerYear[[group]] <- data.frame(Years=years, X=rep(NA, length(years)),
                                                   Y=rep(NA, length(years)), stringsAsFactors=FALSE)
      groupsCentroidsPerYear[[group]][(year-years[1])+1, "X"] <- territoryCentroids[[group]][1]
      groupsCentroidsPerYear[[group]][(year-years[1])+1, "Y"] <- territoryCentroids[[group]][2]
    }
  }
  return(groupsCentroidsPerYear)
}

calculateTerritoryCentroids <- function(territoryCoords,
                                        territoryIDs){
  par(mar=c(0,0,0,0))
  plot(x=NULL, y=NULL, yaxt="n", xaxt="n", bty="n", ylab="",
       xlim=c(territoryCoords[["min"]][1], territoryCoords[["max"]][1]), 
       ylim=c(territoryCoords[["min"]][2], territoryCoords[["max"]][2]), asp=1,
       xlab="")
  legend("top", legend=year, bty="n", cex=2)
  legend("bottom", bty="n",
         legend=
           paste(
             round((territoryCoords[["max"]][1] - territoryCoords[["min"]][1])/1000, digits=1),
             "KM"))
  
  territoryCentroids <- list()
  polygonIDs <- names(territoryCoords)
  for(id in polygonIDs){
    
    if(id == "min" | id == "max"){
      next
    }
    
    territoryCentroids[[territoryIDs[[id]]]] <- c(
      mean(territoryCoords[[id]][, 1]),
      mean(territoryCoords[[id]][, 2]))
    
    points(territoryCoords[[id]], type="l")
    points(x=territoryCentroids[[territoryIDs[[id]]]][1], 
           y=territoryCentroids[[territoryIDs[[id]]]][2], type="p", pch=20)
  }
  
  return(territoryCentroids)
}

getSocialGroupNames <- function(polygonInfo){
  territoryIDs <- list()
  rowNames <- rownames(polygonInfo)
  for(row in 1:nrow(polygonInfo)){
    socialGroupName <- as.character(polygonInfo[row, 1])
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

