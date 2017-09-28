##################
# Load Libraries #
##################

library(maptools) # Read shape file
library(rgeos) # Polygon centroids

###################################
# Read in ROI counties shape file # https://www.townlands.ie/page/download/
###################################

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Paratuberculosis/"

# Read in the shape file
file <- paste(path, "CountyPolygonCoords/ROI_CountyBoundaries/counties.shp", sep="")
countyBorders <- readShapePoly(file) # Generates SpatialPolygonsDataFrame

# Get the coordinates of the counties
countyCoords <- getPolygonCoords(countyBorders)

# Get the county names associated with the coords
countyNames <- getPolygonNames(countyBorders@data, "NAME_TAG")

#########################
# Read in MAP VNTR data #
#########################

# Read in the file
file <- paste(path, "Genotyping data.csv", sep="")
vntrInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE,
                       check.names=FALSE)

# Remove rows with no DNA extraction
vntrInfo <- vntrInfo[vntrInfo[, "Date DNA extraction"] != "", ]

# Count number of samples per county
nSamplesPerCounty <- countNumberSamplesPerCounty(vntrInfo)
maxCount <- max(getValues(nSamplesPerCounty))

#####################
# Plot the counties #
#####################

plotNSamplesPerCounty(countyCoords, countyNames, nSamplesPerCounty)

#################################
# Print out the county polygons #
#################################

writePolygonCoordsToFile(countyCoords, countyNames, path)





#############
# FUNCTIONS #
#############

plotNSamplesPerCounty <- function(countyCoords, countyNames, nSamplesPerCounty){
  par(mar=c(0,0,0,0))
  
  plot(x=NULL, y=NULL, yaxt="n", xaxt="n", ylab="", xlab="", bty="n",
       xlim=c(countyCoords[["min"]][1], countyCoords[["max"]][1]),
       ylim=c(countyCoords[["min"]][2], countyCoords[["max"]][2]))
  
  for(key in names(countyCoords)){
    
    if(key %in% c("min", "max")){
      next
    }
    
    nSamples <- nSamplesPerCounty[[countyNames[[key]]]]
    if(length(nSamples) == 0){
      nSamples <- 0
    }
    
    polygon(countyCoords[[key]], border=rgb(0,0,0, 1), 
            col=rgb(0,0,1,nSamples / maxCount),
            lwd=2)
    
    if(nSamples != 0){
      text(x=mean(countyCoords[[key]][, 1]),
           y=mean(countyCoords[[key]][, 2]),
           labels=paste(countyNames[[key]], " (", nSamples, ")", sep=""),
           cex=0.6, col="red")
    }else{
      text(x=mean(countyCoords[[key]][, 1]),
           y=mean(countyCoords[[key]][, 2]),
           labels=countyNames[[key]], cex=0.6, col="gray50")
    }
  }
}

writePolygonCoordsToFile <- function(countyCoords, countyNames, path){
  for(key in names(countyCoords)){
    
    if(key %in% c("min", "max")){
      next
    }
    
    file <- paste(path, "PolygonCoords_", countyNames[[key]], ".txt", sep="")
    table <- countyCoords[[key]]
    colnames(table) <- c("X", "Y")
    
    write.table(table, file, sep="\t",
                row.names=FALSE, quote=FALSE)
  }
}

getValues <- function(list){
  
  values <- c()
  for(key in names(list)){
    values[length(values) + 1] <- list[[key]]
  }
  
  return(values)
}

invertList <- function(list){
  
  output <- list()
  for(key in names(list)){
    output[[as.character(list[[key]])]] <- key
  }
  
  return(output)
}

countNumberSamplesPerCounty <- function(vntrInfo){
  counties <- list()
  for(row in 1:nrow(vntrInfo)){
    
    county <- strsplit(vntrInfo[row, "Herd Location"], split=" ")[[1]][1]
    if(is.null(counties[[county]]) == TRUE){
      counties[[county]] <- 1
    }else{
      counties[[county]] <- counties[[county]] + 1
    }
  }
  
  return(counties)
}

getPolygonNames <- function(polygonInfo, column){
  names <- list()
  rowNames <- rownames(polygonInfo)
  
  for(row in 1:nrow(polygonInfo)){
    names[[rowNames[row]]] <- as.character(polygonInfo[row, column])
  }
  
  return(names)
}

removeSep <- function(string, sep){
  parts <- strsplit(x=string, split=sep)[[1]]
  return(paste(parts, collapse=""))
}

getPolygonCoords <- function(spatialPolygonsDataFrame){
  polygonCoords <- list()
  
  polygonCoords[["min"]] <- c(NA, NA)
  polygonCoords[["max"]] <- c(NA, NA)
  
  for(i in 1:length(spatialPolygonsDataFrame@polygons)){
    
    
    polygonCoords[[spatialPolygonsDataFrame@polygons[[i]]@ID]] <- 
      spatialPolygonsDataFrame@polygons[[i]]@Polygons[[1]]@coords
    
    rangeX <- range(spatialPolygonsDataFrame@polygons[[i]]@Polygons[[1]]@coords[, 1])
    rangeY <- range(spatialPolygonsDataFrame@polygons[[i]]@Polygons[[1]]@coords[, 2])
    
    if(is.na(polygonCoords[["min"]][1]) == TRUE ||
       rangeX[1] < polygonCoords[["min"]][1]){
      polygonCoords[["min"]][1] <- rangeX[1]
    }
    if(is.na(polygonCoords[["max"]][1]) == TRUE ||
       rangeX[2] > polygonCoords[["max"]][1]){
      polygonCoords[["max"]][1] <- rangeX[2]
    }
    if(is.na(polygonCoords[["min"]][2]) == TRUE ||
       rangeY[1] < polygonCoords[["min"]][2]){
      polygonCoords[["min"]][2] <- rangeY[1]
    }
    if(is.na(polygonCoords[["max"]][2]) == TRUE ||
       rangeY[2] > polygonCoords[["max"]][2]){
      polygonCoords[["max"]][2] <- rangeY[2]
    }
  }
  
  return(polygonCoords)
}
