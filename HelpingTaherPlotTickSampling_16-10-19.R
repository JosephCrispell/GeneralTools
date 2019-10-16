#### Preparation ####

# Load libraries
library(rgdal) # For reading in shape files
library(sp)

# Set the path variable
path <- "/home/josephcrispell/Desktop/"

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

#### Read in the Ireland shape files ####

# Read in the Ireland shape file
file <- paste0(path, "Research/ROI_CountyBoundaries/counties.shp")
ireland <- readOGR(file) # Generates SpatialPolygonsDataFrame

# Get the polygons 
irelandPolygons <- getPolygonCoords(ireland)

#### Read in the tick count data ####

# Read in the tick counts
# Remove location columns as had characters interferring with reading in
file <- paste0(path, "HelpingTaher/Taher_AllSitesAndPathogens_16-10-19.tsv")
tickCounts <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

#### Plot the sampling ####

# Get and set the plotting margins
currentMar <- par()$mar
par(mar=c(0,0,0,0))

# Plot the Ireland county polygons
plotPolygons(irelandPolygons, lwd=1.5)

# Plot the tick counts
points(x=tickCounts$Long, y=tickCounts$Lat, pch=19, cex=2, col=rgb(0,0,0, 0.5))


# Reset the plotting margins
par(mar=currentMar)

#### FUNCTIONS ####

plotPolygons <- function(polygons, xLim=NULL, yLim=NULL, ...){
  
  # Get the axis limits of uk - enough to include ROI
  if(is.null(xLim)){
    limits <- getAxisLimits(polygons)
    xLim <- limits$X
    yLim <- limits$Y
  }
  
  # Create an empty plot
  plot(x=NULL, y=NULL, xlim=xLim, ylim=yLim,
       bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
  
  # Examine each set of polygons
  for(setName in names(polygons)){

    # Examine each polygon within current set
    for(polygonIndex in 1:length(polygons[[setName]])){
      
      # Plot the current polygon
      polygon(polygons[[setName]][[polygonIndex]], ...)
    }
  }
}

getAxisLimits <- function(polygons){
  
  # Initialise list to store the ranges of the X and Y axes
  ranges <- list("X"=c(Inf, -Inf), "Y"=c(Inf, -Inf))
  
  # Examine each of the polygons in the input coords
  for(polygonIndex in 1:length(polygons[[1]])){
    
    # Note the ranges of the X and Y axes of the current polygon
    xRange <- range(polygons[[1]][[polygonIndex]][, 1])
    yRange <- range(polygons[[1]][[polygonIndex]][, 2])
    
    # Update X min if necessary
    if(xRange[1] < ranges$X[1]){
      ranges$X[1] <- xRange[1]
    }
    
    # Update X max if necessary
    if(xRange[2] > ranges$X[2]){
      ranges$X[2] <- xRange[2]
    }
    
    # Update Y min if necessary
    if(yRange[1] < ranges$Y[1]){
      ranges$Y[1] <- yRange[1]
    }
    
    # Update Y max if necessary
    if(yRange[1] > ranges$Y[2]){
      ranges$Y[2] <- yRange[2]
    }
  }
  
  return(ranges)
}

getPolygonCoords <- function(spatialDataFrame){
  
  # Got some good information from: https://stackoverflow.com/questions/29803253/r-extracting-coordinates-from-spatialpolygonsdataframe

  # Note the number of polygon sets - one for each ID
  nSets <- length(spatialDataFrame@polygons)
  
  # Initialise a list to store the IDs and coordinates of each polygon
  output <- list()
  
  # Loop through all the sets of polygons
  for(setIndex in seq_len(nSets)){
    
    # Get the ID of the current polygon set
    id <- spatialDataFrame@polygons[[setIndex]]@ID
    
    # Note the number of polygons in the current set
    nPolygons <- length(spatialDataFrame@polygons[[setIndex]]@Polygons)
    
    # Initialise a list to store the coordinates of each polygon in current set
    polygons <- list()
    
    # Examine each polygon in current set
    for(polygonIndex in seq_len(nPolygons)){
      
      # Get the coordinates for the current polygon
      coords <- spatialDataFrame@polygons[[setIndex]]@Polygons[[polygonIndex]]@coords
      
      # Store the lX and Y coords
      polygons[[polygonIndex]] <- coords
    }
    
    # Store the polygon coordinates
    output[[id]] <- polygons
  }
  
  return(output)
}
