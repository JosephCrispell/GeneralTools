#### Preparation ####

# Shape files downloaded from: https://gadm.org/download_country_v3.html

# Load libraries
library(rgdal) # For reading in shape files
library(sp)

# Set the path variable
path <- "/home/josephcrispell/Desktop/Research/"

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

#### Read in the UK shape file ####

# Read in the shape file
file <- paste(path, "UK_ShapeFile/gadm36_GBR_0.shp", sep="")
uk <- readOGR(file) # Generates SpatialPolygonsDataFrame

# Get the polygons 
ukPolygons <- getPolygonCoords(uk)

#### Read in the ROI shape file

# Read in the shape file
file <- paste(path, "ROI_ShapeFile/gadm36_IRL_0.shp", sep="")
roi <- readOGR(file) # Generates SpatialPolygonsDataFrame

# Get the polygons 
roiPolygons <- getPolygonCoords(roi)

#### Plot the polygons ####

# Open a pdf
file <- paste0(path, "Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/ESM_Figures/Collaborations_", date, ".pdf")
pdf(file)

# Get and set the plotting margins
currentMar <- par()$mar
par(mar=c(0,0,0,0))

# Plot the UK and ROI outlines
plotPolygons(ukPolygons, roiPolygons, xLim=c(-10.5,1.6), yLim=c(50,61), lwd=1.5)

# Plot the locations of interest
# UCD, Glasgow, Roslin, AFBNI, Woodchester & Weybridge
xCoords <- c(-6.2187348, -4.3278438, -3.2020387, -5.8265013, -2.277601, -0.4957386)
yCoords <- c(53.3052183, 55.9096218, 55.8508271, 54.600598, 51.710962, 51.3535598)
for(i in seq_along(xCoords)){
  for(j in seq_along(xCoords)){
    
    # Skip lower half
    if(i > j){
      next
    }
    
    points(x=c(xCoords[i], xCoords[j]), y=c(yCoords[i], yCoords[j]), type="l", lwd=3, col="grey")
  }
}
points(x=xCoords, y=yCoords, pch=19, col="blue", cex=2)

# Plot the UK and ROI outlines
plotPolygons(ukPolygons, roiPolygons, xLim=c(-10.5,1.6), yLim=c(50,61), lwd=1.5)

# Plot the locations of interest
# UCD, Glasgow, Roslin, AFBNI, Woodchester & Weybridge
points(x=-2.277601, y=51.710962, pch=19, col="blue", cex=1)

# Reset the plotting margins
par(mar=currentMar)

# Close the pdf
dev.off()

#### FUNCTIONS ####

plotPolygons <- function(uk, roi, xLim=NULL, yLim=NULL, ...){
  
  # Get the axis limits of uk - enough to include ROI
  if(is.null(xLim)){
    limits <- getAxisLimits(uk)
    xLim <- limits$X
    yLim <- limits$Y
  }
  
  # Create an empty plot
  plot(x=NULL, y=NULL, xlim=xLim, ylim=yLim,
       bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
  
  # Plot each of the UK polygons
  for(polygonIndex in 1:length(uk[[1]])){
      
    # Plot the coordinates
    polygon(uk[[1]][[polygonIndex]], ...)
  }
  
  # Plot each of the ROI polygons
  for(polygonIndex in 1:length(roi[[1]])){
    
    # Plot the coordinates
    polygon(roi[[1]][[polygonIndex]], ...)
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
  # Also from my blog: https://josephcrispell.github.io/BlogPosts/GetPolygonsFromShapeFile_11-10-17/GetPolygonsFromShapeFile_11-10-17.html
  # WHICH NEEDS UPDATED!!!!
  
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
