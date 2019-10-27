# Load the rgdal library
library(rgdal)

# Read in the shape file
shapeFile <- file.path("~", "Desktop", "ROI_CountyBoundaries", "counties.shp")
spatialData <- readOGR(dsn = shapeFile)

# Extract the polygon coords
polygonCoords <- getPolygonCoords(spatialData)

# Plot the polygons
plotPolygons(polygonCoords)

# Add my current location
myLocation <- c(-6.7008791,55.1319208)
arrows(x0=myLocation[1]+0.1, x1=myLocation[1], y0=myLocation[2]+0.3, y1=myLocation[2],
       lwd=2, col="red", length=0.1)
text(x=myLocation[1]+0.1, y=myLocation[2]+0.4, labels="I am here!")

# Extract the polygon information
polygonInfo <- spatialData@data

#### FUNCTIONS ####

plotPolygons <- function(polygonCoords, ...){
  
  # Get and set the plotting margins
  currentMar <- par()$mar
  par(mar=c(0,0,0,0))
  
  # Create an empty plot
  plot(x=NULL, y=NULL, xlim=polygonCoords$BoundaryLimits[1, ], ylim= polygonCoords$BoundaryLimits[2, ],
       bty="n", xlab="", ylab="", xaxt="n", yaxt="n")
  
  # Get the names of the sets of polygons
  setNames <- names(polygonCoords)
  
  # Examine each set of polygons
  for(setID in setNames){
    
    # Ignore the Boundary Limits
    if(setID == "BoundaryLimits"){
      next
    }
    
    # Examine each of the polygons in the current set
    for(polygonIndex in seq_along(polygonCoords[[setID]])){
      
      # Get the coordinates for the current polygon
      coords <- polygonCoords[[setID]][[polygonIndex]]
      
      # Plot the current polygon
      polygon(coords, ...)
    }
  }
  
  # Reset the plotting margins
  par(mar=currentMar)
}

getPolygonCoords <- function(spatialDataFrame){
  
  # Note the number of polygon sets - one for each ID
  nSets <- length(spatialDataFrame@polygons)
  
  # Initialise a list to store the IDs and coordinates of each polygon
  polygonCoords <- list()
  
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
      
      # Store the coordinates
      polygons[[polygonIndex]] <- coords
    }
    
    # Store the polygon coordinates
    polygonCoords[[id]] <- polygons
  }
  
  # Store the boundary coordinates for the polygon
  polygonCoords$BoundaryLimits <- spatialDataFrame@bbox
  
  return(polygonCoords)
}