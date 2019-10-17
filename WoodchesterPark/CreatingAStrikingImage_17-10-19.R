#### Preparation ####

# Load libraries
library(rgdal) # For reading in shape files and converting between projections
library(basicPlotteR) # For setting alpha
library(OpenStreetMap) # Great tutorial here: https://www.r-bloggers.com/the-openstreetmap-package-opens-up/


# Set the path variable
path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/"

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

#### Get badger territories through time ####

# Note the years that territory shape files are available for
years <- c(2000:2011)

# Note the name of the shape for the territories from each year - varies in structure
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

# Read in the polygon coordinates for each social gropu in each year
territoryCoordsInEachYear <- readTerritoryCoordsFromEachYear(years, shapeFileNames)

#### Create a striking image ####

# Plot the badger territories on top of one another
plotBadgerTerritories(territoryCoordsInEachYear, c(2000:2005,2007:2011), scale=FALSE,
                      lwd=3, border=rgb(0,0,0, 0.95))

# Overlay the badger movements

#### FUNCTIONS ####

plotBadgerTerritories <- function(territoryCoordsInEachYear, years, sleep=NULL, file=NULL, map=NULL,
                                  scale=TRUE, ...){
  
  # Open a pdf file if requested
  if(is.null(file) == FALSE){
    pdf(file)
  }
  
  # Get and set the plotting margins
  currentMar <- par("mar")
  par(mar=c(0,0,0,0))
  
  # Check if map provided for plotting
  if(is.null(map) == FALSE){
    plot(map, asp=1)
  }else{
    
    # Create an empty plot
    plot(x=NULL, y=NULL, xlim=territoryCoordsInEachYear$rangeX, ylim=territoryCoordsInEachYear$rangeY,
         bty="n", xaxt="n", yaxt="n", xlab="", ylab="", asp=1)
  }
  
  # Examine each year
  for(year in years){
    
    # Skip 2006 - territories are all messed up
    if(year == 2006){
      next
    }
    
    # Convert year to character to act as key
    year <- as.character(year)
    
    # Pause - if asked
    if(is.null(sleep) == FALSE){
      Sys.sleep(sleep)
    }
    
    # Examine each social group
    for(socialGroup in names(territoryCoordsInEachYear[[year]])){
      
      # Check only one polygon available for current territory
      if(length(territoryCoordsInEachYear[[year]][[socialGroup]]) > 1){
        stop(paste0("More than one polygon available for ", socialGroup))
      }
      
      # Get the coordinates for current social group's polygon
      coords <- territoryCoordsInEachYear[[year]][[socialGroup]][[1]]
      
      # Plot the coordinates
      polygon(coords, ...)
    }
  }
  
  # Add scale
  if(scale){
    axisLimits <- par()$usr
    xLength <- axisLimits[2] - axisLimits[1]
    yLength <- axisLimits[4] - axisLimits[3]
    xPad <- 0.08*xLength
    points(x=c(axisLimits[2] - xPad - 1000, axisLimits[2] - xPad), y=c(axisLimits[3] + 0.1*yLength, axisLimits[3] + 0.1*yLength),
           type="l", lwd=4)
    text(x=axisLimits[2] - xPad - 500, y=axisLimits[3] + 0.07*yLength, labels="1KM", cex=2)
  }
  
  # Reset the plotting margins
  par(mar=currentMar)
  
  # Close the pdf if requested
  if(is.null(file) == FALSE){
    dev.off()
  }
}

readTerritoryCoordsFromEachYear <- function(years, shapeFileNames){
  
  # Initialise list to store the ranges of the X and Y axes
  ranges <- list("X"=c(Inf, -Inf), "Y"=c(Inf, -Inf))
  
  # Initialise a list to store the territory coordinates for each social group in each year
  territoriesForEachYear <- list()
  
  # Examine the territories in each year
  for(i in 1:length(years)){
    
    year <- years[i]
    shapeFileName <- shapeFileNames[i]
    
    # Read in the shape file
    file <- paste(path, "BadgerCaptureData/BadgerTerritoryMarkingData/",
                  "Baitmarking ", year, "/", shapeFileName, sep="")
    territories <- readOGR(file) # Generates SpatialPolygonsDataFrame
    
    # Get the polygon coordinates for each territory
    coords <- getPolygonCoords(territories)
    coords <- changeIDsToSocialGroupNames(coords, territories@data, year)
    
    # Update the X and Y ranges
    ranges <- updateXAndYRanges(ranges, coords)
    
    # Store the polygon coords under the current year
    territoriesForEachYear[[as.character(year)]] <- coords
  }
  
  # Store the range of X and Y axes
  territoriesForEachYear[["rangeX"]] <- ranges$X
  territoriesForEachYear[["rangeY"]] <- ranges$Y
  
  return(territoriesForEachYear)
}

updateXAndYRanges <- function(ranges, coords){
  
  # Examine each of the polygons in the input coords
  for(key in names(coords)){
    for(polygonIndex in seq_len(length(coords[[key]]))){
      
      # Note the ranges of the X and Y axes of the current polygon
      xRange <- range(coords[[key]][[polygonIndex]][, 1])
      yRange <- range(coords[[key]][[polygonIndex]][, 2])
      
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

changeIDsToSocialGroupNames <- function(coords, data, year){
  
  # Note the column that the social group names are kept in
  column <- 1
  if(year == 2007){
    column <- 2
  }
  
  # Initialise a list to store the coords based on social group name rather than ID
  output <- list()
  
  # Examine each polygon in coords
  for(key in names(coords)){
    
    # Convert the ID to a row
    row <- as.numeric(key) + 1
    
    # Convert the ID to a social group name
    socialGroup <- removeSep(as.character(data[row, column]), sep=" ")
    
    # Store the polygon coords under social group name
    output[[socialGroup]] <- coords[[key]]
  }
  
  return(output)
}

removeSep <- function(strings, sep){
  
  # Initialise a vector to store output
  output <- c()
  
  # Examine each of the input strings
  for(i in seq_along(strings)){
    
    # Split the current string into its parts
    parts <- strsplit(x=strings[i], split=sep)[[1]]
    
    # Store the string without sep
    output[i] <- paste(parts, collapse="")
  }
  
  return(output)
}
