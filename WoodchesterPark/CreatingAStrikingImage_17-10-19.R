#### Preparation ####

# Load libraries
library(rgdal) # For reading in shape files and converting between projections
library(basicPlotteR) # For setting alpha
library(OpenStreetMap) # Downloading map
library(sp) # Convert between map projections

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

# Calculate the centroid of each social group in each year
territoryCentroidsInEachYear <- calculateTerritoryCentroidsInEachYear(territoryCoordsInEachYear)

#### Construct badger social group contact network for each year ####

# Read in the consolidated badger data
consolidatedCaptureDataFile <- paste0(path, "BadgerCaptureData/WP_CaptureData_Consolidated_24-06-2019.txt")
captureData <- read.table(consolidatedCaptureDataFile, header=TRUE, check.names=FALSE, sep="\t", 
                          stringsAsFactors=FALSE)

# Build a social group contact network for each year of interest
dispersalNetworksInEachYear <- countDispersalEventsInEachYear(captureData, territoryCoordsInEachYear)

#### Get a WP satellite image ####

# Calculate the plotting area boundaries in latitude and longitudes
coordRanges <- convertXYToLatLongs(x=territoryCoordsInEachYear$rangeX, y=territoryCoordsInEachYear$rangeY)
rangeLat <- coordRanges[, "Latitude"]
rangeLong <- coordRanges[, "Longitude"]

# Get a satellite image of the area
# Maps are initially put in a sperical mercator projection which is the standard for most (all?) map tiling systems
map <- openmap(upperLeft=c(rangeLat[2]+0.0012, rangeLong[2]), 
               lowerRight=c(rangeLat[1]+0.001, rangeLong[1]),
               type="bing")

# Convert from
mapUKGrid <- openproj(map, projection=CRS("+init=epsg:27700"))

#### Create a striking image ####

# Set the plotting margins
par(mar=c(0,0,0,0))

# Note the years we are plotting data for
years <- c(2003)

# Plot the badger territories on top of one another
plotBadgerTerritories(territoryCoordsInEachYear, years, scale=FALSE,
                      lwd=3, border=rgb(249,239,63, maxColorValue=255), 
                      #background=rgb(114,114,114, maxColorValue=255),
                      col=rgb(113,137,211, 50, maxColorValue=255), map=mapUKGrid)

rgb(113,137,211, maxColorValue=255) # cyan
rgb(183,48,121, maxColorValue=255) # pink
rgb(195,59,45, maxColorValue=255) # orange
rgb(50,63,98, maxColorValue=255) # dark blue
rgb(249,239,63, maxColorValue=255) # yellow

# Plot the territory centroids
#addTerritoryCentroids(territoryCentroidsInEachYear, years, pch=19, cex=1.5, col="red")

# Overlay the badger movements
addDispersalEvents(dispersalNetworksInEachYear, territoryCentroidsInEachYear, years, 
                   col=rgb(195,59,45, maxColorValue=255),
                   lty=1, lwd=6, weighted=FALSE)

#### FUNCTIONS ####

convertXYToLatLongs <- function(x, y, ukGrid="+init=epsg:27700"){
  
  # Create variables for holding the coordinate system types
  # see http://www.epsg.org/
  latLong <- "+init=epsg:4326"
  
  # Create a coordinates variable
  coords <- cbind(Easting = x, Northing = y)
  
  # Create a SpatialPointsDataFrame
  spatialDF <- SpatialPoints(coords, proj4string = CRS(ukGrid))
  
  # Convert the Eastings and Northings to Latitude and Longitude
  spatialDFLatLongs <- spTransform(spatialDF, CRS(latLong))
  
  # Create a table to store the converted points
  output <- data.frame(Longitude=spatialDFLatLongs@coords[,"Easting"],
                       Latitude=spatialDFLatLongs@coords[,"Northing"])
  
  # Return the lat longs
  return(output)
}

addDispersalEvents <- function(dispersalNetworksInEachYear, territoryCentroidsInEachYear, years, weighted=TRUE,
                               type="l", lwd=1, lend=2, ...){
  
  # Set the line ending style
  par(lend=lend)
  
  # Examine each year
  for(year in years){
    
    # Skip 2006 - territories are all messed up
    if(year == 2006){
      next
    }
    
    # Convert year to character to act as key
    year <- as.character(year)
    
    # Get the dispersal network for the current year
    dispersalNetwork <- dispersalNetworksInEachYear[[year]]
    
    # Get the territory centroids for the current year
    territoryCentroids <- territoryCentroidsInEachYear[[year]]
  
    # Note the maxium number of dispersal events
    maxNoDispersalEvents <- max(dispersalNetwork)
    
    # Examine the dispersal event network
    for(row in 1:nrow(dispersalNetwork)){
      
      # Note the centroid of the current row social group
      centroidRow <- territoryCentroids[[rownames(dispersalNetwork)[row]]]
      
      for(column in 1:ncol(dispersalNetwork)){
        
        # Skip if no dispersal events between current group
        if(dispersalNetwork[row, column] == 0){
          next
        }
        
        # Note the centroid of the current column social group
        centroidColumn <- territoryCentroids[[rownames(dispersalNetwork)[column]]]
        
        # Add a line to represent the dispersal event - check if want weighted by number of
        if(weighted){
          points(x=c(centroidRow[1], centroidColumn[1]),
                 y=c(centroidRow[2], centroidColumn[2]),
                 lwd=(dispersalNetwork[row, column] / maxNoDispersalEvents) * lwd,
                 type=type, ...)
        }else{
          
          # Check if haven't already connected these social groups
          if(row < column && dispersalNetwork[column, row] != 0){
            next
          }
          
          points(x=c(centroidRow[1], centroidColumn[1]),
                 y=c(centroidRow[2], centroidColumn[2]),
                 type=type, lwd=lwd, ...)
        }
      }
    }
  }
}

initialiseSocialGroupsDispersalNetworkForEachYear <- function(territoryCoordsInEachYear){
  
  # Initialise a list to store the social group dispersal networks for each year
  socialGroupsDispersalNetworkInEachYear <- list()
  
  # Each the data for each year
  for(year in names(territoryCoordsInEachYear)){
    
    # Skip the X and Y ranges
    if(year %in% c("rangeX", "rangeY")){
      next
    }
    
    # Note the social groups present in the current year
    socialGroups <- names(territoryCoordsInEachYear[[year]])
    
    # Initialise a matrix to store the number of dispersal events in the current year
    dispersalEventsCounts <- matrix(0, nrow=length(socialGroups), ncol=length(socialGroups))
    rownames(dispersalEventsCounts) <- socialGroups
    colnames(dispersalEventsCounts) <- socialGroups
    
    # Store the matrix
    socialGroupsDispersalNetworkInEachYear[[year]] <- dispersalEventsCounts
  }
  
  return(socialGroupsDispersalNetworkInEachYear)
}

countDispersalEventsInEachYear <- function(captureData, territoryCoordsInEachYear){
  
  # Initialise a list to store the social group dispersal networks for each year
  socialGroupsDispersalNetworkInEachYear <- 
    initialiseSocialGroupsDispersalNetworkForEachYear(territoryCoordsInEachYear)
  
  # Examine the capture history for each badger
  for(row in 1:nrow(captureData)){
    
    # Get the capture dates
    captureDates <- as.Date(strsplit(captureData[row, "@CaptureDates"], split=";")[[1]], format="%d-%m-%Y")
    
    # Get the social groups at capture
    socialGroups <- strsplit(captureData[row, "@Groups"], split=";")[[1]]
    
    # Skip individual if onlly one capture event
    if(length(captureDates) == 1){
      next
    }
    
    # Record info for any dispersal events for the current badger
    dispersalEvents <- noteDispersalEvents(captureDates, socialGroups, captureData[row, "Tattoo"])
    
    # Skip iof no dispersal events found
    if(length(dispersalEvents) == 0){
      next
    }
    
    # Record the dispersal events in the dispersal networks of each year
    socialGroupsDispersalNetworkInEachYear <- recordDispersalEvents(dispersalEvents, 
                                                                    socialGroupsDispersalNetworkInEachYear)
    
  }
  
  return(socialGroupsDispersalNetworkInEachYear)
}

recordDispersalEvents <- function(dispersalEvents, socialGroupsDispersalNetworkInEachYear){
  
  # Examine each of the dispersal events
  for(eventIndex in seq_along(dispersalEvents)){
    
    # Get the information for the current event
    year <- dispersalEvents[[eventIndex]]$Year
    from <- dispersalEvents[[eventIndex]]$FROM
    to <- dispersalEvents[[eventIndex]]$TO
    
    # Skip if not a year we have territory data for
    if(year %in% names(socialGroupsDispersalNetworkInEachYear) == FALSE){
      next
    }
    
    # Get the dispersal network for the current year
    network <- socialGroupsDispersalNetworkInEachYear[[year]]
    
    # Skip if FROM or TO social groups don't exist within the dispersal network
    if(from %in% rownames(network) == FALSE || to %in% rownames(network) == FALSE){
      next
    }
    
    # Add to the specific count for the current dispersal event
    network[from, to] <- network[from, to] + 1
    socialGroupsDispersalNetworkInEachYear[[year]] <- network
  }
  
  return(socialGroupsDispersalNetworkInEachYear)
}

noteDispersalEvents <- function(captureDates, socialGroups, tattoo){
  
  # Initialise a list to store any dispersal details
  dispersalEvents <- list()
  
  # Examine each capture event
  for(eventIndex in seq_along(captureDates)[-1]){
    
    # Skip if previous or current gropu is unknown
    if(socialGroups[eventIndex-1] == "NA" || socialGroups[eventIndex] == "NA"){
      next
    }
    
    # Skip if previous and current group are the same
    if(socialGroups[eventIndex-1] == socialGroups[eventIndex]){
      next
    }
    
    # Record the information for the current dispersal event
    dispersalEvents[[length(dispersalEvents) + 1]] <- list("Tattoo"=tattoo,
                                                           "FROM"=socialGroups[eventIndex-1],
                                                           "TO"=socialGroups[eventIndex],
                                                           "Year"=format(captureDates[eventIndex], "%Y"),
                                                           "Date"=captureDates[eventIndex])
  }
  
  return(dispersalEvents)
}

getAllSocialGroups <- function(territoryCoordsInEachYear){
  
  # Initialise a vector to store the social group name
  socialGroups <- c()
  
  # Add the social groups from each year into the vector
  for(year in names(territoryCoordsInEachYear)){
    
    socialGroups <- c(socialGroups, names(territoryCoordsInEachYear[[year]]))
  }
  
  # Return the unique names
  return(unique(socialGroups))
}

addTerritoryCentroids <- function(territoryCentroidsInEachYear, years, ...){
  
  # Examine each year
  for(year in years){
    
    # Skip 2006 - territories are all messed up
    if(year == 2006){
      next
    }
    
    # Convert year to character to act as key
    year <- as.character(year)
    
    # Examine each social group
    for(socialGroup in names(territoryCentroidsInEachYear[[year]])){

      # Get the centroid for the current territory in the current year
      centroidCoords <- territoryCentroidsInEachYear[[year]][[socialGroup]]

      # Plot a point at current centroid
      points(x=centroidCoords[1], y=centroidCoords[2], ...)
    }
  }
}

calculateTerritoryCentroidsInEachYear <- function(territoryCoordsInEachYear){
  
  # Initialise a list to stroe the centroid of each social group in each year
  territoryCentroidsInEachYear <- list()
  
  # Examine each of the years
  for(year in names(territoryCoordsInEachYear)){
    
    # Examine each social group present in the current year
    for(socialGroup in names(territoryCoordsInEachYear[[year]])){
      
      # Get the coordinates for the territory of the current social group
      coords <- territoryCoordsInEachYear[[year]][[socialGroup]][[1]]
      
      # Store the centroid for the current territory
      territoryCentroidsInEachYear[[year]][[socialGroup]] <- c(mean(coords[, 1]), mean(coords[, 2]))
    }
  }
  
  return(territoryCentroidsInEachYear)
}

plotBadgerTerritories <- function(territoryCoordsInEachYear, years,
                                  scale=TRUE, background="white", map=NULL, ...){
  
  # Check if map provided for plotting
  if(is.null(map) == FALSE){
    plot(map, asp=1)
  }else{
    
    # Create an empty plot
    plot(x=NULL, y=NULL, xlim=territoryCoordsInEachYear$rangeX, ylim=territoryCoordsInEachYear$rangeY,
         bty="n", xaxt="n", yaxt="n", xlab="", ylab="", asp=1)
  }

  # Plot a background colour
  if(background != "white"){
    
    # Get the axisLimits
    axisLimits <- par()$usr
    
    # Plot rectangle to fill plotting background
    rect(xleft=axisLimits[1], ybottom=axisLimits[3], xright=axisLimits[2], ytop=axisLimits[4], 
         col=background)
  }
  
  # Examine each year
  for(year in years){
    
    # Skip 2006 - territories are all messed up
    if(year == 2006){
      next
    }
    
    # Convert year to character to act as key
    year <- as.character(year)

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
    socialGroup <- toupper(removeSep(as.character(data[row, column]), sep=" "))
    
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
