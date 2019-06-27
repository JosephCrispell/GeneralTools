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

#### Get satellite image of Woodchester Park ####

# Calculate the plotting area boundaries in latitude and longitudes
coordRanges <- convertXYToLatLongs(x=territoryCoordsInEachYear$rangeX, y=territoryCoordsInEachYear$rangeY)
rangeLat <- coordRanges[, "Latitude"]
rangeLong <- coordRanges[, "Longitude"]

# Get a satellite image of the area
# Maps are initially put in a sperical mercator projection which is the standard for most (all?) map tiling systems
map <- openmap(upperLeft=c(rangeLong[1], rangeLat[1]), 
               lowerRight=c(rangeLong[2], rangeLat[2]),
               type="bing")

# Convert from
mapUKGrid <- openproj(map, projection=CRS("+init=epsg:27700"))

#### Plot all territories in single figure ####

# Open a pdf
territoryPlotFile <- paste0(path, "ESM_Figures/BadgerTerritories/BadgerTerritories_", date, ".pdf")
pdf(territoryPlotFile)

# Plot the badger territories
plotBadgerTerritories(territoryCoordsInEachYear, years, lwd=2,
                      border=rgb(1,0,0, 0.5), map=mapUKGrid)

#### Plot prevalence through time ####

# # Read in the badger capture data
# file <- paste0(path, "BadgerCaptureData/WP_CaptureData_Consolidated_24-06-2019.txt")
# captureData <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
# 
# # Count the badgers in each status in each year - USING END OF YEAR
# counts <- countBadgersInEachStatusInEachYear(captureData, yearRange=c(1990, 2010))
# 
# # Plot the counts
# prevalencePlotFile <- paste0(path, "ESM_Figures/BadgerBTBPrevalence_", date, ".pdf")
# plotBadgerPrevalence(counts, brokenDown=TRUE, file=prevalencePlotFile)

# Read in the counts information - badgers in each status in each social group
file <- paste0(path, "BadgerCaptureData/InfectionCategoryCounts_2000-2011_10-08-2017.csv")
counts <- getCountTablesFromFileLinesYears(file)

# Plot Woodchester Map without any territories
plotBadgerTerritories(territoryCoordsInEachYear, years, lwd=2,
                      border=rgb(1,0,0, 0), map=mapUKGrid)

# Plot the territory outlines for 2001 and colour by prevalence
for(year in 2000:2011){

  plotTerritoriesFromYear(territoryCoordsInEachYear, year=year,
                          alphas=calculatePrevalenceInEachGroup(counts, year=year),
                          fill="red", lwd=2, map=mapUKGrid)
}

# Close the PDF file
dev.off()

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
  output <- data.frame(Longitude=spatialDFLatLongs@coords[,"Northing"],
                       Latitude=spatialDFLatLongs@coords[,"Easting"])
  
  # Return the lat longs
  return(output)
}

plotBadgerTerritories <- function(territoryCoordsInEachYear, years, sleep=NULL, file=NULL, map=NULL, ...){
  
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
  axisLimits <- par()$usr
  xLength <- axisLimits[2] - axisLimits[1]
  yLength <- axisLimits[4] - axisLimits[3]
  xPad <- 0.08*xLength
  points(x=c(axisLimits[2] - xPad - 1000, axisLimits[2] - xPad), y=c(axisLimits[3] + 0.1*yLength, axisLimits[3] + 0.1*yLength),
         type="l", lwd=4)
  text(x=axisLimits[2] - xPad - 500, y=axisLimits[3] + 0.07*yLength, labels="1KM", cex=2)
  
  # Reset the plotting margins
  par(mar=currentMar)
  
  # Close the pdf if requested
  if(is.null(file) == FALSE){
    dev.off()
  }
}

convertXYToSphericalMercatorProjection <- function(x, y, grid="+init=epsg:27700"){
  
  # Create variables for holding the coordinate system types
  # see http://www.epsg.org/
  projection <- "+init=epsg:3857"
  
  # Create a coordinates variable
  coords <- cbind(X = x, Y = y)
  
  # Create a SpatialPointsDataFrame
  spatialDF <- SpatialPoints(coords, proj4string = CRS(grid))
  
  # Convert the Eastings and Northings to Latitude and Longitude
  spatialDF<- spTransform(spatialDF, CRS(projection))
  
  # Create a table to store the converted points
  output <- data.frame(X=spatialDF@coords[,"X"], Y=spatialDF@coords[,"Y"])
  
  # Return the lat longs
  return(output)
}

convertTerritoryCoordinatesToSphericalMercator <- function(territoryCoordsInEachYear){
  
  # Initialise a list to store the coordinates as latitude and longitude
  output <- list()
  
  # Convert the range values to latitude and longitude
  coordRanges <- convertXYToSphericalMercatorProjection(x=territoryCoordsInEachYear$rangeX, y=territoryCoordsInEachYear$rangeY)
  output$rangeX <- coordRanges[, "X"]
  output$rangeY <- coordRanges[, "Y"]
  
  # Examine each year
  nYears <- length(territoryCoordsInEachYear) - 2
  for(year in names(territoryCoordsInEachYear)[1:nYears]){
    
    # Initialise a list for the current year
    output[[year]] <- list()
    
    # Examine every social group
    for(socialGroup in names(territoryCoordsInEachYear[[year]])){
      
      # Get the coordinates for the current social gropu in the current year
      coords <- territoryCoordsInEachYear[[year]][[socialGroup]][[1]]
      
      # Convert the coodinates to latitudes and longitudes
      coords <- convertXYToSphericalMercatorProjection(x=coords[, 1], y=coords[, 2])
      
      # Store the latitude and longitudes
      output[[year]][[socialGroup]][[1]] <- coords
    }
  }
  
  return(output)
}

plotTerritoriesFromYear <- function(territoryCoordsInEachYear, year, fill="black", alphas=NULL, file=NULL, map=NULL, ...){
  
  # Open a pdf file if requested
  if(is.null(file) == FALSE){
    pdf(file)
  }
  
  # Get and set the plotting margins
  currentMar <- par("mar")
  par(mar=c(0,0,0,0))
  
  # Convert the year to character
  year <- as.character(year)
  
  # Check if map provided for plotting
  if(is.null(map) == FALSE){
    plot(map, asp=1)
  }else{
    
    # Create an empty plot
    plot(x=NULL, y=NULL, xlim=territoryCoordsInEachYear$rangeX, ylim=territoryCoordsInEachYear$rangeY,
         bty="n", xaxt="n", yaxt="n", xlab="", ylab="", asp=1)
  }
  
  # Add plot label
  axisLimits <- par()$usr
  xLength <- axisLimits[2] - axisLimits[1]
  yLength <- axisLimits[4] - axisLimits[3]
  text(x=axisLimits[2] - 0.1*xLength, y=axisLimits[4] - 0.1*yLength, labels=year, cex=2)
  
  # Examine each social group
  for(socialGroup in names(territoryCoordsInEachYear[[year]])){
    
    # Get the coordinates for current social group's polygon
    coords <- territoryCoordsInEachYear[[year]][[socialGroup]][[1]]
    
    # Convert the social group to upper case
    socialGroup <- toupper(socialGroup)
    
    # Plot the coordinates - check if transparency wanted
    if(is.null(alphas[[socialGroup]])){
      polygon(coords, ...)
    }else{
      polygon(coords, col=setAlpha(fill, alphas[[socialGroup]]), ...)
    }
  }
  
  # Add scale
  axisLimits <- par()$usr
  xLength <- axisLimits[2] - axisLimits[1]
  yLength <- axisLimits[4] - axisLimits[3]
  xPad <- 0.08*xLength
  points(x=c(axisLimits[2] - xPad - 1000, axisLimits[2] - xPad), y=c(axisLimits[3] + 0.1*yLength, axisLimits[3] + 0.1*yLength),
         type="l", lwd=4)
  text(x=axisLimits[2] - xPad - 500, y=axisLimits[3] + 0.07*yLength, labels="1KM", cex=2)
  
  # Open a pdf file if requested
  if(is.null(file) == FALSE){
    dev.off()
  }
  
  # Reset the plotting margins
  par(mar=currentMar)
}

calculatePrevalenceInEachGroup <- function(counts, year){
  
  # Get the number of badgers in each status in each group for the current year
  nBadgers <- getNumberInEachStatusInEachGroupForGivenYear(counts, year)
  
  # Remove the NA group
  nBadgers <- nBadgers[nBadgers$Group != "NA", ]
  
  # Calculate total number in each group
  nBadgers$Total <- nBadgers$Negative + nBadgers$Exposed + nBadgers$Excretor + nBadgers$Superexcretor
  
  # Calculate the prevalence within group
  nBadgers$Prevalence <- (nBadgers$Exposed + nBadgers$Excretor + nBadgers$Superexcretor) / nBadgers$Total

  # Create a list reporting prevalence of each group
  prevalence <- list()
  for(row in 1:nrow(nBadgers)){
    prevalence[[nBadgers[row, "Group"]]] <- ifelse(is.nan(nBadgers[row, "Prevalence"]), 0, nBadgers[row, "Prevalence"])
  }
  
  return(prevalence)
}

getNumberInEachStatusInEachGroupForGivenYear <- function(counts, year){
  
  # Get a vector of the social group names
  socialGroups <- colnames(counts$Negative)[-which(colnames(counts$Negative) %in% c("Years", "Total"))]
  
  # Identify the row of the current year
  row <- which(counts$Negative$Years == year)
  
  # Initialise a data frame to store the number of badgers in each status
  output <- data.frame(matrix(0, nrow=length(socialGroups), ncol=5))
  colnames(output) <- c("Group", "Negative", "Exposed", "Excretor", "Superexcretor")
  output$Group <- socialGroups
  rownames(output) <- socialGroups
  
  # Examine the counts for each status
  for(status in c("Negative", "Exposed", "Excretor", "Superexcretor")){
    
    # Get the counts for the current year
    values <- as.numeric(counts[[status]][row, socialGroups])
    
    # Store the counts
    output[socialGroups, status] <- values
  }
  
  return(output)
}

getCountTablesFromFileLinesYears <- function(fileName){
  
  connection <- file(fileName, open='r')
  fileLines <- readLines(connection)
  close(connection)
  
  counts <- list()
  counts[["Negative"]] <- data.frame(Years=rep(0, length(fileLines) - 1), stringsAsFactors=FALSE)
  counts[["Exposed"]] <- data.frame(Years=rep(0, length(fileLines) - 1), stringsAsFactors=FALSE)
  counts[["Excretor"]] <- data.frame(Years=rep(0, length(fileLines) - 1), stringsAsFactors=FALSE)
  counts[["Superexcretor"]] <- data.frame(Years=rep(0, length(fileLines) - 1), stringsAsFactors=FALSE)
  
  # Get the social group names from the header
  socialGroups <- strsplit(fileLines[1], "\t")[[1]][-1]
  
  # Read in the tables
  for(row in 2:length(fileLines)){
    cols <- strsplit(fileLines[row], "\t")[[1]]
    
    counts[["Negative"]][row-1, "Years"] <- as.numeric(cols[1])
    counts[["Exposed"]][row-1, "Years"] <- as.numeric(cols[1])
    counts[["Excretor"]][row-1, "Years"] <- as.numeric(cols[1])
    counts[["Superexcretor"]][row-1, "Years"] <- as.numeric(cols[1])
    
    for(col in 2:length(cols)){
      
      parts <- strsplit(cols[col], ":")[[1]]
      counts[["Negative"]][row-1, col] <- parts[1]
      counts[["Exposed"]][row-1, col] <- parts[2]
      counts[["Excretor"]][row-1, col] <- parts[3]
      counts[["Superexcretor"]][row-1, col] <- parts[4]
      
    }
  }
  
  # Add the column names
  colnames(counts[["Negative"]]) <- c("Years", socialGroups)
  colnames(counts[["Exposed"]]) <- c("Years", socialGroups)
  colnames(counts[["Excretor"]]) <- c("Years", socialGroups)
  colnames(counts[["Superexcretor"]]) <- c("Years", socialGroups)
  
  # Add total column
  counts[["Negative"]]$Total <- calculateTotalForEachYear(counts[["Negative"]])
  counts[["Exposed"]]$Total <- calculateTotalForEachYear(counts[["Exposed"]])
  counts[["Excretor"]]$Total <- calculateTotalForEachYear(counts[["Excretor"]])
  counts[["Superexcretor"]]$Total <- calculateTotalForEachYear(counts[["Superexcretor"]])
  
  return(counts)
}

calculateTotalForEachYear <- function(table){
  
  # Initialise a vector to store the sum of counts for each year
  sums <- c()
  
  # Examine each year
  for(row in seq_len(nrow(table))){
    
    sums[row] <- sum(as.numeric(table[row, 2:ncol(table)]))
  }
  
  return(sums)
}

plotBadgerPrevalence <- function(counts, file=NULL){
  
  # Open a pdf file if requested
  if(is.null(file) == FALSE){
    pdf(file)
  }
  
  # Get and set the plotting margins
  currentMar <- par("mar")
  par(mar=c(5.1,5.1,1,1))
  
  # Plot the negative counts
  plot(x=counts$Negative$Years, y=counts$Negative$Total, las=1, type="o", pch=19, lwd=2, bty="n",
       ylim=c(0, max(counts$Negative$Total)), xlab="Year", ylab="", cex.axis=1.5,
       cex.lab=1.75, xpd=TRUE)
  mtext(side=2, text="Number Captured", cex=1.75, line=3.5)

  # Plot counts in each infection category
  points(x=counts$Exposed$Years, y=counts$Exposed$Total,
         type="o", col="brown", pch=19, lwd=2)
  points(x=counts$Excretor$Years, y=counts$Excretor$Total,
         type="o", col="green", pch=19, lwd=2)
  points(x=counts$Superexcretor$Years, y=counts$Superexcretor$Total,
         type="o", col="blue", pch=19, lwd=2)

  # Add legend
  legend("topright", legend=c("Negative", "Exposed", "Excretor", "Super excretor"),
         text.col=c("black", "brown", "green", "blue"), bty="n", cex=1.5)

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
