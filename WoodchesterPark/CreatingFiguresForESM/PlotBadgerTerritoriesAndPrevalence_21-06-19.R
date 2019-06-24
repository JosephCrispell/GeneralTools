#### Preparation ####

# Load libraries
library(rgdal) # For reading in shape files

# Set the path variable
path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/"

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

#### Plot badger territories through time ####

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

# Plot the badger territories
territoryPlotFile <- paste0(path, "ESM_Figures/BadgerTerritories_", date, ".pdf")
plotBadgerTerritories(territoryCoordsInEachYear, years, file=territoryPlotFile, lwd=3, border=rgb(0,0,0, 1))


#### Plot prevalence through time ####

# Read in the badger capture data
file <- paste0(path, "BadgerCaptureData/WP_CaptureData_Consolidated_24-06-2019.txt")
captureData <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

# Count the badgers in each status in each year - USING END OF YEAR
counts <- countBadgersInEachStatusInEachYear(captureData, yearRange=c(1990, 2010))

# Plot the counts
prevalencePlotFile <- paste0(path, "ESM_Figures/BadgerBTBPrevalence_", date, ".pdf")
plotBadgerPrevalence(counts, brokenDown=TRUE, file=prevalencePlotFile)

#### FUNCTIONS ####

countBadgersInEachStatusInEachYear <- function(captureData, yearRange){
  
  # Create a vector years we want counts for
  years <- yearRange[1]:yearRange[2]

  # Initialise a table to store the counts
  counts <- data.frame(matrix(0, nrow=length(years), ncol=5))
  colnames(counts) <- c("Year", "Negative", "Exposed", "Excretor", "Superexcretor")
  counts$Year <- years

  # Examine the capture data for each year
  for(row in seq_len(nrow(captureData))){
    
    # Get the capture dates
    captureDates <- strsplit(captureData[row, "@CaptureDates"], split=";")[[1]]
    captureDates <- as.Date(captureDates, format="%d-%m-%Y")

    # Get the disease status dates
    statuses <- strsplit(captureData[row, "@Statuses"], split=";")[[1]]
    
    # Calculate the period spent in each disease status
    statusPeriods <- calculatePeriodSpentInEachStatus(captureDates, statuses)
    
    # Examine each status and the period spent it
    for(status in c("Negative", "Exposed", "Excretor", "Superexcretor")){
      
      # Are dates available for the current status
      if(length(statusPeriods[[status]]) < 2){
        next
      }
      
      # Examine each of our years of interest
      for(yearIndex in seq_along(years)){

        # Create a date to mark end of currnt year
        endOfYear <- paste0(years[yearIndex], "-12-31")
        
        # Check if current year in in period for current status
        if(endOfYear > statusPeriods[[status]][1] && endOfYear < statusPeriods[[status]][2]){

          # Add to count for current status in current year
          counts[yearIndex, status] <- counts[yearIndex, status] + 1
        }
      }
    }
  }
  
  # Calculate the population total
  counts$Total <- counts$Negative + counts$Exposed + counts$Excretor + counts$Superexcretor
  
  return(counts)
}

calculatePeriodSpentInEachStatus <- function(captureDates, statuses){
  
  # Initialise a list to store the periods spent in each status
  periods <- list("Negative"=c(), "Exposed"=c(), "Excretor"=c(), "Superexcretor"=c())
  
  # Store the first date
  status <- statuses[1]
  periods[[status]][1] <- as.character(captureDates[1])

  # Examine the capture dates and disease statuses
  for(i in seq_along(statuses)){
    
    # Skip first index
    if(i == 1){
      next;
    }
    
    # Check if status has changed
    if(status != statuses[i]){
      
      # Set end date of period status
      periods[[status]][2] <- as.character(captureDates[i] - 1)
      
      # Set start date of new status
      status <- statuses[i]
      periods[[status]][1] <- as.character(captureDates[i])
    }
  }
  
  # Store the last date captured
  periods[[status]][2] <- as.character(captureDates[length(captureDates)])
  
  return(periods)
}

plotBadgerPrevalence <- function(captureCounts, xLim=NULL, brokenDown=FALSE, file=NULL){
  
  # Open a pdf file if requested
  if(is.null(file) == FALSE){
    pdf(file)
  }
  
  # Get and set the plotting margins
  currentMar <- par("mar")
  par(mar=c(5.1,5.1,1,1))
  
  # Check if X axis (years) is limited
  if(is.null(xLim) == FALSE){
    captureCounts <- captureCounts[captureCounts$Year >= xLim[1] & captureCounts$Year <= xLim[2], ]
  }
  
  # Plot the total counts
  plot(x=captureCounts$Year, y=captureCounts$Negative, las=1, type="o", pch=19, lwd=2, bty="n",
       ylim=c(0, max(captureCounts$Negative)), xlab="Year", ylab="", cex.axis=1.5,
       cex.lab=1.75, xpd=TRUE)
  mtext(side=2, text="Number Captured", cex=1.75, line=3.5)

  # Add infected counts - check if wanyt broekn into categories
  if(brokenDown){
    
    # Plot counts in each infection category
    points(x=captureCounts$Year, y=captureCounts$Exposed,
           type="o", col="brown", pch=19, lwd=2)
    points(x=captureCounts$Year, y=captureCounts$Excretor,
           type="o", col="green", pch=19, lwd=2)
    points(x=captureCounts$Year, y=captureCounts$Superexcretor,
           type="o", col="blue", pch=19, lwd=2)

    # Add legend
    legend("topright", legend=c("Negative", "Exposed", "Excretor", "Super excretor"),
           text.col=c("black", "brown", "green", "blue"), bty="n", cex=1.5)

  }else{
    
    # Plot total infected
    points(x=captureCounts$Year, y=captureCounts$Exposed + captureCounts$Excretor + captureCounts$`Super excretor`,
           type="o", col="green", pch=19)
    
    # Add legend
    legend("topright", legend=c("Negative", "Positive"),
           text.col=c("black", "green"), bty="n", cex=1.5)
  }
  
  # Reset the plotting margins
  par(mar=currentMar)
  
  # Close the pdf if requested
  if(is.null(file) == FALSE){
    dev.off()
  }
}

plotBadgerTerritories <- function(territoryCoordsInEachYear, years, sleep=NULL, file=NULL, ...){
  
  # Open a pdf file if requested
  if(is.null(file) == FALSE){
    pdf(file)
  }
  
  # Get and set the plotting margins
  currentMar <- par("mar")
  par(mar=c(0,0,0,0))
  
  # Create an empty plot
  plot(x=NULL, y=NULL, xlim=territoryCoordsInEachYear$rangeX, ylim=territoryCoordsInEachYear$rangeY,
       bty="n", xaxt="n", yaxt="n", xlab="", ylab="", asp=1)
  
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
