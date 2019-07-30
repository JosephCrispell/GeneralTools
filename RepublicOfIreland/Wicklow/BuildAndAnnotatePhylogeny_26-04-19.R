#### Load libraries ####

library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(rgdal) # Convert X and Y to lat longs and reading in shape files
library(OpenStreetMap) # Great tutorial here: https://www.r-bloggers.com/the-openstreetmap-package-opens-up/

#### Load sampling information ####

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

# Create a path variable
#path <- "/home/josephcrispell/Desktop/Research/RepublicOfIreland/Mbovis/"
path <- "J:\\WGS_Wicklow\\"

# Read in table that links original sequence ID to aliquot IDs
file <- paste0(path, "Mbovis_SamplingInfo_17-07-18.tsv")
linkTable <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep="\t")

# Read in the isolate metadata
file <- paste0(path, "IsolateSpeciesAndYear_26-04-19.csv")
metadata <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep=",")

# Read in the FASTA file
fastaFile <- paste(path, "vcfFiles/sequences_Prox-10_19-03-2019.fasta", sep="")
nSites <- getNSitesInFASTA(fastaFile)

# Read in the coverage information
coverageFile <- paste0(path, "vcfFiles/isolateCoverageSummary_DP-20_19-03-2019.txt")
coverage <- read.table(coverageFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)

#### Load the spatial data ####

### Read in the badger shape file
shapeFile <- paste0(path, "Raw_locations__all_species_May_2019\\badgers.shp")
badgerShapeFile <- readOGR(dsn=shapeFile)

# Extract the location data - X and Y coords
badgerLocations <- data.frame("Aliquot"=as.character(badgerShapeFile$ALIQUOT),
                              "X"=badgerShapeFile$X_COORD,
                              "Y"=badgerShapeFile$Y_COORD, stringsAsFactors=FALSE)

# Convert the X and Y coords to longitudes and latitudes
latLongs <- convertXYToLatLongs(x=badgerLocations$X, y=badgerLocations$Y)
badgerLocations$Lat <- latLongs$Latitude
badgerLocations$Long <- latLongs$Longitude

### Read in the deer shape file
shapeFile <- paste0(path, "Raw_locations__all_species_May_2019\\deer.shp")
deerShapeFile <- readOGR(dsn = shapeFile)

# Store the latitude and longitude
deerLocations <- data.frame("Aliquot"=as.character(deerShapeFile$ALIQUOT), 
                            "Lat"=deerShapeFile$LAT,
                            "Long"=deerShapeFile$LONG)

# Convert longitudes and latitudes to X and Y coordinates
coords <- convertLatLongsToXY(latitudes=deerLocations$Lat, longitudes=deerLocations$Long)
deerLocations$X <- coords$X
deerLocations$Y <- coords$Y

### Read in the cattle herd shape file
shapeFile <- paste0(path, "Raw_locations__all_species_May_2019\\farms_common.shp")
cattleShapeFile <- readOGR(dsn = shapeFile)

# Extract the polygon coords
landParcelCoords <- getPolygonCoords(cattleShapeFile)

# Extract the sampld herd information
herdInfo <- cattleShapeFile@data

#### Build phylogeny ####

# Build a phylogeny using RAxML
tree <- runRAXML(fastaFile, date="26-04-19", path, alreadyRun=TRUE, outgroup="\\>Ref-1997")

# Remove NI isolates and Reference
tree <- drop.tip(tree, tree$tip.label[grepl(tree$tip.label, pattern=">Ref-1997|>182-MBovis|>161-MBovis")])

# Edit the tip labels
tree$tip.label <- editTipLabels(tree$tip.label)

# Get the tip information (species and sampling date)
tipInfo <- getTipInfo(tree$tip.label, metadata, linkTable, coverage)

# Convert branch lengths to SNPs
tree$edge.length <- tree$edge.length * nSites

#### Identify and remove duplicates ####

# Count how many times each aliquot ID appears
aliquotCounts <- table(tipInfo$Aliquot)

# Identify the duplicated aliquots
duplicated <- names(aliquotCounts)[aliquotCounts > 1]

# Get the tip information for the duplicates
duplicatedTipInfo <- tipInfo[tipInfo$Aliquot %in% duplicated, ]

# Pick which tip to keep - note I checked phylogeny and these tips are identical on the phylogeny :-)
remove <- pickDuplicatedTipToRemove(duplicatedTipInfo)

# Drop one of each duplicated tip
tree <- drop.tip(tree, remove)

# Update the tip information
tipInfo <- getTipInfo(tree$tip.label, metadata, linkTable, coverage)

#### Plot the phylogeny ####

# Open an output PDF
outputPlotFile <- paste0(path, "MbovisAnnotatedPhylogeny_", date, ".pdf")
pdf(outputPlotFile)

# Get and set the margins
currentMar <- par()$mar
par(mar=c(4,0,0,10))

# Plot the phylogeny
plot.phylo(tree, show.tip.label=FALSE, edge.color="dimgrey", edge.width=4)

# Add tips coloured by species
tipShapesAndColours <- list("Badger"=c("red", 19), "Cow"=c("blue", 17), "Deer"=c("black", 15))
tiplabels(pch=getTipShapeOrColourBasedOnSpecies(tipInfo, tipShapesAndColours, which="shape"),
          col=getTipShapeOrColourBasedOnSpecies(tipInfo, tipShapesAndColours, which="colour"), cex=1.25)

# Add scale bar
addScaleBar(2)

# Add species legend
addSpeciesLegend(tipShapesAndColours)

# Add temporal sampling plot
plotTipSamplingDates(tipInfo)

# Reset the margins
par(mar=currentMar)

# Close the output pdf
dev.off()

#### Plot the sampling locations ####

# Add tip locations to tip info for badgers and deer
tipInfo <- addBadgerAndDeerLocationsToTipInfo(tipInfo, badgerLocations, deerLocations)

# Note the herd IDs associated with each cow
tipInfo <- addHerdIDsToTipInfo(tipInfo, herdInfo)

# Calculate the centroids of each herd's field
herdCentroids <- calculateLandParcelCentroids(landParcelCoords)

# Get the bounds of all the spatial data


# Get the bounds off all the spatial data and convert to lat longs
bounds <- calculateBoundsForAllSpatialData(badgerShapeFile@bbox, deerShapeFile@bbox,
                                           cattleShapeFile@bbox)
boundsLongsLats <- convertXYToLatLongs(x=bounds[1, ],
                                       y=bounds[2, ])

# Get a satellite image of the area
map <- getSatelliteImage(upperLeft=c(boundsLongsLats[2, 1], boundsLongsLats[1, 2]), 
                         lowerRight=c(boundsLongsLats[1, 1], boundsLongsLats[2, 2]))

# Open an output PDF
outputPlotFile <- paste0(path, "SamplingLocations_", date, ".pdf")
pdf(outputPlotFile)

# Plot the satellite image
plot(map)

# Add the cattle herd land parcels
plotHerdLandParcels(tipInfo, landParcelCoords, herdCentroids, col=rgb(0,1,0, 0.1),
                    border=rgb(0,0,0,1), plotHerds=FALSE, plotLines=FALSE)
plotHerdLandParcels(tipInfo, landParcelCoords, herdCentroids, 
                    connectingLineColour=rgb(0,0,1, 0.25), plotPolygons=FALSE,
                    plotHerds=FALSE)
plotHerdLandParcels(tipInfo, landParcelCoords, herdCentroids, plotPolygons=FALSE,
                    plotLines=FALSE)

# Plot the badger and deer sampling locations
points(tipInfo$X, tipInfo$Y, 
     pch=ifelse(tipInfo$Species == "Badger", 21, 22),
     bg=ifelse(tipInfo$Species == "Badger", rgb(1,0,0,0.75), rgb(0,0,0, 0.75)),
     col="white", xpd=TRUE, bty="n", yaxt="n", xaxt="n", cex=1)

dev.off()






#### Plot phylogeny linked to spatial locations ####



#### Make a note of the aliquot IDs that we don't have sequence data for yet ####

# Get the metadata for the aliquots that haven't been sequenced
notSequencedInfo <- metadata[metadata$Aliquot %in% tipInfo$Aliquot == FALSE, ]

# Print to file - just aliquot IDs
notSequencedFile <- paste0(path, "AliquotsNotSequenced_", date, ".csv")
write.table(notSequencedInfo[, "Aliquot"], file=notSequencedFile, quote=FALSE, sep=",", row.names=FALSE)

#### FUNCTIONS - Satellite image ####

getSatelliteImage <- function(upperLeft, lowerRight){
  
  # Get a satellite image of the area
  # Maps are initially put in a sperical mercator projection which is the standard for most (all?) map tiling systems
  map <- openmap(upperLeft,lowerRight, type="bing")
  
  # Convert from
  mapROIGrid <- openproj(map, projection=CRS("+init=epsg:29903"))
  
  return(mapROIGrid)
}

#### FUNCTIONS - Shape files ####

calculateBoundsForAllSpatialData <- function(badgerBounds, deerBounds, cattleBounds){
  
  # Initialise a dataframe to store the overall bounds
  bounds <- cattleBounds
  
  # Calculate the bounds across the badger, deer and cattle locations
  bounds[1, ] <- range(badgerBounds[1, ], deerBounds[1, ], cattleBounds[1, ])
  bounds[2, ] <- range(badgerBounds[2, ], deerBounds[2, ], cattleBounds[2, ])
  
  return(bounds)
}

calculateLandParcelCentroids <- function(landParcelCoords){
  
  # Initialise a list to store the centroids associated with each herd
  herdCentroids <- list()
  
  # Examine each of the herd
  for(herd in names(landParcelCoords)){
    
    # Initialise two vectors to store the X and Y coordinates of each field's centroid
    centroidsX <- c()
    centroidsY <- c()
    
    # Examine the field polygons for the current herd
    for(fieldIndex in seq_len(length(landParcelCoords[[herd]]))){
      
      # Get the polygon coordinates for the current field
      coords <- landParcelCoords[[herd]][[fieldIndex]]
      
      # Calculate the centroid for the current field
      centroidsX[fieldIndex] <- mean(coords$X)
      centroidsY[fieldIndex] <- mean(coords$Y)
    }
    
    # Calculate the overall centroid for the current herd
    overallX <- mean(centroidsX)
    overallY <- mean(centroidsY)
    
    # Store the calculated centroids
    herdCentroids[[herd]] <- list(
      "X"=overallX,
      "Y"=overallY,
      "CentroidXs"=centroidsX,
      "CentroidYs"=centroidsY)
  }
  
  return(herdCentroids)
}

plotHerdLandParcels <- function(tipInfo, landParcelCoords, herdCentroids,
                                connectingLineColour, plotPolygons=TRUE,
                                plotLines=TRUE, plotHerds=TRUE, ...){
  
  # Plot the cattle herd land parcels
  for(row in seq_len(nrow(tipInfo))){
    
    # Skip badgers and deer
    if(tipInfo[row, "Species"] != "Cow"){
      next
    }
    
    # Get herd code
    herdCode <- tipInfo[row, "HerdCode"]
    
    # Skip if herd code not available - send warning
    if(is.na(herdCode)){
      
      # Send warning message
      warning(paste0("No herd code available for aliquot: ", tipInfo[row, "Aliquot"]))
      next
    }
    
    # Examine the land parcel polygons for the current herd
    for(fieldIndex in seq_len(length(landParcelCoords[[herdCode]]))){
      
      # Plot the current polygon
      if(plotPolygons){
        polygon(landParcelCoords[[herdCode]][[fieldIndex]][, c("X", "Y")], ...)
      }
      
      # Connect current field to overall herd centre
      if(plotLines){
        points(x=c(herdCentroids[[herdCode]]$X, 
                   herdCentroids[[herdCode]]$CentroidXs[fieldIndex]),
               y=c(herdCentroids[[herdCode]]$Y, 
                   herdCentroids[[herdCode]]$CentroidYs[fieldIndex]),
               type="l", col=connectingLineColour)
      }
    }
    
    # Plot a symbol at the current herd's overall centre
    if(plotHerds){
      points(x=herdCentroids[[herdCode]]$X, y=herdCentroids[[herdCode]]$Y,
             pch=24, bg=rgb(0,0,1, 0.75), col="white")
    }
  }
}

addHerdIDsToTipInfo <- function(tipInfo, herdInfo){
  
  # Add an empty herd code column
  tipInfo$HerdCode <- NA
  
  # Examine each of the tips
  for(row in seq_len(nrow(tipInfo))){
    
    # Skip badgers and deer
    if(tipInfo[row, "Species"] != "Cow"){
      next
    }
    
    # Get the current tips aliquot
    aliquot <- tipInfo[row, "Aliquot"]
    
    # Get the row in the herd info table for the current aliquot
    rowInHerdInfo <- which(herdInfo$ALIQUOT == aliquot)
    
    # Store the herd ID
    if(length(rowInHerdInfo) > 0){
      tipInfo[row, "HerdCode"] <- as.character(herdInfo[rowInHerdInfo, "HERD_CODE"] - 1)
    }
  }
  
  return(tipInfo)
}

addBadgerAndDeerLocationsToTipInfo <- function(tipInfo, badgerLocations, deerLocations){
  
  # Add extra columns to store the location information
  tipInfo$X <- NA
  tipInfo$Y <- NA
  tipInfo$Longitude <- NA
  tipInfo$Latitude <- NA
  
  # Examine each of tips
  for(row in seq_len(nrow(tipInfo))){
    
    # Check if current tip associated with badger
    if(tipInfo[row, "Species"] == "Badger"){
    
      # Get the current tips aliquot
      aliquot <- tipInfo[row, "Aliquot"]
      
      # Get the row in the locations table for the current aliquot
      locationRow <- which(badgerLocations$Aliquot == aliquot)
      
      # Store the sampling location for the current badger
      tipInfo[row, "X"] <- badgerLocations[locationRow, "X"]
      tipInfo[row, "Y"] <- badgerLocations[locationRow, "Y"]
      tipInfo[row, "Longitude"] <- badgerLocations[locationRow, "Long"]
      tipInfo[row, "Latitude"] <- badgerLocations[locationRow, "Lat"]
    
    # Check if current tip associated with deer
    }else if(tipInfo[row, "Species"] == "Deer"){
      
      # Get the current tips aliquot
      aliquot <- tipInfo[row, "Aliquot"]
      
      # Get the row in the locations table for the current aliquot
      locationRow <- which(deerLocations$Aliquot == aliquot)
      
      # Store the sampling location for the current badger
      tipInfo[row, "X"] <- deerLocations[locationRow, "X"]
      tipInfo[row, "Y"] <- deerLocations[locationRow, "Y"]
      tipInfo[row, "Longitude"] <- deerLocations[locationRow, "Long"]
      tipInfo[row, "Latitude"] <- deerLocations[locationRow, "Lat"]
    }
  }
  
  return(tipInfo)
}

getPolygonCoords <- function(spatialDataFrame){
  
  # Got some good information from: https://stackoverflow.com/questions/29803253/r-extracting-coordinates-from-spatialpolygonsdataframe
  # Also from my blog: https://josephcrispell.github.io/BlogPosts/GetPolygonsFromShapeFile_11-10-17/GetPolygonsFromShapeFile_11-10-17.html
  # WHICH NEEDS UPDATED!!!!
  
  
  # Note the number of polygon sets - one for each ID
  nSets <- length(cattleShapeFile@polygons)
  
  # Initialise a list to store the IDs and coordinates of each polygon
  output <- list()
  
  # Loop through all the sets of polygons
  for(setIndex in seq_len(nSets)){
    
    # Get the ID of the current polygon set
    id <- cattleShapeFile@polygons[[setIndex]]@ID
    
    # Note the number of polygons in the current set
    nPolygons <- length(cattleShapeFile@polygons[[setIndex]]@Polygons)
    
    # Initialise a list to store the coordinates of each polygon in current set
    polygons <- list()
    
    # Examine each polygon in current set
    for(polygonIndex in seq_len(nPolygons)){
      
      # Get the coordinates for the current polygon
      coords <- cattleShapeFile@polygons[[setIndex]]@Polygons[[polygonIndex]]@coords
      colnames(coords) <- c("X", "Y")
      
      # Convert them to latitude and longitude
      latLongs <- convertXYToLatLongs(x=coords[, 1], y=coords[, 2])
      
      # Store the latitude and longitudes and X and Y coords
      polygons[[polygonIndex]] <- cbind(latLongs, coords)
    }
    
    # Store the polygon coordinates
    output[[id]] <- polygons
  }
  
  return(output)
}

convertXYToLatLongs <- function(x, y, grid="+init=epsg:29903"){
  
  # Create variables for holding the coordinate system types
  # see http://www.epsg.org/
  latLong <- "+init=epsg:4326"
  
  # Create a coordinates variable
  coords <- cbind(X = x, Y = y)
  
  # Create a SpatialPointsDataFrame
  spatialDF <- SpatialPoints(coords, proj4string = CRS(grid))
  
  # Convert the X and Ys to Longitudes and Latitudes
  spatialDFLatLongs <- spTransform(spatialDF, CRS(latLong))
  
  # Create a table to store the converted points
  output <- data.frame(Latitude=spatialDFLatLongs@coords[,"Y"],
                       Longitude=spatialDFLatLongs@coords[,"X"])
  
  # Return the lat longs
  return(output)
}

convertLatLongsToXY <- function(latitudes, longitudes, latLong="+init=epsg:4326"){
  
  # Create variables for holding the coordinate system types
  # see http://www.epsg.org/
  grid="+init=epsg:29903"
  
  # Create a coordinates variable
  coords <- cbind(Longitude=longitudes, Latitude=latitudes)
  
  # Create a SpatialPointsDataFrame
  spatialDF <- SpatialPoints(coords, proj4string = CRS(latLong))
  
  # Convert the latitude and longitudes to X and Y coordinates
  spatialDFLatLongs <- spTransform(spatialDF, CRS(grid))
  
  # Create a table to store the converted points
  output <- data.frame(X=spatialDFLatLongs@coords[,"Longitude"],
                       Y=spatialDFLatLongs@coords[,"Latitude"])
  
  # Return the X and Y coordinates
  return(output)
}

#### FUNCTIONS - Phylogeny ####

pickDuplicatedTipToRemove <- function(duplicatedTipInfo){
  
  # Initialise an array to store the tips to remove
  remove <- c()
  
  # Examine the tips associate with each aliquot
  for(aliquot in unique(duplicatedTipInfo$Aliquot)){
    
    # Get the tip info for the tips associated with the current aliquot
    tipInfo <- duplicatedTipInfo[duplicatedTipInfo$Aliquot == aliquot, ]
    
    # Pick the tip with the lowest coverage to remove
    remove <- c(remove, tipInfo[which(tipInfo$Coverage < max(tipInfo$Coverage)), "ID"])
  }
  
  return(remove)
}

addSpeciesLegend <- function(tipShapesAndColours, cex=1){
  
  # Get the plotting region dimensions = x1, x2, y1, y2 
  # (coordinates of bottom left and top right corners)
  dimensions <- par("usr")
  xLength <- dimensions[2] - dimensions[1]
  yLength <- dimensions[4] - dimensions[3]

  # Set the Y position
  yPos <- -0.08*yLength
  
  # Set the X start
  xStart <- 0.1*xLength
  
  # Set the x spacing
  xSpace <- 0.15*xLength
  
  # Loop through the tip options
  species <- names(tipShapesAndColours)
  for(i in seq_along(species)){
    
    # Plot a point for the current species
    points(x=xStart+(i*xSpace), y=yPos, pch=as.numeric(tipShapesAndColours[[species[i]]][2]),
           col=tipShapesAndColours[[species[i]]][1], xpd=TRUE, cex=cex)
    
    # Add a label
    text(x=xStart+(i*xSpace), y=yPos, labels=species[i], pos=4, xpd=TRUE, cex=cex)
  }
}

plotTipSamplingDates <- function(tipInfo){
  
  # Get the dimensions of the current plot
  dimensions <- par("usr")
  
  # Calculate the amount of space left to right of tree
  # par("plt") gives you proportion of device width taken up by plot: https://www.rstudio.com/wp-content/uploads/2016/10/how-big-is-your-graph.pdf
  lengthOfXAxis <- (dimensions[2] - dimensions[1]) / par("plt")[2]
  
  # Get all the information about the last phylogenetic tree plotted
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  
  # Add the tip coordinates into the tip information table
  tipInfo$X <- lastPP$xx[seq_len(nrow(tipInfo))]
  tipInfo$Y <- lastPP$yy[seq_len(nrow(tipInfo))]
  
  # Calculate the start and end of the plotting region on the X axis
  start <- max(lastPP$xx) + (0.05 * lengthOfXAxis)
  end <- lengthOfXAxis - (0.05 * lengthOfXAxis)
  
  # Convert the tip dates and calculate range
  tipInfo$Date <- as.Date(tipInfo$Date)
  dateRange <- range(tipInfo$Date, na.rm=TRUE)
  
  # Examine each of the tips
  for(row in seq_len(nrow(tipInfo))){
    
    # Skip tips with no date
    if(is.na(tipInfo[row, "Date"])){
      next
    }
    
    # Calculate the x position for the current tips date
    datePos <- getXPositionOfDate(start, end, tipInfo[row, "Date"], dateRange)
    
    # Plot the current date
    points(x=datePos, y=tipInfo[row, "Y"], pch=19, col="black", xpd=TRUE)
    
    # Plot a faint line to link date to tip
    points(x=c(tipInfo[row, "X"] + (0.01 * lengthOfXAxis), datePos), y=c(tipInfo[row, "Y"], tipInfo[row, "Y"]),
           type="l", lty=2, col=rgb(0,0,0, 0.25), xpd=TRUE)
  }
  
  # Add X axis bar
  points(x=c(getXPositionOfDate(start, end, as.Date("2014-03-01"), dateRange), 
             getXPositionOfDate(start, end, as.Date("2015-06-01"), dateRange)),
         y=c(-2,-2), xpd=TRUE, type="l", lwd=2)
  
  # Add ticks to x axis
  points(x=c(getXPositionOfDate(start, end, as.Date("2014-03-01"), dateRange),
             getXPositionOfDate(start, end, as.Date("2014-03-01"), dateRange)),
         y=c(-2, -3), xpd=TRUE, type="l", lwd=2)
  points(x=c(getXPositionOfDate(start, end, as.Date("2014-06-01"), dateRange),
             getXPositionOfDate(start, end, as.Date("2014-06-01"), dateRange)),
         y=c(-2, -3), xpd=TRUE, type="l", lwd=2)
  points(x=c(getXPositionOfDate(start, end, as.Date("2014-09-01"), dateRange),
             getXPositionOfDate(start, end, as.Date("2014-09-01"), dateRange)),
         y=c(-2, -3), xpd=TRUE, type="l", lwd=2)
  points(x=c(getXPositionOfDate(start, end, as.Date("2014-12-01"), dateRange),
             getXPositionOfDate(start, end, as.Date("2014-12-01"), dateRange)),
         y=c(-2, -3), xpd=TRUE, type="l", lwd=2)
  points(x=c(getXPositionOfDate(start, end, as.Date("2015-03-01"), dateRange),
             getXPositionOfDate(start, end, as.Date("2015-03-01"), dateRange)),
         y=c(-2, -3), xpd=TRUE, type="l", lwd=2)
  points(x=c(getXPositionOfDate(start, end, as.Date("2015-06-01"), dateRange),
             getXPositionOfDate(start, end, as.Date("2015-06-01"), dateRange)),
         y=c(-2, -3), xpd=TRUE, type="l", lwd=2)
  
  # Add tick labels
  text(x=c(getXPositionOfDate(start, end, as.Date("2014-03-01"), dateRange),
           getXPositionOfDate(start, end, as.Date("2014-06-01"), dateRange),
           getXPositionOfDate(start, end, as.Date("2014-09-01"), dateRange),
           getXPositionOfDate(start, end, as.Date("2014-12-01"), dateRange),
           getXPositionOfDate(start, end, as.Date("2015-03-01"), dateRange),
           getXPositionOfDate(start, end, as.Date("2015-06-01"), dateRange)),
       y=c(-4, -4, -4, -4, -4, -4),
       labels = c("MAR 14", "JUN 14", "SEP 14", "DEC 14", "MAR 15", "JUN 15"), 
       xpd=TRUE, cex=0.4)
}

getXPositionOfDate <- function(xStart, xEnd, date, dateRange){
  
  # Calculate the length of the X axis
  xLength <- xEnd - xStart
  
  # Calculate number of days in date range
  nDays <- as.numeric(dateRange[2] - dateRange[1])
  
  # Calculate the distance on the x axis for one day
  xLengthForOneDay <- xLength / nDays
  
  # Calculate number of days between date and start of dates
  nDaysToCurrentDate <- as.numeric(date - dateRange[1])
  
  # Calculate X position of current date
  xPosition <- xStart + (nDaysToCurrentDate * xLengthForOneDay)
  
  return(xPosition)
}

addScaleBar <- function(scaleSize, cex=1){
  
  # Get the plotting region dimensions = x1, x2, y1, y2 
  # (coordinates of bottom left and top right corners)
  dimensions <- par("usr")
  xLength <- dimensions[2] - dimensions[1]
  yLength <- dimensions[4] - dimensions[3]
  
  # Add Scale bar
  xPad <- 0.1 * xLength

  points(x=c(dimensions[1] + xPad, dimensions[1] + xPad + scaleSize), 
         y=c(dimensions[3] - (0.01 * yLength), dimensions[3] - (0.01 * yLength)),
         type="l", lwd=3, xpd=TRUE)
  if(scaleSize == 1){
    text(x=dimensions[1] + xPad + (0.5*scaleSize), y=dimensions[3] - (0.04 * yLength),
         labels=paste0("~ ", scaleSize, " SNP"), cex=cex, xpd=TRUE)
  }else{
    text(x=dimensions[1] + xPad + (0.5*scaleSize), y=dimensions[3] - (0.04 * yLength),
         labels=paste0("~ ", scaleSize, " SNPs"), cex=cex, xpd=TRUE)
  }
}

getTipShapeOrColourBasedOnSpecies <- function(tipInfo, tipShapesAndColours, which){
  
  # Initialise a vector to store the shapes or colours
  output <- c()
  
  # Examine each tip
  for(row in seq_len(nrow(tipInfo))){
    
    # Check if Species available
    if(is.na(tipInfo[row, "Species"])){
      
      # Check if wanting shape or colour
      if(which == "shape"){
        output[row] <- as.numeric(tipShapesAndColours[["NA"]][2])
      }else if(which == "colour"){
        output[row] <- tipShapesAndColours[["NA"]][1]
      }else{
        cat(paste("Error! Option", which, "not recognised!"))
      }
    
    # If species available assign appropriate colour or shape
    }else{
      
      # Check if wanting shape or colour
      if(which == "shape"){
        output[row] <- as.numeric(tipShapesAndColours[[tipInfo[row, "Species"]]][2])
      }else if(which == "colour"){
        output[row] <- tipShapesAndColours[[tipInfo[row, "Species"]]][1]
      }else{
        cat(paste("Error! Option", which, "not recognised!"))
      }
    }
  }
  
  return(output)
}

getTipInfo <- function(tipLabels, metadata, linkTable, coverage){
  
  # Initialise a dataframe to store the tip information
  tipInfo <- data.frame(ID=tipLabels, Species=NA, Date=NA, Aliquot=NA, Coverage=NA, stringsAsFactors=FALSE)
  
  # Examine each of the tips
  for(index in seq_along(tipLabels)){
    
    # Initialise variables to store the tip's information
    aliquotCode <- NA
    species <- NA
    date <- NA
    quality <- NA
    
    # Check if the current tip is associated with the original dataset
    if(grepl(tipLabels[index], pattern="-MBovis")){
      
      # Get the sequence number from the curren tip label
      sequenceNumber <- strsplit(tipLabels[index], split="-")[[1]][1]
      
      # Find the row in the link table
      row <- which(linkTable$Isolate.Code == sequenceNumber)
      
      # Get the current tips aliquot code
      if(length(row) != 0){
        aliquotCode <- linkTable[row, "Aliquot"]
      }else{
        cat(paste("Error for old batch. Couldn't find sequence number: ", sequenceNumber, " ", tipLabels[index], "\n"))
      }
    
    # Get information from most recent sequencing run
    }else{
      
      # Get the second part of the tip label - looks like an aliuot label without 00s
      aliquotCodePart <- strsplit(tipLabels[index], split="-")[[1]][2]
      
      # Find row that matches above part
      row <- which(grepl(metadata$Aliquot, pattern=aliquotCodePart))
      
      # Get the full aliquot code
      if(length(row) == 1){
        aliquotCode <- metadata[row, "Aliquot"]
      }else{
        cat(paste("Error for new batch. Couldn't find aliquot part: ", aliquotCodePart, " (found ", length(row), " matches)\n\n"))
      }
    }
    
    # If aliquot code available then get the tip species and sampling year
    if(is.na(aliquotCode) == FALSE){
      
      # Get the row in the metadata table for the current aliquot code
      row <- which(metadata$Aliquot == aliquotCode)
      
      # Note the species and convert multiple "Cow labels to single "Cow label
      species <- metadata[row, "Species"]
      if(species %in% c("Heifer", "Steer", "Calf", "Bull", "Bovine")){
        species <- "Cow"
      }
      
      # Note the sampling date
      if(metadata[row, "Received.Test.date"] != ""){
        date <- as.character(as.Date(metadata[row, "Received.Test.date"], format="%d/%m/%Y"))
      }
    }
    
    # Check species isn't nothing
    if(is.na(species) == FALSE && species == ""){
      species <- NA
    }
    
    # Get the isolate's genome coverage information
    tipInfo[index, "Coverage"] <- coverage[grepl(coverage$IsolateID, pattern=tipLabels[index]), "PercentageCoverage"]
    
    # Store the current tips information
    tipInfo[index, "Species"] <- species
    tipInfo[index, "Date"] <- date
    tipInfo[index, "Aliquot"] <- aliquotCode
  }
  
  return(tipInfo)
}

editTipLabels <- function(tipLabels){
  
  # Initialise a vector to store the new labels
  output <- c()
  
  # Examine each tip label
  for(label in tipLabels){
    
    # Remove the ">" prefix
    label <- substr(label, 2, nchar(label))
    
    # Split the label and retain first part
    label <- strsplit(label, split="_")[[1]][1]
    
    # Store the new label
    output[length(output) + 1] <- label
  }
  
  return(output)
}

getNSitesInFASTA <- function(fastaFile){
  
  # Open a connection to a file to read (open="r")
  connection <- file(fastaFile, open="r")
  
  # Get first line of file
  firstLine <- readLines(connection, n=1)
  
  # Close file connection
  close(connection)
  
  # Get the number of sites used in the FASTA file from first line
  nSites <- as.numeric(strsplit(firstLine, " ")[[1]][2])
  
  return(nSites)
}

runRAXML <- function(fastaFile, date, path, nBootstraps=100, nThreads=6, alreadyRun=FALSE, outgroup=NULL, model="GTRCAT"){
  
  # Create a directory for the output file
  directory <- paste(path, "RAxML_", date, sep="")
  suppressWarnings(dir.create(directory))
  
  # Set the Working directory - this will be where the output files are dumped
  setwd(directory)
  
  # Build analysis name
  analysisName <- paste("RaxML-R_", date, sep="")
  
  # Check if already Run and just want to retrieve tree
  if(alreadyRun == FALSE){
    
    # Build the command
    seeds <- sample(1:100000000, size=2, replace=FALSE) # For parsimony tree and boostrapping
    
    if(is.null(outgroup)){
      command <- paste("raxmlHPC", 
                       " -f a", # Algorithm: Rapid boostrap inference
                       " -N ", nBootstraps,
                       " -T ", nThreads,
                       " -m ", model, " -V", # -V means no rate heterogenity
                       " -p ", seeds[1], " -x ", seeds[2], # Parsimony and boostrapping seeds
                       " -n ", analysisName,
                       " -s ", fastaFile, sep="")
    }else{
      command <- paste("raxmlHPC", 
                       " -f a", # Algorithm: Rapid boostrap inference
                       " -N ", nBootstraps,
                       " -T ", nThreads,
                       " -m ", model, " -V", # -V means no rate heterogenity
                       " -p ", seeds[1], " -x ", seeds[2], # Parsimony and boostrapping seeds
                       " -n ", analysisName,
                       " -s ", fastaFile, 
                       " -o ", outgroup, sep="")
    }
    
    system(command, intern=TRUE)
  }
  
  # Get the tree and read it in
  treeBS <- getTreeFileWithSupportValues(analysisName)
  
  return(treeBS)
}

getTreeFileWithSupportValues <- function(analysisName){
  
  # Get files in current working directory
  files <- list.files()
  
  # Select the tree file with BS support values
  treeBSFile <- files[grepl(files, pattern=paste("RAxML_bipartitions[.]", analysisName, sep="")) == TRUE]
  
  # Open the file
  treeBS <- read.tree(treeBSFile)
  
  return(treeBS)
}
