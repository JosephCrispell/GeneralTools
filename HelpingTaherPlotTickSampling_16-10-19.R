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
# Removed location columns as had characters interferring with reading in
tickCountsFile <- paste0(path, "HelpingTaher/Taher_AllSitesAndPathogens_16-10-19.tsv")
tickCounts <- read.table(tickCountsFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Calculate prevalence
tickCounts$Prevalence <- tickCounts$No..of.qPCR.positives / tickCounts$Ticks.screened
tickCounts$Prevalence[is.nan(tickCounts$Prevalence)] <- 0

#### Plot the sampling ####

# Open a pdf
file <- paste0(substr(tickCountsFile, 1, nchar(tickCountsFile)-4), ".pdf")
pdf(file, width=14)

# Get and set the plotting margins
currentMar <- par()$mar
par(mar=c(0,0,2,0), mfrow=c(1,2))

# Plot the polygons from Ireland
plotPolygons(irelandPolygons, lwd=1.5, border=rgb(0.4,0.4,0.4, 0.1), col=rgb(0,0,1, 0.1),
             main="Sites sampled for presence of ticks from 2017 to 2019", label="a")

# Plot the locations that were sampled
points(x=tickCounts$Long, y=tickCounts$Lat, cex=1, pch=19, 
       col=ifelse(tickCounts$No..of.collected.ticks > 0, rgb(1,0,0, 0.4), rgb(0,0,0, 0.4)))

# Add a legend
legend("right", legend=c("Present", "Absent"), pch=19, col=c(rgb(1,0,0, 0.4), rgb(0,0,0, 0.4)),
       bty="n")

# Plot the Ireland county polygons
plotPolygons(irelandPolygons, lwd=1.5, border=rgb(0.4,0.4,0.4, 0.1), col=rgb(0,0,1, 0.1),
             main="Prevelance of Borrelia at screened sites", label="b")

# Plot the tick counts
addPoints(xCoords=tickCounts$Long, yCoords=tickCounts$Lat, pch=21,
          cex=(tickCounts$Ticks.screened/max(tickCounts$Ticks.screened))*4,
          bg=rgb(1,0,0, tickCounts$Prevalence),
          col="black", lty=2)

# Add a scale for number of ticks
legend("right", title="No. ticks screened", legend=c(100, 75, 50, 25, 10), pch=21, 
       col="black", bty="n", pt.cex=(c(100, 75, 50, 25, 10)/max(tickCounts$Ticks.screened))*4)

# Add a scale for prevalence
legend("bottomright", title="Prevalence         ", legend=paste0(c(50, 40, 30, 20, 10, 0), "%"), pch=21, 
       pt.bg=c(rgb(1,0,0, 0.5), rgb(1,0,0, 0.4), rgb(1,0,0, 0.3), 
            rgb(1,0,0, 0.2), rgb(1,0,0, 0.1), rgb(1,0,0, 0)),
       bty="n", pt.cex=2, col="black")

# Reset the plotting margins
par(mar=currentMar, mfrow=c(1,1))

# Close the pdf
dev.off()

#### FUNCTIONS ####

plotPolygons <- function(polygons, xLim=NULL, yLim=NULL, main="", asp=NULL, label=NULL, ...){
  
  # Get the axis limits of uk - enough to include ROI
  if(is.null(xLim)){
    
    xLim <- polygons$bbox[1, ]
    yLim <- polygons$bbox[2, ]
  }
  
  # Create an empty plot
  if(is.null(asp)){
    plot(x=NULL, y=NULL, xlim=xLim, ylim=yLim,
         bty="n", xaxt="n", yaxt="n", xlab="", ylab="", main=main)
  }else{
    plot(x=NULL, y=NULL, xlim=xLim, ylim=yLim,
         bty="n", xaxt="n", yaxt="n", xlab="", ylab="", asp=asp, main=main)
  }
  
  # Add a label if requested
  if(is.null(label) == FALSE){
    mtext(label, side=3, line=1, at=axisLimits[1], cex=2.5)
  }
  
  # Examine each set of polygons
  for(setID in names(polygons)){

    # Skip the bounding box
    if(setID == "bbox"){
      next
    }
    
    # Examine each polygon within current set
    for(polygonIndex in 1:length(polygons[[setID]])){
      
      # Plot the current polygon
      polygon(polygons[[setID]][[polygonIndex]], ...)
    }
  }
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
  
  # Store the bounds
  output[["bbox"]] <- spatialDataFrame@bbox
  
  return(output)
}
