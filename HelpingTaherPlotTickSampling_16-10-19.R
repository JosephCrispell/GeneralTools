#### Preparation ####

# Load libraries
library(rgdal) # For reading in shape files
library(sp)
library(basicPlotteR)

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

#### Read in the tick count data BORRELIA ####

# Read in the tick counts
# Removed location columns as had characters interferring with reading in
tickCountsFileBorrelia <- paste0(path, "HelpingTaher/Taher_Borrelia_16-10-19.tsv")
tickCountsBorrelia <- read.table(tickCountsFileBorrelia, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Calculate prevalence
tickCountsBorrelia$Prevalence <- tickCountsBorrelia$No..of.qPCR.positives / tickCountsBorrelia$Ticks.screened
tickCountsBorrelia$Prevalence[is.nan(tickCountsBorrelia$Prevalence)] <- 0

#### Plot the sampling for BORRELIA ####

# Open a pdf
file <- paste0(substr(tickCountsFileBorrelia, 1, nchar(tickCountsFileBorrelia)-4), ".pdf")
pdf(file, width=14)

# Get and set the plotting margins
currentMar <- par()$mar
par(mar=c(0,0,2,0), mfrow=c(1,2))

# Plot the polygons from Ireland
plotPolygons(irelandPolygons, lwd=1.5, border=rgb(0.4,0.4,0.4, 0.1), col=rgb(0,0,1, 0.1),
             main="Sites sampled for presence of ticks from 2017 to 2019", label="a")

# Plot the locations that were sampled
points(x=tickCountsBorrelia$Long, y=tickCountsBorrelia$Lat, cex=1, pch=19, 
       col=ifelse(tickCountsBorrelia$No..of.collected.ticks > 0, rgb(1,0,0, 0.4), rgb(0,0,0, 0.4)))

# Add a legend
legend("right", legend=c("Present", "Absent"), pch=19, col=c(rgb(1,0,0, 0.4), rgb(0,0,0, 0.4)),
       bty="n")

# Plot the Ireland county polygons
plotPolygons(irelandPolygons, lwd=1.5, border=rgb(0.4,0.4,0.4, 0.1), col=rgb(0,0,1, 0.1),
             main="Prevelance of Borrelia at screened sites", label="b")

# Plot the tick counts
screenedSitesBorrelia <- tickCountsBorrelia[tickCountsBorrelia$Ticks.screened > 0, ]
addPoints(xCoords=screenedSitesBorrelia$Long, yCoords=screenedSitesBorrelia$Lat, pch=21,
          cex=(screenedSitesBorrelia$Ticks.screened/max(screenedSitesBorrelia$Ticks.screened))*4,
          bg=rgb(1,0,0, screenedSitesBorrelia$Prevalence),
          col="black", lty=2)

# Add a scale for number of ticks
legend("right", title="No. ticks screened", legend=c(100, 75, 50, 25, 10), pch=21, 
       col="black", bty="n", pt.cex=(c(100, 75, 50, 25, 10)/max(screenedSitesBorrelia$Ticks.screened))*4)

# Add a scale for prevalence
legend("bottomright", title="Prevalence         ", legend=paste0(c(50, 40, 30, 20, 10, 0), "%"), pch=21, 
       pt.bg=c(rgb(1,0,0, 0.5), rgb(1,0,0, 0.4), rgb(1,0,0, 0.3), 
            rgb(1,0,0, 0.2), rgb(1,0,0, 0.1), rgb(1,0,0, 0)),
       bty="n", pt.cex=2, col="black")

# Reset the plotting margins
par(mar=currentMar, mfrow=c(1,1))

# Close the pdf
dev.off()

#### Read in the tick counts for other pathogens ####

# Read in the tick counts
# Removed location columns as had characters interferring with reading in
tickCountsFileOtherPathogens <- paste0(path, "HelpingTaher/Taher_OtherPathogens_07-11-19.tsv")
tickCountsOtherPathogens <- read.table(tickCountsFileOtherPathogens, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Calculate prevalence for the phagocytophilum
tickCountsOtherPathogens$PrevalencePhagocytophilum <- tickCountsOtherPathogens$No..of.A..phagocytophilum.positives / 
  tickCountsOtherPathogens$Ticks.screened
tickCountsOtherPathogens$PrevalencePhagocytophilum[is.nan(tickCountsOtherPathogens$PrevalencePhagocytophilum)] <- 0

# Calculate prevalence for the phagocytophilum
tickCountsOtherPathogens$PrevalenceDivergens <- tickCountsOtherPathogens$No..of.B..divergens.positives / 
  tickCountsOtherPathogens$Ticks.screened
tickCountsOtherPathogens$PrevalenceDivergens[is.nan(tickCountsOtherPathogens$PrevalenceDivergens)] <- 0


#### Plot the sampling for Other Pathogens ####

# Open a pdf
file <- paste0(substr(tickCountsFileOtherPathogens, 1, nchar(tickCountsFileOtherPathogens)-4), ".pdf")
pdf(file, width=14, height=14)

# Get and set the plotting margins
currentMar <- par()$mar
par(mar=c(0,0,2,0), mfrow=c(2,2))

# Plot the polygons from Ireland
plotPolygons(irelandPolygons, lwd=1.5, border=rgb(0.4,0.4,0.4, 0.1), col=rgb(0,0,1, 0.1),
             main="Sites sampled for presence of ticks from 2017 to 2019", label="a")

# Plot the locations that were sampled
points(x=tickCountsOtherPathogens$Long, y=tickCountsOtherPathogens$Lat, cex=1, pch=19, 
       col=ifelse(tickCountsOtherPathogens$No..of.collected.ticks > 0, rgb(1,0,0, 0.4), rgb(0,0,0, 0.4)))

# Add a legend
legend("right", legend=c("Present", "Absent"), pch=19, col=c(rgb(1,0,0, 0.4), rgb(0,0,0, 0.4)),
       bty="n")

# Plot the Ireland county polygons
plotPolygons(irelandPolygons, lwd=1.5, border=rgb(0.4,0.4,0.4, 0.1), col=rgb(0,0,1, 0.1),
             main="Prevelance of A. phagocytophilum at screened sites", label="b")

# Plot the tick counts
screenedSitesOtherPathogens <- tickCountsOtherPathogens[tickCountsOtherPathogens$Ticks.screened > 0, ]
addPoints(xCoords=screenedSitesOtherPathogens$Long, yCoords=screenedSitesOtherPathogens$Lat, pch=21,
          cex=(screenedSitesOtherPathogens$Ticks.screened/max(screenedSitesOtherPathogens$Ticks.screened))*4,
          bg=rgb(1,0,0, screenedSitesOtherPathogens$PrevalencePhagocytophilum),
          col="black", lty=2)

# Add a scale for number of ticks
legend("right", title="No. ticks screened", legend=c(100, 75, 50, 25, 10), pch=21, 
       col="black", bty="n", pt.cex=(c(100, 75, 50, 25, 10)/max(screenedSitesOtherPathogens$Ticks.screened))*4)

# Add a scale for prevalence
legend("bottomright", title="Prevalence         ", legend=paste0(c(50, 40, 30, 20, 10, 0), "%"), pch=21, 
       pt.bg=c(rgb(1,0,0, 0.5), rgb(1,0,0, 0.4), rgb(1,0,0, 0.3), 
               rgb(1,0,0, 0.2), rgb(1,0,0, 0.1), rgb(1,0,0, 0)),
       bty="n", pt.cex=2, col="black")

# Plot the Ireland county polygons
plotPolygons(irelandPolygons, lwd=1.5, border=rgb(0.4,0.4,0.4, 0.1), col=rgb(0,0,1, 0.1),
             main="Prevelance of B. divergens at screened sites", label="c")

# Plot the tick counts
screenedSitesOtherPathogens <- tickCountsOtherPathogens[tickCountsOtherPathogens$Ticks.screened > 0, ]
addPoints(xCoords=screenedSitesOtherPathogens$Long, yCoords=screenedSitesOtherPathogens$Lat, pch=21,
          cex=(screenedSitesOtherPathogens$Ticks.screened/max(screenedSitesOtherPathogens$Ticks.screened))*4,
          bg=rgb(1,0,0, screenedSitesOtherPathogens$PrevalenceDivergens),
          col="black", lty=2)

# Add a scale for number of ticks
legend("right", title="No. ticks screened", legend=c(100, 75, 50, 25, 10), pch=21, 
       col="black", bty="n", pt.cex=(c(100, 75, 50, 25, 10)/max(screenedSitesOtherPathogens$Ticks.screened))*4)

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
    axisLimits <- par()$usr
    xLength <- axisLimits[2] - axisLimits[1]
    mtext(label, side=3, line=-1, at=axisLimits[1]+(0.10*xLength), cex=2.5)
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
