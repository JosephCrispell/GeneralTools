#### Load libraries ####

library(ape) # Reading and phylogeny
library(rgdal) # Convert X and Y to lat longs and reading in shape files

#### Read in the sample data ####

# Set the path 
path <- "J:\\WGS_Monaghan\\"

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

# Read in the animal IDs
file <- paste0(path, "EartagsAndSettIDs_08-08-19.csv")
animalTags <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Read in the cattle herd IDs
file <- paste0(path, "Animal_HerdIDs_04-07-19.csv")
herdIDs <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Read in the FASTA file
fastaFile <- paste0(path, "sequences_Prox-10_30-07-2019.fasta")
nSites <- getNSitesInFASTA(fastaFile)

# Read in the coverage information
coverageFile <- paste0(path, "Quality/isolateCoverageSummary_DP-20_30-07-2019.txt")
coverage <- read.table(coverageFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Read in the test summaries table
file <- paste0(path, "AHCS_data/Joe_20190731_TestSummaryOnly.csv")
cattleTestInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE,
                             check.names=FALSE)

# Read in the badger capture information
file <- paste0(path, "Badger_data/TBL_captured_badgers_2019_07_19.txt")
badgerCaptureData <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Read in the sett capture event information
file <- paste0(path, "Badger_data/TBL_capture_events_2019_07_19.txt")
settCaptureEventData <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

#### Build the phylogeny ####

# Build a phylogeny using RAxML
tree <- runRAXML(fastaFile, date="30-07-19", path, alreadyRun=TRUE, outgroup="\\>Ref-1997")

# Remove Reference
tree <- drop.tip(tree, tree$tip.label[grepl(tree$tip.label, pattern=">Ref-1997")])

# Edit the tip labels
tree$tip.label <- editTipLabels(tree$tip.label)

# Get the tip information (species and sampling date)
tipInfo <- getTipInfo(tree$tip.label, animalTags, coverage, herdIDs)

# Convert branch lengths to SNPs
tree$edge.length <- tree$edge.length * nSites

# Add in slaughter and capture dates
tipInfo <- addSlaughterOrCaptureDates(tipInfo, cattleTestInfo, badgerCaptureData, settCaptureEventData)

#### Plot the phylogeny ####

# Open a PDF
pdf(paste0(path, "MbovisAnnotatedPhylogeny_", date, ".pdf"))

# Plot the phylogeny
plot.phylo(tree, show.tip.label=FALSE, edge.color=rgb(0,0,0, 0.5), edge.width=4)

# Add tips coloured by species
tipShapesAndColours <- list("BADGER"=c(rgb(1,0,0,1), 19), "COW"=c(rgb(0,0,1,1), 17), "NA"=c("grey", 15))
tiplabels(pch=getTipShapeOrColourBasedOnSpecies(tipInfo, tipShapesAndColours, which="shape"),
          col=getTipShapeOrColourBasedOnSpecies(tipInfo, tipShapesAndColours, which="colour"), cex=1.25)

# Add scale bar
addScaleBar(20)

# Add species legend
addSpeciesLegend(tipShapesAndColours)

# Close the pdf
dev.off()


#### Plot the temporal sampling ####

# Open a PDF
pdf(paste0(path, "TemporalSampling_", date, ".pdf"))

plotTemporalSamplingRange(tipInfo)

# Close the pdf
dev.off()

#### Read in the spatial data ####

# Get the badger sett locations
shapeFile <- paste0(path, "Badger_data/badgers.shp")
badgerShapeFile <- readOGR(dsn=shapeFile)
settLocations <- data.frame("SettID"=badgerShapeFile@data$SETT_NO,
                            "X"=badgerShapeFile@coords[, 1],
                            "Y"=badgerShapeFile@coords[, 2])

# Read in the cattle land parcels
shapeFile <- paste0(path, "Cattle_data/herds_mn_july_2019.shp")
cattleShapeFile <- readOGR(dsn = shapeFile)

# Extract the polygon coords
landParcelCoords <- getPolygonCoords(cattleShapeFile)

# Extract the sampld herd information
herdInfo <- cattleShapeFile@data

# Calculate the herd centroids
landParcelCentroids <- calculateLandParcelCentroids(landParcelCoords)


#### Plot the sampling locations ####

#### FUNCTIONS - spatial data ####

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
      centroidsX[fieldIndex] <- mean(coords[, "X"])
      centroidsY[fieldIndex] <- mean(coords[, "Y"])

    }

    # Store the calculated centroids
    herdCentroids[[herd]] <- list(
      "X"=mean(centroidsX),
      "Y"=mean(centroidsY),
      "CentroidXs"=centroidsX,
      "CentroidYs"=centroidsY)
  }
  
  return(herdCentroids)
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
    
    # Convert the ID to a herd ID
    id <- as.character(spatialDataFrame@data[as.numeric(id) + 1, "HERD_NO"])
    
    # Note the number of polygons in the current set
    nPolygons <- length(spatialDataFrame@polygons[[setIndex]]@Polygons)
    
    # Initialise a list to store the coordinates of each polygon in current set
    polygons <- list()
    
    # Examine each polygon in current set
    for(polygonIndex in seq_len(nPolygons)){
      
      # Get the coordinates for the current polygon
      coords <- spatialDataFrame@polygons[[setIndex]]@Polygons[[polygonIndex]]@coords
      colnames(coords) <- c("X", "Y")

      # Store the latitude and longitudes and X and Y coords
      polygons[[polygonIndex]] <- coords
    }
    
    # Store the polygon coordinates
    output[[id]] <- polygons
  }
  
  return(output)
}

#### FUNCTIONS - plotting temporal sampling ####

plotTemporalSamplingRange <- function(tipInfo){
  
  # Convert the date format
  tipInfo$DateAtSlaughterOrCapture <- as.Date(tipInfo$DateAtSlaughterOrCapture)
  
  # Calculate the range in the sampling dates
  dateRange <- range(tipInfo$DateAtSlaughterOrCapture, na.rm=TRUE)
  
  # Get and set the plotting margins
  currentMar <- par()$mar
  par(mar=c(4.1, 5.1, 4.1, 4.4))
  
  # Define the locations of the sampling ranges on the Y axis
  yLocations <- c(0.15, 0.35)
  speciesShapes <- c(19, 19)
  
  # Create an empty plot
  plot(x=NULL, y=NULL, ylim=c(0, 0.5), xlim=dateRange,
       bty="n", ylab="", yaxt="n", xlab="", xaxt="n",
       main="Number of samples available in time",
       cex.main=2)
  
  # Add an X axis
  at <- as.Date(c("2018-05-15", "2018-08-15", "2018-11-15", "2019-02-15", "2019-05-15"),  format="%Y-%m-%d")
  axis(side=1, at=at, labels=FALSE, xpd=TRUE)
  axis(side=1, at=at, labels=format(at, "%d-%m-%Y"), tick=FALSE, cex.axis=1,
       line=1, xpd=TRUE)
  
  # Add y axis
  species <- c("BADGER", "COW")
  axis(side=2, at=yLocations, labels=species, las=1, tick=FALSE, cex.axis=1.25,
       line=-0.5)
  
  # Examine the sampling of each species
  for(speciesIndex in seq_along(species)){
    
    # Get the sampling dates for the current species
    dates <- tipInfo[tipInfo$Species == species[speciesIndex], "DateAtSlaughterOrCapture"]
    
    # Remove NAs
    dates <- dates[is.na(dates) == FALSE]
    
    # Order dates
    dates <- sort(dates)
    
    # Plot the dates
    points(x=dates, y=rep(yLocations[speciesIndex], length(dates)),
           type="o", pch=speciesShapes[speciesIndex], col=rgb(0,0,0, 0.25), cex=5)
    
    # Note the number of samples available
    text(x=dateRange[2]-(0.001*(dateRange[2] - dateRange[1])),
         y=yLocations[speciesIndex], adj=0, xpd=TRUE,
         labels=paste0("     ", length(dates), "/",
                       length(which(tipInfo$Species == species[speciesIndex]))),
         cex=1.25)
  }
  
  # Reset the plotting margins
  par(mar=currentMar)
}

#### FUNCTIONS - linking to metadata ####

addSlaughterOrCaptureDates <- function(tipInfo, cattleTestInfo, badgerCaptureData, settCaptureEventData){
  
  # Add a date into the tip info table
  tipInfo$DateAtSlaughterOrCapture <- NA
  
  # Examine each row of tip information
  for(row in seq_len(nrow(tipInfo))){
    
    # Skip if no animal ID available
    if(is.na(tipInfo[row, "AnimalID"])){
      next
    }
    
    # Check if cow
    if(tipInfo[row, "Species"] == "COW"){
      
      # Try to find test information row for current cow
      testRow <- which(cattleTestInfo[, "Animal No."] == tipInfo[row, "AnimalID"])
      
      # Check if found
      if(length(testRow) == 0){
        warning("Not able to find testing information for cattle eartag: ", tipInfo[row, "AnimalID"])
      
      # Look at slaughter dates
      }else{
        
        # Get the unique slaughter dates for the current animal
        slaughterDates <- unique(as.Date(cattleTestInfo[testRow, "PM (Date of SL)"], format="%d/%m/%Y"))
        
        # Check if more than unique date available
        if(length(slaughterDates) > 1){
          warning("Multiple slaughter dates available for cow: ", tipInfo[row, "AnimalID"])
        }
        
        # Store a slaughter date
        tipInfo[row, "DateAtSlaughterOrCapture"] <- as.character(slaughterDates[1])
      }
    
    # Get the capture date for the badger
    }else{
      
      # Try and find badger in the capture event table
      badgersRow <- which(badgerCaptureData$BADGER_ID == tipInfo[row, "AnimalID"])
      
      # Skip badger if didn't find a record
      if(length(badgersRow) == 0){
        warning(paste0("Unable to find capture data for badger: ", tipInfo[row, "AnimalID"]))
        next
      }
      
      # Store the sett ID
      tipInfo[row, "HerdOrSettID"] <- badgerCaptureData[badgersRow, "SETT_ID"]
      
      # Find the sett capture event information
      eventRow <- which(settCaptureEventData$CAPTURE_BLOCK_EVENT == badgerCaptureData[badgersRow, "CAPTURE_BLOCK_EVENT"])
      
      # Skip if weren't able to find capture event row
      if(length(eventRow) == 0){
        warning(paste0("Unable to find information for sett capture event:: ", badgerCaptureData[badgersRow, "CAPTURE_BLOCK_EVENT"]))
        next
      }
      
      # Estimate the date of capture for the current badger
      tipInfo[row, "DateAtSlaughterOrCapture"] <- estimateCaptureDate(settCaptureEventData[eventRow, "DATE_COMMENCED"],
                                                                      settCaptureEventData[eventRow, "DATE_COMPLETED"])
    }
  }
}

estimateCaptureDate <- function(dateTrappingCommenced, dateTrappingCompleted){
  
  # Convert the date strings to date objects
  dateTrappingCommenced <- as.Date(strsplit(dateTrappingCommenced, split=" ")[[1]][1], format="%d-%b-%y")
  dateTrappingCompleted <- as.Date(strsplit(dateTrappingCompleted, split=" ")[[1]][1], format="%d-%b-%y")
  
  # Get the middle date
  middleDate <- dateTrappingCommenced + ((dateTrappingCompleted - dateTrappingCommenced)/2)
  
  return(as.character(middleDate))
}

getTipInfo <- function(tipLabels, animalTags, coverage, herdIDs){
  
  # Initialise a dataframe to store the tip information
  tipInfo <- data.frame("Tip"=tipLabels,"Aliquot"=NA, "AnimalID"=NA, "Species"=NA,
                        "Coverage"=NA, "HerdOrSettID"=NA, 
                        stringsAsFactors=FALSE)
  
  # Examine each of the tips
  for(index in seq_along(tipLabels)){
    
    # Build an aliquot code for the current isolate
    aliquot <- paste0("TB19-", paste(rep(0, 6-nchar(tipLabels[index])), collapse=""), tipLabels[index])
    tipInfo[index, "Aliquot"] <- aliquot
    
    # Find the row in the coverage information for the current tip
    coverageRow <- which(grepl(coverage$IsolateID, pattern=paste0("^", tipLabels[index])))
    
    # Check that row was found
    if(length(coverageRow) == 0){
      warning("Unable to find coverage information for: ", tipLabels[index],
              "\tAliquot: ", aliquot)
    }else if(length(coverageRow) > 1){
      warning("Multiple entries in coverage information for: ", tipLabels[index],
              "\tAliquot: ", aliquot)
    }else{
      # Store the coverage information
      tipInfo[index, "Coverage"] <- coverage[coverageRow, "PercentageCoverage"]
    }
    
    # Find the row in the sample information table for the current tip
    tagRow <- which(animalTags$Aliquot == aliquot)
    
    # Check that row was found
    if(length(tagRow) == 0){
      warning("Unable to find sampling information for: ", tipLabels[index],
              "\tAliquot: ", aliquot)
    }else if(length(tagRow) > 1){
      warning("Multiple entries in sampling information for: ", tipLabels[index],
              "\tAliquot: ", aliquot)
    }else{
      # Store the animal ID
      tipInfo[index, "AnimalID"] <- animalTags[tagRow, "Animal.ID"]
      
      # Store the species of the current tip
      tipInfo[index, "Species"] <- ifelse(grepl(tipInfo[index, "AnimalID"], pattern="^RR"),
                                          "BADGER", "COW")
    }
    
    # Skip badgers
    if(is.na(tipInfo[index, "Species"]) == FALSE && 
       tipInfo[index, "Species"] == "BADGER"){
      next
    }
    
    # Find the row in the herd IDs for the current tip (note we'll ignore badgers)
    herdRow <- which(herdIDs$Aliquot == aliquot)
    
    # Check that row was found
    if(length(herdRow) == 0){
      warning("Unable to find herd for: ", tipLabels[index],
              "\tAliquot: ", aliquot)
    }else if(length(herdRow) > 1){
      warning("Multiple entries in herd IDs table for: ", tipLabels[index],
              "\tAliquot: ", aliquot)
    }else{
      # Store the animal ID
      tipInfo[index, "HerdOrSettID"] <- herdIDs[herdRow, "Herd_ID"]
    }
  }
  
  return(tipInfo)
}

#### FUNCTIONS - plotting phylogeny ####

editTipLabels <- function(tipLabels){
  
  # Initialise a vector to store the new labels
  output <- c()
  
  # Examine each tip label
  for(label in tipLabels){
    
    # Remove the ">" prefix
    label <- substr(label, 2, nchar(label))
    
    # Split the label and retain first part
    label <- strsplit(label, split="_")[[1]][1]
    
    # Remove a trailing "p"
    label <- gsub("p", "", label)
    
    # Store the new label
    output[length(output) + 1] <- label
  }
  
  return(output)
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
  xSpace <- 0.2*xLength
  
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

addScaleBar <- function(scaleSize, cex=1){
  
  # Get the plotting region dimensions = x1, x2, y1, y2 
  # (coordinates of bottom left and top right corners)
  dimensions <- par("usr")
  xLength <- dimensions[2] - dimensions[1]
  yLength <- dimensions[4] - dimensions[3]
  
  # Add Scale bar
  xPad <- 0.4 * xLength
  
  points(x=c(dimensions[1] + xPad, dimensions[1] + xPad + scaleSize), 
         y=c(dimensions[3] + (0.01 * yLength), dimensions[3] + (0.01 * yLength)),
         type="l", lwd=3, xpd=TRUE)
  if(scaleSize == 1){
    text(x=dimensions[1] + xPad + (0.5*scaleSize), y=dimensions[3] - (0.02 * yLength),
         labels=paste0("~ ", scaleSize, " SNV"), cex=cex, xpd=TRUE)
  }else{
    text(x=dimensions[1] + xPad + (0.5*scaleSize), y=dimensions[3] - (0.02 * yLength),
         labels=paste0("~ ", scaleSize, " SNVs"), cex=cex, xpd=TRUE)
  }
}

#### FUNCTIONS - Phylogeny ####

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
        stop(paste("Option", which, "not recognised!"))
      }
      
      # If species available assign appropriate colour or shape
    }else{
      
      # Check if wanting shape or colour
      if(which == "shape"){
        output[row] <- as.numeric(tipShapesAndColours[[tipInfo[row, "Species"]]][2])
      }else if(which == "colour"){
        output[row] <- tipShapesAndColours[[tipInfo[row, "Species"]]][1]
      }else{
        stop(paste("Option", which, "not recognised!"))
      }
    }
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
