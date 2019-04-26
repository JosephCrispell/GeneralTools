#### Load libraries ####

library(ape)
library(phangorn)
library(phytools)
library(geiger)

#### Load data ####

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

# Create a path variable
path <- "/home/josephcrispell/Desktop/Research/RepublicOfIreland/Mbovis/"

# Read in table that links original sequence ID to aliquot IDs
file <- paste0(path, "Mbovis_CattleSamplingInfo_17-07-18.tsv")
linkTable <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep="\t")

# Read in the isolate metadata
file <- paste0(path, "IsolateSpeciesAndYear_26-04-19.csv")
metadata <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep=",")

# Read in the FASTA file
fastaFile <- paste(path, "vcfFiles/sequences_Prox-10_19-03-2019.fasta", sep="")
nSites <- getNSitesInFASTA(fastaFile)

#### Build phylogeny ####

# Build a phylogeny using RAxML
tree <- runRAXML(fastaFile, date="26-04-19", path, alreadyRun=TRUE, outgroup="\\>Ref-1997")

# Remove NI isolates and Reference
tree <- drop.tip(mlTree, mlTree$tip.label[grepl(mlTree$tip.label, pattern=">Ref-1997|>182-MBovis|>161-MBovis")])

# Edit the tip labels
tree$tip.label <- editTipLabels(tree$tip.label)

# Get the tip information (species and sampling date)
tipInfo <- getTipInfo(tree$tip.label, metadata, linkTable)

# Convert branch lengths to SNPs
tree$edge.length <- tree$edge.length * nSites

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
tipShapesAndColours <- list("Badger"=c("red", 19), "Cow"=c("blue", 17), "Deer"=c("black", 15), "NA"=c("grey", 18))
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

#### FUNCTIONS ####

addSpeciesLegend <- function(tipShapesAndColours){
  
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
           col=tipShapesAndColours[[species[i]]][1], xpd=TRUE)
    
    # Add a label
    text(x=xStart+(i*xSpace), y=yPos, labels=species[i], pos=4, xpd=TRUE)
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
    text(x=dimensions[1] + xPad + (0.5*scaleSize), y=dimensions[3] - (0.03 * yLength),
         labels=paste0("~ ", scaleSize, " SNP"), cex=cex, xpd=TRUE)
  }else{
    text(x=dimensions[1] + xPad + (0.5*scaleSize), y=dimensions[3] - (0.03 * yLength),
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

getTipInfo <- function(tipLabels, metadata, linkTable){
  
  # Initialise a dataframe to store the tip information
  tipInfo <- data.frame(ID=tipLabels, Species=NA, Date=NA, stringsAsFactors=FALSE)
  
  # Examine each of the tips
  for(index in seq_along(tipLabels)){
    
    # Initialise variables to store the tip's information
    aliquotCode <- NA
    species <- NA
    date <- NA
    
    # Check if the current tip is associated with the original dataset
    if(grepl(tipLabels[index], pattern="-MBovis")){
      
      # Get the sequence number from the curren tip label
      sequenceNumber <- strsplit(tipLabels[index], split="-")[[1]][1]
      
      # Find the row in the link table
      row <- which(linkTable$Isolate.Code == sequenceNumber)
      
      # Get the current tips aliquot code
      if(length(row) != 0){
        aliquotCode <- linkTable[row, "Aliquot"]
      }else if(sequenceNumber %in% c(14, 23)){
        species <- "Deer"
      }else{
        cat(paste("Error for old batch. Couldn't find sequence number: ", sequenceNumber, "\n\n"))
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
      if(species %in% c("Heifer", "Steer", "Calf", "Bull")){
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
    
    # Store the current tips information
    tipInfo[index, "Species"] <- species
    tipInfo[index, "Date"] <- date
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
