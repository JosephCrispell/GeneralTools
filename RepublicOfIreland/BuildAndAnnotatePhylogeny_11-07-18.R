#### Load packages ####
library(ape)
library(phangorn)
library(phytools)
library(geiger)

#### Read FASTA ####

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

# Create a path variable
path <- "/home/josephcrispell/Desktop/Research/RepublicOfIreland/Mbovis/"

# Read in the FASTA file
fastaFile <- paste(path, "vcfFiles/sequences_Prox-10_09-01-2019.fasta", sep="")
sequencesDNAbin <- read.dna(fastaFile,
                            format = "fasta", skip=1) # skip first line - I added this line into FASTA: nSequences length

#### Read isolate data ####

# Read in the table
file <- paste(path, "Mbovis_CattleSamplingInfo_17-07-18.tsv", sep="")
sampleInfo <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep="\t")

# Format the date columns
sampleInfo$Born <- as.Date(sampleInfo$Born, format="%d/%m/%Y")
sampleInfo$Died <- as.Date(sampleInfo$Died, format="%d/%m/%Y")
sampleInfo$Test.Pos.Date[sampleInfo$Test.Pos.Date == "Abbatoir check"] <- NA
sampleInfo$Test.Pos.Date <- as.Date(sampleInfo$Test.Pos.Date, format="%d/%m/%Y")
sampleInfo$Moved.in[sampleInfo$Moved.in == ""] <- NA
sampleInfo$Moved.in <- as.Date(sampleInfo$Moved.in, format="%d/%m/%Y")

#### Build tree ####

# Build a phylogeny using RAxML
mlTree <- runRAXML(fastaFile, date, path, alreadyRun=TRUE)

# Remove NI isolates and Reference
mlTree <- drop.tip(mlTree, mlTree$tip.label[grepl(mlTree$tip.label, pattern=">Ref-1997|>182-MBovis|>161-MBovis")])

# Parse the isolate IDs
mlTree$tip.label <- parseIDs(mlTree$tip.label)

# Change the branch lengths to SNPs
mlTree$edge.length <- mlTree$edge.length * getNSitesInFASTA(fastaFile)

#### Plot the tree with sampling info ####

# Open a pdf
pdf(paste(path, "MBovisAnnotatedPhylogeny_09-09-19.pdf", sep=""))

# Set the margins
par(mai=c(0.75,0,0,3))

# Plot the phylogeny
plot.phylo(mlTree, edge.color="grey", edge.width=2, label.offset=0.5, 
           tip.color=ifelse(mlTree$tip.label %in% c("49", "35"), rgb(0,0,0, 0.5), "black"))

# Bootstrap values
#bootstrapValues <- as.numeric(mlTree$node.label)/100
#bootstrapValues[is.na(bootstrapValues)] <- 0
#nodelabels(pch=20, frame="none", col=rgb(0,0,0, bootstrapValues))

# Add tip shapes
tiplabels(pch=19, col=ifelse(mlTree$tip.label %in% c("14", "23"), "pink", "black"))

# Get the tip coordinates
tipCoordinates <- getTipCoordinates(mlTree$tip.label)

# Get the axis limits
axisLimits <- par("usr")

# Add a scale bar
addScaleBar(length=2, axisLimits)

# Add vertical bars indicating the sampling area
addSamplingAreaBars(tipCoordinates, axisLimits)

# Plot the animal lifespans
plotLifespans(tipCoordinates, axisLimits, sampleInfo, yearAxisGap=4,
              tickLabelCex=0.75)

# Add tip colour legend
legend(x=axisLimits[2] - (0.3 * (axisLimits[2] - axisLimits[1])),
       y=axisLimits[3] + (0.1 * (axisLimits[4] - axisLimits[3])),
       legend=c("Same Herd", "Deer"), pch=19, text.col=c(rgb(0,0,0, 0.5), "black"), col=c("white", "pink"),
       cex=0.75)

dev.off()

#############
# FUNCTIONS #
#############

runRAXML <- function(fastaFile, date, path, nBootstraps=100, nThreads=6, alreadyRun=FALSE){
  
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
    model <- "GTRCAT" # No rate heterogenity
    seeds <- sample(1:100000000, size=2, replace=FALSE) # For parsimony tree and boostrapping
    
    command <- paste("raxmlHPC", 
                     " -f a", # Algorithm: Rapid boostrap inference
                     " -N ", nBootstraps,
                     " -T ", nThreads,
                     " -m ", model, " -V", # -V means no rate heterogenity
                     " -p ", seeds[1], " -x ", seeds[2], # Parsimony and boostrapping seeds
                     " -n ", analysisName,
                     " -s ", fastaFile, sep="")
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

addSamplingAreaBars <- function(tipCoordinates, axisLimits, propPadLeft=0.2){
  
  # Calculate the space available for the area bar
  lengthOfXAxis <- (axisLimits[2] - axisLimits[1]) / par("plt")[2]
  spaceLeft <- lengthOfXAxis - axisLimits[2]
  spaceAvailable <- propPadLeft * spaceLeft
  
  # Calculate the xPosition of the middle of the space available
  xPosition <- axisLimits[2]# + (0.5 * spaceAvailable)
  
  # Assign a colour to each area
  areaColours <- assignColoursToAreas(sampleInfo$Area, c("red", "blue", "cyan", "green", "brown", "darkorchid4", "magenta"))
  
  # Set the ends of the lines to be square
  currentLend <- par("lend")
  par(lend="square")
  
  # Examine each isolate
  for(row in 1:nrow(sampleInfo)){
    
    # Get the Y position of the current sample's tip
    yPosition <- tipCoordinates[[as.character(sampleInfo[row, "Isolate.Code"])]][2]
    
    # Add a vertical line for the current tip
    points(x=c(xPosition, xPosition), y=c(yPosition - 0.25, yPosition + 0.25),
           lwd=8, col=areaColours[[sampleInfo[row, "Area"]]], xpd=TRUE, type="l")
  }
  
  # Reset line ends
  par(lend=currentLend)
  
  # Calculate space available below phylogeny
  lengthOfYAxis <- (axisLimits[4] - axisLimits[3]) / (1 - par("plt")[3])
  spaceAtBottom <- lengthOfYAxis - (axisLimits[4] - axisLimits[3])
  
  # Add a label
  text(x=xPosition, y=axisLimits[3] - (0.65 * spaceAtBottom), labels="Area", xpd=TRUE)
}

assignColoursToAreas <- function(areas, colours){
  
  # Get the unique areas
  uniqueAreas <- unique(areas)
  
  # Initialise a list to store each area's colour
  areaColours <- list()
  
  # Assign a colour to each area
  for(i in 1:length(uniqueAreas)){
    
    if(i < length(colours)){
      areaColours[[uniqueAreas[i]]] <- colours[i]
    }else{
      cat("ERROR!! More areas than there are colours! Will begin recycling them\n")
      areaColours[[uniqueAreas[i]]] <- colours[i %% length(colours)]
    }
  }
  
  return(areaColours)
}

plotLifespans <- function(tipCoordinates, axisLimits, sampleInfo, yearAxisGap, propPadRight=0.1, propPadLeft=0.2,
                          tickLabelCex){
  
  # Calculate the amount of space left to right of tree
  # par("plt") gives you proportion of device width taken up by plot: https://www.rstudio.com/wp-content/uploads/2016/10/how-big-is-your-graph.pdf
  lengthOfXAxis <- (axisLimits[2] - axisLimits[1]) / par("plt")[2]
  spaceLeft <- lengthOfXAxis - axisLimits[2]
  
  # Calculate the range of the dates
  dateRange <- range(c(sampleInfo$Born, sampleInfo$Died, sampleInfo$Test.Pos.Date), na.rm=TRUE)

  # Note the amount of space available for a day
  proportionToBeUsed <- 1 - propPadRight - propPadLeft
  daySpace <- (proportionToBeUsed * spaceLeft) / as.numeric(dateRange[2] - dateRange[1])
  
  # Calculate sizes of left and right pads
  leftPadSize <- propPadLeft * spaceLeft
  rightPadSize <- propPadRight * spaceLeft
  
  # Set the ends of the lines to be square
  currentLend <- par("lend")
  par(lend="square")
  
  # Plot the lifespan of sample cow
  for(row in 1:nrow(sampleInfo)){
    
    # Get the Y position of the current sample's tip
    yPosition <- tipCoordinates[[as.character(sampleInfo[row, "Isolate.Code"])]][2]
    
    # Calculate the X start position
    xStart <- axisLimits[2] + leftPadSize + (as.numeric(sampleInfo[row, "Born"] - dateRange[1]) * daySpace)
    xEnd <- xStart + (as.numeric(sampleInfo[row, "Died"] - sampleInfo[row, "Born"]) * daySpace)
    
    # Plot the current animals lifespan
    points(x=c(xStart, xEnd), y=c(yPosition, yPosition), type="l", lwd=4, xpd=TRUE)
    
    # Add test positive date
    if(is.na(sampleInfo[row, "Test.Pos.Date"]) == FALSE){
      xTest <- axisLimits[2] + leftPadSize + (as.numeric(sampleInfo[row, "Test.Pos.Date"] - dateRange[1]) * daySpace)
      points(x=xTest, y=yPosition, col="red", pch=18, xpd=TRUE)
    }else{
      points(x=xEnd, y=yPosition, col="red", pch=18, xpd=TRUE)
    }
    
    # Add grey line up until move IN date
    if(is.na(sampleInfo[row, "Moved.in"]) == FALSE){
      xEnd <- xStart + (as.numeric(sampleInfo[row, "Moved.in"] - sampleInfo[row, "Born"]) * daySpace)
      points(x=c(xStart, xEnd), y=c(yPosition, yPosition), type="l", lwd=4, xpd=TRUE, col="grey")
    }
  }
  
  # Reset the ends of the lines
  par(lend=currentLend)
  
  # Calculate space available below lifespan plot
  lengthOfYAxis <- (axisLimits[4] - axisLimits[3]) / (1 - par("plt")[3])
  spaceAtBottom <- lengthOfYAxis - (axisLimits[4] - axisLimits[3])
  
  # Define start year dates within date range
  startOfYears <- defineYearStartsWithinRange(dateRange, yearGap=yearAxisGap)
  
  # Add year label for each of the start of year dates
  for(i in 1:length(startOfYears)){
    
    # Calculate the position for the current start of year date
    xPosition <- axisLimits[2] + leftPadSize + (as.numeric(startOfYears[i] - dateRange[1]) * daySpace)
    yPosition <- axisLimits[3] - (0.3 * spaceAtBottom)
    
    # Get the year from the current date
    year <- strsplit(as.character(startOfYears[i]), split="-")[[1]][1]
    
    # Add the year label
    text(x=xPosition, y=yPosition, labels=year, xpd=TRUE, cex=tickLabelCex)
    
    # Add tick
    points(x=c(xPosition, xPosition), y=c(axisLimits[3] - (0.1*spaceAtBottom), axisLimits[3] - (0.2*spaceAtBottom)), type="l", xpd=TRUE)
  }
  
  # Add axis bar
  start <- xPosition <- axisLimits[2] + leftPadSize + (as.numeric(startOfYears[1] - dateRange[1]) * daySpace)
  end <- xPosition <- axisLimits[2] + leftPadSize + (as.numeric(startOfYears[length(startOfYears)] - dateRange[1]) * daySpace)
  points(x=c(start, end),
         y=c(axisLimits[3] - (0.1*spaceAtBottom), axisLimits[3] - (0.1*spaceAtBottom)),
         type="l", xpd=TRUE)
  
  # Add axis label
  text(x=axisLimits[2] + leftPadSize + (0.5 * proportionToBeUsed * spaceLeft), y=axisLimits[3] - (0.65*spaceAtBottom), labels="Years", xpd=TRUE)
}

defineYearStartsWithinRange <- function(dateRange, yearGap){
  
  # Get the start year
  startYear <- as.numeric(strsplit(as.character(dateRange[1]), split="-")[[1]][1])
  
  # Get the end year
  endYear <- as.numeric(strsplit(as.character(dateRange[2]), split="-")[[1]][1])
  
  # Define start of year dates within this range
  startOfYears <- c()
  for(year in seq((startYear+1),(endYear-1), yearGap)){
    
    startOfYears[length(startOfYears) + 1] <- paste("01/01/", year, sep="")
  }
  startOfYears <- as.Date(startOfYears, format="%d/%m/%Y")
  
  return(startOfYears)
}

addScaleBar <- function(length, axisLimits){
  
  # Calculate space available below phylogeny
  lengthOfYAxis <- (axisLimits[4] - axisLimits[3]) / (1 - par("plt")[3])
  spaceAtBottom <- lengthOfYAxis - (axisLimits[4] - axisLimits[3])
  
  # Calculate the position for the scale bar
  xHalfway <- 0.5 * (axisLimits[2] - axisLimits[1])
  yPosition <- axisLimits[3] - (0.7 * spaceAtBottom)
  
  # Add the scale bar
  lines(x=c(xHalfway - (0.5*length), xHalfway + (0.5*length)), y=c(yPosition, yPosition), lwd=2, xpd=TRUE)
  
  # Add scale bar title
  if(length == 1){
    text(x=xHalfway, y=axisLimits[3] - (0.5 * spaceAtBottom), 
         labels=paste(length, "SNP"), xpd=TRUE)
  }else{
    text(x=xHalfway, y=axisLimits[3] - (0.5 * spaceAtBottom), 
         labels=paste(length, "SNPs"), xpd=TRUE)
  }
  
}

getTipCoordinates <- function(tipLabels){
  
  # Get all the information about the last phylogenetic tree plotted
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  
  # Create an empty list to store the coordinates
  tips <- list()
  
  # Examine each of the tip labels - order must match tree$tip.labels of plotted tree
  for(i in 1:length(tipLabels)){
    
    # Get and store the coordinates for the current tip label
    # Note that you are converting them to actual coordinates within the plotting window
    tips[[tipLabels[i]]] <- c(lastPP$xx[i], lastPP$yy[i])
  }
  
  return(tips)
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

parseIDs <- function(tipLabels){
  
  # Get the sample ID from the first column ("-" delmited), and remove ">" from start
  newLabels <- c()
  for(i in 1:length(tipLabels)){
    
    newLabels[i] <- strsplit(tipLabels[i], split="-")[[1]][1]
    newLabels[i] <- substr(newLabels[i], 2, nchar(newLabels[i]))
  }
  
  return(newLabels)
}
