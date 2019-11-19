#### Preparation ####

# Load libraries
library(ape) # Reading and phylogeny
library(geiger) # For the tips function in used in defineBranchColoursOfClade
library(basicPlotteR) # Because I made it! :-)

# Set the path 
path <- file.path("~", "Desktop", "Research", "EdgeArea_UK")

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

#### Read in sequencing data ####

# Read in the sample metadata
sampleInfoFile <- file.path(path, "SampleInformation_28-05-19.csv")
sampleInfo <- read.table(sampleInfoFile, header=TRUE, sep=",", stringsAsFactors=FALSE)
sampleInfo <- sampleInfo[sampleInfo$Sample.Ref != "", ]
rownames(sampleInfo) <- sampleInfo$Sample.Ref

# Read in the FASTA file
fastaFile <- file.path(path, "vcfFiles", "sequences_Prox-10_31-10-2019.fasta")
sequences <- read.dna(fastaFile, format="fasta", as.character=TRUE)
nSites <- ncol(sequences)

# Plot the FASTA file
fastaPlotFile <- paste0(substr(fastaFile, 0, nchar(fastaFile) - 6), ".pdf")
plotFASTA(sequences, pdfFileName=fastaPlotFile, pdfHeight=14, pdfWidth=21, lineForSequenceNames=-3)

# Calculate sequencing quality in the FASTA file
# Note average proportion Ns may be higher because of the extreme coverage in one genome - reaching areas with no coverage in other genomes
propNs <- data.frame("ID"=rownames(sequences), "ProportionNs"=sapply(1:nrow(sequences), calculatePropNsInSequence, sequences),
                     stringsAsFactors=FALSE)
rownames(propNs) <- parseTipLabels(propNs$ID)

#### Build the phylogeny ####

# Build a phylogeny using RAxML
tree <- runRAXML(fastaFile, date="31-10-19", file.path(path, "vcfFiles", ""), outgroup="\\>Ref-1997")

# Remove Reference
tree <- drop.tip(tree, tree$tip.label[grepl(tree$tip.label, pattern=">Ref-1997")])

# Convert the branch lengths to SNPs
if(mean(tree$edge.length) < 1){
  tree$edge.length <- tree$edge.length * nSites
}

# Parse the tip labels
tree$tip.label <- parseTipLabels(tree$tip.label)

#### Get the tip information ####

# Get the tip sampling and quality information
tipInfo <- getTipInfo(tree$tip.label, sampleInfo, propNs)

# Count number of tips that we have no data for
warning(length(which(is.na(tipInfo$Sample.Ref))), " sequence labels not found in sample information!")

#### Plot the phylogeny ####

# Extract the densely sample clade
node <- 121
clade <- extract.clade(tree, node=node)

# Open an output pdf
pdf(file.path(path, paste0("SummaryPlots_EDGE_", date, ".pdf")))

# Plot the full phylogeny and highlight the densely sampled clade
branchColours <- defineBranchColoursOfClade(tree, nodeDefiningClade=node, colour="red", defaultColour="black")
plot.phylo(tree, show.tip.label=FALSE, edge.width=2, edge.color=branchColours)
addSNPScale(position="bottom", size=50, lineWidth=4, cex=1.5)

# Plot the densely sample clade
plot.phylo(clade, show.tip.label=FALSE, edge.color=rgb(0.1,0.1,0.1, 1), edge.width=2)
addSNPScale(position="bottom", size=1, lineWidth=4, cex=1.5)

# Plot the spatial sampling
plotSamplingLocations(tipInfo$Mapx, tipInfo$Mapy, bty="n", xaxt="n", yaxt="n", xlab="", ylab="", pch=19,
                      asp=1, main="Sampling locations", 
                      col=ifelse(tipInfo$Sample.Ref %in% clade$tip.label, rgb(1,0,0, 0.25), rgb(0,0,0, 0.25)))
legend("right", legend=c("all", "in clade"), text.col=c(rgb(0.1,0.1,0.1, 0.75), rgb(1,0,0, 0.75)), bty="n")

# Plot the temporal sampling
plotMultipleHistograms(distributions=list(tipInfo$Year, tipInfo[tipInfo$Sample.Ref %in% clade$tip.label, "Year"]),
                       colours=c(rgb(0.1,0.1,0.1, 0.5), rgb(1,0,0, 0.5)), nBins=25, las=1, xaxt="n", xlab="Year",
                       main="Sampling dates")
yearRange <- range(tipInfo$Year, na.rm=TRUE)
axis(side=1, at=seq(from=yearRange[1], to=yearRange[2]+2, by=2), xpd=TRUE)
legend("top", legend=c("all", "in clade"), text.col=c(rgb(0.1,0.1,0.1, 0.75), rgb(1,0,0, 0.75)), bty="n")

# Close the output pdf
dev.off()

#### FUNCTIONS ####

defineBranchColoursOfClade <- function(tree, nodeDefiningClade,
                                       colour, defaultColour){
  branchColours <- rep(defaultColour, dim(tree$edge)[1])
  clade <- tips(tree, node=nodeDefiningClade)
  branchesInClades <- which.edge(tree, clade)
  branchColours[branchesInClades] <- colour
  
  return(branchColours)
}

plotSamplingLocations <- function(xCoords, yCoords, scaleSize=10000, xPad=0.75, yPad=0.1, ...){
  
  # Plot the sampling locations
  plot(x=xCoords, y=yCoords, ...)
  
  # Note the size of the axes
  axisLimits <- par()$usr
  xAxisSize <- axisLimits[2] - axisLimits[1]
  yAxisSize <- axisLimits[4] - axisLimits[3]
  
  # Note the start and end of scale box
  xStart <- axisLimits[1] + (xAxisSize*xPad)
  xEnd <- xStart + scaleSize
  yStart <- axisLimits[3] + (yAxisSize*yPad)
  yEnd <- yStart + scaleSize
  
  # Add a scale box
  rect(xleft=xStart, ybottom=yStart, xright=xEnd, ytop=yEnd)
  
  # Add label for scale box
  label <- bquote(paste(.((scaleSize/1000)*(scaleSize/1000)), "km")^2)
  text(x=xStart+(0.5*scaleSize), y=yStart, labels=label, pos=1)
}

getTipInfo <- function(tipLabels, sampleInfo, propNs){
  
  # Set the row names of the sample info table
  rownames(sampleInfo) <- sampleInfo$Sample.Ref
  
  # Create the tipInfo table
  tipInfo <- sampleInfo[tipLabels, c("Sample.Ref", "Mapx", "Mapy", "Eartag", "CPHH", "Breakdown.ID", "County", "Year", "TYPING.RESULT")]
  rownames(tipInfo) <- tipLabels
  
  # Add the proportion Ns to tipInfo table
  rownames(propNs) <- parseTipLabels(propNs$ID)
  tipInfo$ProportionNs <- propNs[tipLabels, "ProportionNs"]
  
  return(tipInfo)
}

parseTipLabels <- function(tipLabels){
  
  # Check if already parsed
  if(grepl(tipLabels[1], pattern="_") == FALSE){
    return(tipLabels)
  }
  
  # Initialise a vector to store the parsed labels
  parsedLabels <- c()
  
  # Examine eahc tip label
  for(index in seq_along(tipLabels)){
    
    # Remove the ">" from the start of the tip label
    label <- tipLabels[index]
    if(grepl(tipLabels[1], pattern="^>")){
      label <- substr(tipLabels[index], 2, nchar(tipLabels[index]))
    }
    
    # Split the label and select first part
    label <- strsplit(label, split="_")[[1]][1]
    
    # If it doesn't start with "AF" change "-" to "/"
    if(grepl(label, pattern="^AF") == FALSE){
      
      label <- gsub(pattern="-", replacement="/", label)
    }
    
    # Store the parsed label
    parsedLabels[index] <- label
  }
  
  return(parsedLabels)
}

calculatePropNsInSequence <- function(sequenceIndex, sequences){
  
  # Count the number of Ns in the sequence
  nucleotideCounts <- table(sequences[sequenceIndex, ])
  names(nucleotideCounts) <- toupper(names(nucleotideCounts))
  
  # Note the number of Ns
  numberMissing <- 0
  if("N" %in% names(nucleotideCounts)){
    numberMissing <- nucleotideCounts[["N"]]
  }
  
  return(numberMissing / ncol(sequences))
}

runRAXML <- function(fastaFile, date, path, nBootstraps=100, nThreads=10, outgroup=NULL, model="GTRCAT"){
  
  # Note the directory for the RAxML output files
  directory <- paste(path, "RAxML_", date, sep="")
  
  # Check if analyses already run
  alreadyRun <- dir.exists(directory)
  
  # If not already run, create output directory for RAxML
  if(alreadyRun == FALSE){
    suppressWarnings(dir.create(directory))
  }
  
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
