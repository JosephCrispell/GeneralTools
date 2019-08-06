#### Load libraries ####

#### Read in the sample data ####

# Set the path 
#path <- "/home/josephcrispell/Desktop/Research/RepublicOfIreland/Mbovis/Monaghan/"
path <- "J:\\WGS_Monaghan\\"

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

# Read in the sample information
file <- paste0(path, "SampleInformation_04-07-19.csv")
sampleInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Read in the FASTA file
fastaFile <- paste(path, "Fastqs_29-07-19/vcfFiles/sequences_Prox-10_30-07-2019.fasta", sep="")
nSites <- getNSitesInFASTA(fastaFile)

# Read in the coverage information
coverageFile <- paste0(path, "Fastqs_29-07-19/vcfFiles/isolateCoverageSummary_DP-20_30-07-2019.txt")
coverage <- read.table(coverageFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)

#### Build the phylogeny ####

# Build a phylogeny using RAxML
tree <- runRAXML(fastaFile, date="30-07-19", path, alreadyRun=TRUE, outgroup="\\>Ref-1997")

# Remove Reference
tree <- drop.tip(tree, tree$tip.label[grepl(tree$tip.label, pattern=">Ref-1997")])

# Edit the tip labels
tree$tip.label <- editTipLabels(tree$tip.label)

# Get the tip information (species and sampling date)
tipInfo <- getTipInfo(tree$tip.label, sampleInfo, coverage)

# Convert branch lengths to SNPs
tree$edge.length <- tree$edge.length * nSites


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

getTipInfo <- function(tipLabels, sampleInfo, coverage){
  
  # Initialise a dataframe to store the tip information
  tipInfo <- data.frame(ID=tipLabels, Species=NA, Coverage=NA, stringsAsFactors=FALSE)
  
  # Examine each of the tips
  for(index in seq_along(tipLabels)){
    
    # Build an aliquot code for the current isolate
    aliquot <- paste0("TB19-", paste(rep(0, 6-nchar(tipLabels[index])), collapse=""), tipLabels[index])
    
    # Find the row in the sample information table for the current tip
    sampleInfoRow <- which(sampleInfo$Aliquot == aliquot)
    
    # Check that row was found
    if(length(sampleInfoRow) == 0){
      warning("Unable to find sampling information for: ", tipLabels[index], "\tAliquot: ", aliquot)
      next
    }else if(length(sampleInfoRow) > 1){
      warning("Multiple entries in sampling information for: ", tipLabels[index], "\tAliquot: ", aliquot)
      next
    }
    
    # Store the species of the current tip
    tipInfo[index, "Species"] <- sampleInfo[sampleInfoRow, "Species"]
    
    # Find the row in the coverage information for the current tip
    coverageRow <- which(grepl(coverage$IsolateID, pattern=paste0("^", tipLabels[index])))
    
    # Check that row was found
    if(length(coverageRow) == 0){
      warning("Unable to find coverage information for: ", tipLabels[index], "\tAliquot: ", aliquot)
      next
    }else if(length(coverageRow) > 1){
      warning("Multiple entries in coverage information for: ", tipLabels[index], "\tAliquot: ", aliquot)
      next
    }
    
    # Store the coverage information
    tipInfo[index, "Coverage"] <- coverage[coverageRow, "PercentageCoverage"]
  }
  
  return(tipInfo)
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
