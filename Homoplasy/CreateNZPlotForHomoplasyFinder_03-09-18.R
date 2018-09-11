#### Load libraries ####
library("devtools")
install_github("JosephCrispell/homoplasyFinder")
library(homoplasyFinder)
library(addTextLabels)

#### Run HomoplasyFinder ####

# Note the input files
path <- "/home/josephcrispell/Desktop/Research/Homoplasy/DataForTesting/"
fastaFile <- paste0(path, "sequences_withRef_Prox-10_14-06-16_NZ.fasta")
treeFile <- paste0(path, "mlTree_withRef_14-06-16_NZ.tree")

# Run HomoplasyFinder
inconsistentPositions <- runHomoplasyFinderInJava(treeFile, fastaFile, path)

# Get the annotated tree
annotatedTree <- readAnnotatedTree(path)

#### Change the fasta position annotations on tree ####

# Get the genome positions of each position in FASTA
fastaPositions <- read.table(paste0(path, "fastaPositions_Prox-10_14-06-16_NZ.txt"), header=TRUE, sep="\t")

# Change labels on annotated tree
annotatedTree <- changeFastaPositionsToGenomePositions(annotatedTree, inconsistentPositions, fastaPositions)

#### Plot the annotated tree ####

# Get genome positions of inconsistent positions
truePositions <- fastaPositions$Position[inconsistentPositions]

# Plot the annotated tree
pdf(paste0(path, "HomoplasyFinder_AnnotatedTree_11-09-18.pdf"))
plotAnnotatedTree(annotatedTree, inconsistentPositions, fastaFile, addScale=FALSE, truePositions, nodeLabelCex=0.7)
dev.off()

#### FUNCTIONS ####

changeFastaPositionsToGenomePositions <- function(annotatedTree, inconsistentPositions, fastaPositions){
  
  # Get genome positions of inconsistent positions
  genomePositions <- fastaPositions$Position[inconsistentPositions]
  
  # Create list: key=fastaPosition, value:genomePosition
  positions <- list()
  for(i in seq_along(inconsistentPositions)){
    positions[[as.character(inconsistentPositions[i])]] <- genomePositions[i]
  }
  
  # Replace node labels with genome positions
  for(i in seq_along(annotatedTree$node.label)){
    
    # Skip blank labels
    if(annotatedTree$node.label[i] == ""){
      next
    }
    
    # Replace fasta position with genome position
    annotatedTree$node.label[i] <- positions[[annotatedTree$node.label[i]]]
  }
  
  return(annotatedTree)
}

plotAnnotatedTree <- function(tree, inconsistentPositions, fastaFile, addScale=TRUE, truePositions, nodeLabelCex){
  
  # Get the current plotting margins - so that we can revert to these once finished
  marginSettings <- par("mar")
  
  # Set the plotting margins
  par(mar=c(0,0,1,0.5))
  
  # Plot the phylogeny
  # Note plotting invisible tip labels - these provide space for plotting alignment
  ape::plot.phylo(tree, show.tip.label=TRUE, type="phylogram", align.tip.label=TRUE, 
                  tip.color=rgb(0,0,0,0), main="Homoplasies Identified", cex=1)
  
  # Add a scale bar if requested
  addScaleBar(addScale)
  
  # Add internal node labels indicating the inconsistent positions associated with each internal node
  addInternalNodeLabels(tree, cex=nodeLabelCex)
  
  # Add an alignment detailing the nucleotides for each inconsistent position
  addInconsistentPositionAlignment(tree, fastaFile, inconsistentPositions, truePositions)
  
  # Revert the previous margin settings
  par("mar"= marginSettings)
}

getTipCoordinates <- function(tipLabels){
  lastPP <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  tips <- list()
  for(i in seq_along(tipLabels)){
    tips[[as.character(tipLabels[i])]] <- c(lastPP$xx[i], lastPP$yy[i])
  }
  
  return(tips)
}

getInternalNodeCoordinates <- function(nTips){
  
  # Get all the information from the most recent plot.phylo() call
  lastPP <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  
  # Initialise a matrix to store X and Y coordinates of the internal nodes
  coordinates <- matrix(NA, nrow=length(lastPP$xx) - nTips, ncol=2)
  for(i in (nTips + 1):length(lastPP$xx)){
    coordinates[(i - nTips), 1] <- lastPP$xx[i]
    coordinates[(i - nTips), 2] <- lastPP$yy[i]
  }
  
  return(coordinates)
}

addInternalNodeLabels <- function(tree, cex){
  
  # Get the coordinates of the internal nodes
  internalNodeCoords <- getInternalNodeCoordinates(length(tree$tip.label))
  
  # Note the indices of labels to plot
  indices <- c()
  for(i in 1:length(tree$node.label)){
    
    # Skip internal nodes with no label
    if(tree$node.label[i] == ""){
      next
    }
    
    # Add the current index
    indices[length(indices) + 1] <- i
  }
  
  # Add labels for each internal node
  addTextLabels(xCoords=internalNodeCoords[indices, 1], yCoords=internalNodeCoords[indices, 2],
                labels=tree$node.label[indices], cex=cex, col.label="white", col.background=rgb(0,0,0, 0.75), col.line="red")
}

getTipSequences <- function(fastaFile){
  
  # Read in the fasta file
  sequences <- ape::read.dna(fastaFile, format="fasta", as.character=TRUE)
  
  # Create a list recording each sequence
  output <- list()
  for(row in rownames(sequences)){
    output[[row]] <- sequences[row, ]
  }
  
  return(output)
}

getMaxCoordinates <- function(tipCoordinates){
  
  max <- c(NA, NA)
  
  for(key in names(tipCoordinates)){
    
    if(is.na(max[1]) == TRUE || tipCoordinates[[key]][1] > max[1]){
      max[1] <- tipCoordinates[[key]][1]
    }
    
    if(is.na(max[2]) == TRUE || tipCoordinates[[key]][2] > max[2]){
      max[2] <- tipCoordinates[[key]][2]
    }
  }
  
  return(max)
}

addScaleBar <- function(addScale){
  
  # Get the axis Limits
  axisLimits <- par("usr")
  
  # Add Scale bar
  if(addScale == TRUE){
    xLength <- axisLimits[2] - axisLimits[1]
    yLength <- axisLimits[4] - axisLimits[3]
    points(x=c(axisLimits[1] + 0.05*xLength, axisLimits[1] + 0.15*xLength),
           y=c(axisLimits[3]+0.05*yLength, axisLimits[3]+0.05*yLength), type="l", lwd=3)
    text(x=axisLimits[1] + 0.1*xLength, y=axisLimits[3] +0.07*yLength,
         labels=round(0.1*xLength, digits=1), cex=1, xpd=TRUE)
  }
}

addInconsistentPositionAlignment <- function(tree, fastaFile, inconsistentPositions, truePositions){
  
  # Get the axis Limits
  axisLimits <- par("usr")
  
  # Read in the nucleotide sequences
  tipSequences <- getTipSequences(fastaFile)
  
  # Note the locations of the tips
  tipCoordinates <- getTipCoordinates(tree$tip.label)
  maxCoords <- getMaxCoordinates(tipCoordinates)
  
  # Calculate width of space for nucleotide
  charWidth <- (axisLimits[2] - maxCoords[1]) / length(inconsistentPositions)
  
  # Set nucleotide colours
  nucleotideColours <- list("a"="red", "c"="blue", "g"="cyan", "t"="orange")
  
  # Plot FASTA alignment beside tree
  for(positionIndex in seq_along(inconsistentPositions)){
    
    # Examine each of the tips
    for(tipIndex in seq_along(tree$tip.label)){
      
      # Get XY coordinates for the current tip
      xy <- tipCoordinates[[tree$tip.label[tipIndex]]]
      
      # Get the sequence for the current tip
      sequence <- tipSequences[[tree$tip.label[tipIndex]]]
      
      # Plot a polygon for the current tip's nucleotide at the current homoplasy's position
      polygon(x=c(maxCoords[1] + ((positionIndex-1) * charWidth),
                  maxCoords[1] + ((positionIndex-1) * charWidth),
                  maxCoords[1] + (positionIndex * charWidth),
                  maxCoords[1] + (positionIndex * charWidth)),
              y=c(xy[2]-0.5, xy[2] + 0.5, xy[2] + 0.5, xy[2]-0.5),
              col=nucleotideColours[[sequence[inconsistentPositions[positionIndex]]]],
              border=rgb(0,0,0,0), xpd=TRUE)
      
    }
  }
  
  # Add separator lines
  for(positionIndex in seq_along(inconsistentPositions)){
    points(x=c(maxCoords[1] + ((positionIndex-1) * charWidth),
               maxCoords[1] + ((positionIndex-1) * charWidth)),
           y=c(axisLimits[3], axisLimits[4]),
           type="l", col="white", xpd=TRUE)
  }
  
  # Calculate the length of homoplasy position label
  charHeight <- strheight("o") * (max(nchar(truePositions)) / 4)
  maxWidth <- max(strwidth(truePositions))
  
  # Note the positions of each homoplasy
  for(positionIndex in seq_along(truePositions)){
    text(x=maxCoords[1] + ((positionIndex-1) * charWidth) + (0.5 * charWidth),
         y=maxCoords[2] + charHeight,
         labels=truePositions[positionIndex],
         cex=0.5, srt=90, xpd=TRUE)
  }
}
