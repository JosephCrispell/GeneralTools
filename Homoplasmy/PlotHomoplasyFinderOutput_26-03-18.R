#~~~~~~~~~~~~~~~~~~~#
#### Preparation ####
#~~~~~~~~~~~~~~~~~~~#

# Packages
library(ape)

# Set the path
#path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/NewZealand/NewAnalyses_12-05-16/MLTree/"
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/vcfFiles/"

#~~~~~~~~~~~~~~~~~~~~#
#### Read in data ####
#~~~~~~~~~~~~~~~~~~~~#

# Read in the tree file
treeFile <- paste(path, "mlTree_27-03-18.tree", sep="")
tree <- read.tree(file=treeFile)

# Read in the HomoplasyFinder output
homoplasyFinderOutputFile <- paste(path, "homoplasyReport_27-03-18.txt", sep="")
results <- read.table(homoplasyFinderOutputFile, header=TRUE, sep="\t", comment.char="@")

# Read in the FASTA file - only keep homoplasy site information
fastaFile <- paste(path, "sequences_Prox-10_24-03-2018.fasta", sep="")
sequences <- readFASTA(fastaFile)
sequences <- keepPositions(sort(results$Position), sequences)
sequences <- parseNewZealandSequenceNames(sequences)

#~~~~~~~~~~~~~~~~~~~~~~~~#
#### Plot the results ####
#~~~~~~~~~~~~~~~~~~~~~~~~#

file <- paste(strsplit(homoplasyFinderOutputFile, split="[.]")[[1]][1], ".pdf", sep="")
pdf(file)

plotTreeAndHomoplasySites(tree, sequences, cex=1, sort(results$Position))

dev.off()


#############
# FUNCTIONS #
#############

plotTreeAndHomoplasySites <- function(tree, sequences, cex, positions){
  
  # Set the plotting margins
  par(mar=c(0,0,1,0.5))
  
  # Plot the phylogeny
  plot.phylo(tree, show.tip.label=TRUE, type="phylogram", align.tip.label=TRUE, tip.color=rgb(0,0,0,0),
             main="Homoplasies Identified")
  
  # Note the locations of the tips
  tipCoordinates <- getTipCoordinates(tree$tip.label)
  maxCoords <- maxCoordinates(tipCoordinates)
  
  # Calculate width of space for nucleotide
  axisLimits <- par("usr")
  charWidth <- (axisLimits[2] - maxCoords[1]) / length(sequences[[1]])
  
  # Set nucleotide colours
  nucleotideColours <- list("A"="red", "C"="blue", "G"="cyan", "T"="orange")
  
  # Plot FASTA alignment beside tree
  ids <- names(sequences)
  for(i in 1:length(ids)){
    
    # Get sequence
    sequence <- sequences[[ids[i]]]
    
    # Get XY coordinates fo tip
    xy <- tipCoordinates[[ids[i]]]
    
    # Plot alignment
    for(i in 1:length(sequence)){
      if(sequence[i] != "N"){
        polygon(x=c(maxCoords[1] + ((i-1) * charWidth),
                    maxCoords[1] + ((i-1) * charWidth),
                    maxCoords[1] + (i * charWidth),
                    maxCoords[1] + (i * charWidth)),
                y=c(xy[2], xy[2] + 1, xy[2] + 1, xy[2]),
                col=nucleotideColours[[sequence[i]]],
                border=rgb(0,0,0,0), xpd=TRUE)
      }
    }
  }
  
  # Add separator lines
  for(i in 1:length(sequences[[1]])){
    points(x=c(maxCoords[1] + ((i-1) * charWidth),
               maxCoords[1] + ((i-1) * charWidth)),
           y=c(axisLimits[3], axisLimits[4]),
           type="l", col="white")
  }
  
  # Note the positions of each homoplasy
  charHeight <- strheight("A")
  for(i in 1:length(positions)){
    text(x=maxCoords[1] + ((i-1) * charWidth) + (0.5 * charWidth),
         y=maxCoords[2] + charHeight,
         labels=positions[i], cex=0.5, srt=90, xpd=TRUE)
  }
  
  # Reset plotting margins
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

maxCoordinates <- function(tipCoordinates){
  
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

parseNewZealandSequenceNames <- function(sequences){
  
  output <- list()
  for(id in names(sequences)){
    
    parts <- strsplit(id, "_")[[1]]
    newId <- parts[1]
    if(grepl(id, pattern="#") == TRUE){
      newId <- paste(parts[1], parts[2], sep="_")
    }
    
    output[[newId]] <- sequences[[id]]
    
  }
  
  return(output)
}

keepPositions <- function(positions, sequences){
  
  for(id in names(sequences)){
    sequences[[id]] <- sequences[[id]][positions]
  }
  
  return(sequences)
}

getTipCoordinates <- function(tipLabels){
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  tips <- list()
  for(i in 1:length(tipLabels)){
    tips[[as.character(tipLabels[i])]] <- c(lastPP$xx[i], lastPP$yy[i])
  }
  
  return(tips)
}

readFASTA <- function(fastaFile){
  
  # Open the file and store lines
  connection <- file(fastaFile, "r")
  fileLines <- readLines(connection)
  close(connection)
  
  # Initialise a list to store the sequences
  sequences <- list()
  
  # Examine each of the file lines one by one
  skip <- 0
  for(i in 1:length(fileLines)){
    
    if(i == 1 && grepl(fileLines[i], pattern="^>") == FALSE){
      skip <- 1
      next
      
    # Check if line starts with ">"
    }else if(grepl(fileLines[i], pattern="^>") == TRUE){
      
      # Store the previous sequence
      if(i != 1 + skip){
        sequences[[name]] <- strsplit(sequence, split="")[[1]]
      }
      
      # Initialise variables to store sequence and id
      name <- substr(fileLines[i], 2, stop=nchar(fileLines[i]))
      sequence <- ""
    }else{
      
      sequence <- paste(sequence, fileLines[i], sep="")
    }
  }
  
  # Store the last sequence
  sequences[[name]] <- strsplit(sequence, split="")[[1]]
  
  return(sequences)
}
