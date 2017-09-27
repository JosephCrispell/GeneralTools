#################
# Load Packages #
#################

library(gplots) # Heatmap

################
# Read in data #
################

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/InvestigatingCattleMislabelling/Spoligotyping/"

# Read in the Spoligotype conversion table
file <- paste(path, "SpoligotypeConversionTable_17-01-17.txt", sep="")
conversionTable <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                              colClasses="character")

# Read in the FASTA file containing the spacer primers
file <- paste(path, "Xia2016-43SpacerSequences-25bp.fasta", sep="")
sequences <- readFastaFile(file)

# Note the reference spacer binary code
referenceBinaryCode <- getSpoligotypeBinaryCode(
  "1100101000001000011100110111101111111100000")

# Note the spacer binary codes of spoligotypes 9 and 17
spoligotypes <- list(
  "9" = getSpoligotypeBinaryCode(
    "1101101000001110111111111111111111111100000"),
  "17" = getSpoligotypeBinaryCode(
    "1101101000001110111111111000000111111100000")
)

# Note the spacer counts for the badger and cattle isolates
isolateSpacerCounts <- c(315, 318, 0, 317, 315, 0, 317, 0, 0, 0, 0, 0, 315, 316, 313, 0, 316, 316, 318, 319, 320, 315, 317, 312, 318, 63, 62, 63, 63, 61, 63, 302, 314, 314, 315, 311, 315, 314, 0, 0, 0, 0, 0)

# Open a PDf
file <- paste(path, "SpoligotypeSpacers_30-05-17.pdf", sep="")
pdf(file)

########################################
# Examine the frequency of each spacer #
########################################

# Count the occurrence of each spacer in the conversion table
sum <- countSpacerOccurrenceInConversionTable(conversionTable)

# Plot spacer info
plotSpacers(sum, "Spacers present in database")
plotSpacers(referenceBinaryCode, "Spacers present in reference")
plotSpacers(isolateSpacerCounts, "Spacers present in isolates")

# Compare spacers present in reference to those in isolates
spacerCounts <- list(
  "Sum" = sum,
  "Reference" = referenceBinaryCode,
  "Isolates" = isolateSpacerCounts
)
plotSpacersComparison(spacerCounts, "Compare spacer counts in database, reference, and isolates")

# Compare spoligotypes 9 and 17
plotSpacersComparison(spoligotypes, "Compare spoligotypes 9 and 17")

#####################################################
# Examine the genetic distances between the spacers #
#####################################################

# Calculate the genetic distances
geneticDistances <- buildGeneticDistanceMatrix(sequences)

# Plot a heatmap showing the genetic distances between the spacer's primers
plotHeatmap(geneticDistances, "Genetic distances between spacers")

dev.off()

#############
# FUNCTIONS #
#############

plotSpacersComparison <- function(listOfArrays, title){
  
  par(mar=c(1,3,4,0)) # Defaults: 5.1,4.1,4.1,2.1 bottom, left, top, right
  
  # Divide all values by the max
  names <- names(listOfArrays)
  for(key in names){
    listOfArrays[[key]] <- normalise(listOfArrays[[key]])
  }

  # Create the plot
  plot(x=NULL, y=NULL, xlim=c(0.5,43.5), ylim=c(0,1), yaxt="n", xaxt="n", ylab="",
       xlab="", bty="n", main=title)
  
  # Add lower axis
  axis(side=1, at=1:43, cex.axis=0.5, las=2, tick=FALSE, line=-1.5)
  
  # Add coloured bars to represent each spacer
  for(spacerIndex in 1:43){
    
    for(i in 1:length(names)){
      
      array <- listOfArrays[[names[i]]]
      
      polygon(x=c(spacerIndex - 0.5, spacerIndex - 0.5,
                  spacerIndex + 0.5, spacerIndex + 0.5), 
              y=c((1 / length(names)) * (i-1), (1 / length(names)) * i,
                  (1 / length(names)) * i, (1 / length(names)) * (i-1)),
              col=rgb(0,0,0, array[spacerIndex]), border=NA)
    }
  }
  
  # Add lines to separate each spacer
  for(i in 1:43){
    lines(x=c(i - 0.5, i - 0.5), y=c(0,1), col="red")
    lines(x=c(i + 0.5, i + 0.5), y=c(0,1), col="red")
  }
  
  # Add lines to separate
  for(i in 1:(length(names)-1)){
    lines(x=c(0,44), y=c((1 / length(names)) * i, (1 / length(names)) * i),
          lwd=2, col="white")
  }
  
  # Label the spacer counts
  for(i in 1:length(names)){
    
    text(x=-6, y=((1 / length(names)) * i) - ((1 / length(names))/2) ,
         labels=names[i], xpd=TRUE, pos=4)
  }
  
  par(mar=c(5.1,4.1,4.1,2.1))
}

plotHeatmap <- function(matrix, title){
  colBreaks <- seq(0, 1, by=0.1)
  
  # Plot the heatmap
  heatmap.2(matrix, # matrix is the input data
            
            # Show values
            cellnote=matrix,
            notecex=0.35,
            notecol="black",
            
            # Create the colour scale 
            col=colorpanel(n=length(colBreaks)-1, low="white", high="red"),
            breaks=colBreaks,
            
            # Turn off a density plot
            density.info="none", 
            
            # Turn off the trace
            trace="none",
            
            # Column Labels
            labCol=colnames(matrix),
            cexCol=1, # Change the size of the column labels
            srtCol=90, # Set the angle of the column labels (degrees from horizontal)
            offsetCol=-0.85, # Set size of space between column labels and heatmap
            
            # Change the size of the margins around the plot: c(column space, row space)
            margins = c(5, 5), 
            
            # Row labels
            labRow=rownames(matrix),
            cexRow=1, # Change the size of the Row labels
            offsetRow=0,
            
            # Make sure the order of the rows and columns is changed
            Rowv=TRUE, Colv=TRUE,
            
            # Don't plot any dendogram
            dendrogram="none",
            
            # Set up the Key
            key=FALSE, # Turn the key OFF
            
            # Say where to plot each part of the heatmap
            #     4     3
            #     2     1
            # 1. Heatmap
            # 2. Row Dendrogram
            # 3. Column Dendrogram
            # 4. Key
            lmat=rbind(4:3, 2:1),
            
            # Set the size of the spaces for output plots:
            #     Colour Key      |   Column Dendrogram
            #     -------------------------------------
            #     Row Dendrogram  |   Heatmap
            #
            # Note that these will be affected by the width and height you set for the PDF
            lhei=c(1, 10), # c(row1Width, row2Width)
            lwid=c(0.1, 10), # c(column1Width, column2Width)
            
            # Note that the input matrix is not symmetric
            symm = FALSE
  )
  
  text(x=18, y=1.14, labels=title, xpd=TRUE,
       cex=1)
}

getSpacerNumbers <- function(array){
  
  numbers <- c()
  for(i in 1:length(array)){
    numbers[i] <- substring(array[i], 7)
  }
  
  return(numbers)
}

buildGeneticDistanceMatrix <- function(sequences){
  
  # Get the sequence ids
  ids <- getSpacerNumbers(names(sequences))
  fullIDs <- names(sequences)
  
  # Initialise a matrix to store the genetic distances
  distances <- matrix(nrow=length(ids), ncol=length(ids))
  colnames(distances) <- ids
  rownames(distances) <- ids
  
  # Calculate the genetic distances
  for(i in 1:length(ids)){
    for(j in 1:length(ids)){
      
      # Skip self-comparisons and only use the upper triangle
      if(i >= j){
        next
      }
      
      # Calculate the genetic distance
      distance <- geneticDistance(sequences[[fullIDs[i]]],
                                  sequences[[fullIDs[j]]], proportion=TRUE)
      distances[i, j] <- distance
      distances[j, i] <- distance
    }
  }
  
  return(distances)
}

geneticDistance <- function(a, b, proportion){
  
  dist <- 0
  
  for(i in 1:length(a)){
    
    if(a[i] != "N" && b[i] != "N" && a[i] != b[i]){
      dist <- dist + 1
    }
  }
  
  if(proportion == TRUE){
    dist <- dist / length(a)
  }
  
  return(dist)
}

readFastaFile <- function(fileName){
  
  # Read in the lines in the file
  connection <- file(fileName, open="r")
  fileLines <- readLines(connection)
  close(connection)
  
  # Store each sequence
  sequences <- list()
  for(i in 2:length(fileLines)){
    
    # Have we reached another sequence?
    if(grepl(x=fileLines[i], pattern=">") == TRUE){
      
      # Store the previous sequence
      if(i > 2){
        sequences[[name]] <- strsplit(sequence, split="")[[1]]
      }
      
      # Start the new sequence
      name <- substring(fileLines[i], 2)
      sequence <- ""
    }else{
      sequence <- paste(sequence, fileLines[i], sep="")
    }
  }
  
  # Store the last sequence
  sequences[[name]] <- strsplit(sequence, split="")[[1]]
  
  return(sequences)
}

getSpoligotypeBinaryCode <- function(string){
  
  return(as.numeric(strsplit(string, split="")[[1]]))
}

countSpacerOccurrenceInConversionTable <- function(conversionTable){
  sum <- rep(0, 43)
  for(row in 1:nrow(conversionTable)){
    
    binaryCode <- getSpoligotypeBinaryCode(conversionTable[row, "Binary"])
    
    sum <- sum + binaryCode
  }
  
  return(sum)
}

plotSpacers <- function(array, title){
  
  par(mar=c(1,0,4,0)) # Defaults: 5.1,4.1,4.1,2.1 bottom, left, top, right
  
  # Divide all values by the max
  array <- normalise(array)
  
  # Create the plot
  plot(x=NULL, y=NULL, xlim=c(0.5,43.5), ylim=c(0,1), yaxt="n", xaxt="n", ylab="",
       xlab="", bty="n", main=title)
  
  # Add lower axis
  axis(side=1, at=1:43, cex.axis=0.5, las=2, tick=FALSE, line=-1.5)
  
  # Add coloured bars to represent each spacer
  for(i in 1:length(array)){
    polygon(x=c(i - 0.5, i - 0.5, i + 0.5, i + 0.5), 
            y=c(0, 1, 1, 0), col=rgb(0,0,0, array[i]), border=NA)
  }
  
  # Add lines to separate each spacer
  for(i in 1:43){
    lines(x=c(i - 0.5, i - 0.5), y=c(0,1), col="red")
    lines(x=c(i + 0.5, i + 0.5), y=c(0,1), col="red")
  }
  
  par(mar=c(5.1,4.1,4.1,2.1))
}

normalise <- function(array){
  
  return(array / max(array))
}
