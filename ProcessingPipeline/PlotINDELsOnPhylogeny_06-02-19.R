#### Load libraries ####
library(ape)

#### Set the path ####

path <- "/home/josephcrispell/Desktop/Research/RepublicOfIreland/Mbovis/WorkingOnIndels/"

#### Read in the INDEL presence/absence matrix ####

# Read in the presence/absence matrix
file <- paste0(path, "indelInfo_presenceAbsence_7-2-2019.txt")
presenceAbsence <- read.table(file, header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE, na.strings="F")

# Remove INDELs are found in all or none
presenceAbsence <- removeINDELSThatShowBelowMinDifference(presenceAbsence, minDiff=2)

# Remove INDEL IDs and make them the row names
rownames(presenceAbsence) <- presenceAbsence[, 1]
presenceAbsence <- presenceAbsence[, -1]

# Parse the sequence names in the columns
colnames(presenceAbsence) <- parseSequenceLabels(colnames(presenceAbsence))

# Get the transpose of the presence/absence matrix
presenceAbsence <- t(presenceAbsence)

#### Read in the phylogeny ####

# Read in the RAxML tree file
treeFile <- paste0(path, "RAxML_bipartitions.RaxML-R_07-02-19")
tree <- read.tree(treeFile)

# Note the order of the tips in the tree and remove reference
tips <- parseSequenceLabels(tree$tip.label)
tips <- tips[which(tips != "Ref-1997")]

#### Plot heatmap ####

plotHeatmap(presenceAbsence)
plotHeatmap(presenceAbsence, order=tips)

library(gplots)
heatmap.2(presenceAbsence,
          Colv=FALSE,
          dendrogram="row",
          trace="none",
          labCol="",
          na.color="black")

#### FUNCTIONS ####

parseSequenceLabels <- function(labels){
  
  # Initialise a vector to store the parsed labels
  parsedLabels <- c()
  
  # Examine each of the labels
  for(i in seq_along(labels)){
    
    # Remove the ">" if present at start of current label
    parsed <- labels[i]
    if(startsWith(parsed, ">")){
      parsed <- substr(parsed, 2, nchar(parsed))
    }
    
    # Extract sequence ID
    parsedLabels[i] <- strsplit(parsed, split="_")[[1]][1]
  }
  
  return(parsedLabels)
}

sortSequencesBySimilarity <- function(presenceAbsence){
  
  # Calculate the distance between the sequences based on the presence/absence of INDELS
  distances <- dist(presenceAbsence, method = "euclidean") # Calculates distances between rows
  
  # Conduct hierarchical clustering
  clusters <- hclust(distances)
  
  return(clusters$order)
}

plotHeatmap <- function(presenceAbsence, order=NULL){
  
  # Order the sequences in presence/absence table by similarity if no order specified
  if(is.null(order)){
    order <- sortSequencesBySimilarity(presenceAbsence)
    presenceAbsence <- presenceAbsence[order, ]
  }else{
    presenceAbsence <- presenceAbsence[order, ]
  }
  
  # Set the margins
  par(mar=c(0,5,0,0)) # par(mar=c(5.1, 4.1, 4.1, 2.1))
  
  # Create an empty plot
  plot(x=NULL, y=NULL, xlim=c(0, ncol(presenceAbsence)), ylim=c(0, nrow(presenceAbsence)),
       yaxt="n", xaxt="n", bty="n", ylab="", xlab="")
  
  # Plot each of the cells in the presence/absence matrix
  for(row in seq_len(nrow(presenceAbsence))){
    
    for(col in 1:ncol(presenceAbsence)){
      
      # Define the colour - red=present, white=absent
      colour <- "white"
      if(is.na(presenceAbsence[row, col]) == FALSE && presenceAbsence[row, col] == 1){
        colour <- "red"
      }else if(is.na(presenceAbsence[row, col]) == FALSE && presenceAbsence[row, col] == 0){
        colour <- "blue"
      }
      
      # Plot the current cell
      rect(xleft=col-1, ybottom=row-1, xright=col, ytop=row, col=colour, border=rgb(0,0,0, 0))
    }
  }
  
  # Add the sequence labels
  labels <- rownames(presenceAbsence)
  mtext(text=labels, side=2, at=seq(from=0.5, to=length(labels)), las=1, line=-1)
  
  # Add grid
  for(y in 0:nrow(presenceAbsence)){
    lines(x=c(0, ncol(presenceAbsence)), y=c(y,y))
  }
  for(x in 0:ncol(presenceAbsence)){
    lines(x=c(x, x), y=c(0,nrow(presenceAbsence)))
  }
  
  # Reset the margins
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

removeINDELSThatShowBelowMinDifference <- function(presenceAbsence, minDiff){
  
  # Initialise a vector to store the indices of the rows (INDELS) that aren't different across the isolates
  rowsToIgnore <- c()
  
  # Examine each row of the presence/absence table
  for(row in seq_len(nrow(presenceAbsence))){
    
    # Create a variable to record if difference observed
    nDifferences <- 0
    
    # Create a variable to record the first non-NA value
    firstValue <- NA
    
    # Examine whether current INDEL is present/absent in each isolate - compare to first non-NA isolate
    for(value in presenceAbsence[row, 2:ncol(presenceAbsence)]){
      
      # Check if current value is non-NA and we haven't encountered a non-NA value
      if(is.na(firstValue) && is.na(value) == FALSE){
        firstValue <- value
        next
      }
      
      # Compare current value to the first - if it isn't an NA. Record if different and finish
      if(is.na(value) == FALSE && value != firstValue){
        nDifferences <- nDifferences + 1
      }
    }
    
    # Check if no difference was found
    if(nDifferences < minDiff){
      rowsToIgnore[length(rowsToIgnore) + 1] <- row
    }
  }
  
  # Remove the INDELs without a difference
  if(length(rowsToIgnore) > 0){
    
    # Report how many INDELs will be removed
    cat(paste0("Removed ", length(rowsToIgnore), " (of ", nrow(presenceAbsence),
               ") INDELs that showed less than ", minDiff, " difference(s) across the sequences"))
    
    # Remove the INDELS without differences
    presenceAbsence <- presenceAbsence[-rowsToIgnore, ]
  }else{
    cat("No INDELs without differences across the sequeces were found!")
  }
  
  return(presenceAbsence)
}
