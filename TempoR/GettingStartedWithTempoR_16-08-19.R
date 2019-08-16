#### libraries ####

library(ape) # Reading in phylogeny
library(adePhylo)

#### Read in tree and get dates ####

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

# Create a path variable
path <- "/home/josephcrispell/Desktop/Research/tempoR/"

# Read in the tree file
tree <- read.tree(paste0(path, "example_dated.tree"))

# Get dates from tip labels
tipDates <- getTipDatesFromTipLabels(tree$tip.label, dateFormat="%Y", dateColumn=3)

#### Find best root for phylogeny ####

# Get the best rooted tree
bestRootedTreeInfo <- identifyBestRootForSamplingDates(tree, tipDates, plot=TRUE,
                                                       las=1, bty="n", pch=19, col=rgb(0,0,0, 0.5))

# Create a root-to-tip plot
plotTipToRootDistances(bestRootedTreeInfo$tipDates, bestRootedTreeInfo$tipToRootDistances, 
                       lineOfBestFirstColour=rgb(1,0,0, 0.75), lineOfBestFitType=2,
                       las=1, bty="n", pch=19, col=rgb(0,0,0, 0.5))

#### Plot root-to-tip distances versus sampling date ####

#### FUNCTIONS ####

identifyBestRootForSamplingDates <- function(tree, tipDates, plot=FALSE, ...){
  
  # Add the tip dates to the tip labels
  tree$tip.label <- addDatesToTipLabels(tree$tip.label, tipDates)
  
  # Note each internal numbers internal node
  nTips <- length(tree$tip.label)
  internalNodeIndices <- seq(from=nTips+1, to=nTips+tree$Nnode)
  
  # Root the tree on each of the internal nodes
  rSquaredValues <- sapply(internalNodeIndices, FUN=rootTreeAndFitLinearModel, tree)

  # Note the index of the best internal node to root by
  indexOfBest <- which.max(rSquaredValues)

  # Create the best rooted tree
  reRootedTree <- root(tree, node=internalNodeIndices[indexOfBest])
  tipDates <- as.Date(getLastPartOfStrings(reRootedTree$tip.label, sep="_"), format="%Y-%m-%d")
  reRootedTree$tip.label <- removeDatesFromTipLabels(reRootedTree$tip.label)
  
  # Calculate the distances from the tips to the root
  tipToRootDistances <- sapply(1:nTips, FUN=calculateTipToRootDistance, reRootedTree)
  
  # Create summary plots if requested
  if(plot){
    
    # Plot the rSquared values values calculated for each internal node
    plot(x=internalNodeIndices, y=rSquaredValues, ylab="rSq from tip-to-root vs. sampling dates", xlab="Internal node index", 
         main="Change in rSq with different rooted trees",...)
    
    # Plot of tip to root distances versus sampling time
    plotTipToRootDistances(tipDates, tipToRootDistances, ...)
  }
  
  # Create a list to store the rooted tree, rSq value, tip dates
  bestTreeInfo <- list("tree"=reRootedTree, "rSquared"=rSquaredValues[indexOfBest], "tipDates"=tipDates,
                       "tipToRootDistances"=tipToRootDistances)
  
  # Return the rooted tree with the best Rsq value
  return(bestTreeInfo)
}

plotTipToRootDistances <- function(tipDates, tipToRootDistances, main="Examining temporal signature", 
                                   lineOfBestFirstColour=rgb(0,0,0,0.75), lineOfBestFitWidth=2, lineOfBestFitType=1, ...){
  
  # Fit a linear model
  linearModel <- lm(tipToRootDistances ~ tipDates)
  summary <- summary(linearModel)
  
  # Plot the tip-to-root distances versus sampling date
  plot(x=tipDates, y=tipToRootDistances, ylab="Tip-to-root distance", xlab="Sampling date",
       main=main, ...)
  
  # Add a line of best fit
  abline(linearModel, col=lineOfBestFirstColour, lwd=lineOfBestFitWidth, lty=lineOfBestFitType)
  
  # Add a legend to report rSq value
  legend("topleft", 
         legend=c(paste("R^2 = ", round(summary$adj.r.squared, 2)),
                  paste("p-value = ", round(anova(linearModel)$Pr[[1]], 2))),
         bty='n')
}

removeDatesFromTipLabels <- function(tipLabels){
  
  # Intialise a vector to store the tipLabels without dates
  labelsWithoutDates <- c()
  
  # Examine each tip label
  for(index in seq_along(tipLabels)){
    
    # Split the current label into its parts
    parts <- strsplit(tipLabels[index], split="_")
    
    # Create a label without date
    labelsWithoutDates[index] <- paste(parts[1:(length(parts)-1)], collapse="_")
  }
  
  return(labelsWithoutDates)
}

addDatesToTipLabels <- function(tipLabels, dates){
  
  # Intialise a vector to store the tipLabels with dates
  labelsWithDates <- c()
  
  # Examine each tip label
  for(index in seq_along(tipLabels)){
    labelsWithDates[index] <- paste0(tipLabels[index], "_", dates[index])
  }
  
  return(labelsWithDates)
}

rootTreeAndFitLinearModel <- function(nodeIndex, tree){

  # Re-root the tree
  reRootedTree <- root(tree, node=nodeIndex)
  
  # Calculate the distance of each tip to the root
  tipToRootDistances <- sapply(1:length(tree$tip.label), FUN=calculateTipToRootDistance, reRootedTree)
  
  # Get the dates from the tip labels
  # NOTE: dates will be in default format
  tipDates <- as.Date(getLastPartOfStrings(reRootedTree$tip.label, sep="_"))
  
  # Calculate the corrlation between the sampling dates
  linearModel <- lm(tipToRootDistances ~ tipDates)
  summaryOfLinearModel <- summary(linearModel)
  
  # Note the direction of the relationship (negative or positive)
  direction <- sign(summaryOfLinearModel$coefficients[2, 1])
  
  return(direction*summaryOfLinearModel$adj.r.squared)
}

getLastPartOfStrings <- function(strings, sep){
  
  # Intialise avector to store the last part of each string
  ends <- c()
  
  # Examine each string
  for(index in seq_along(strings)){
    
    # Split the current string into its parts
    parts <- strsplit(strings[index], split=sep)[[1]]
    
    # Store the end of the current string
    ends[index] <- parts[length(parts)]
  }
  
  return(ends)
}

calculateTipToRootDistance <- function(tipIndex, tree){
  
  # Identify the index of the root
  rootIndex <- length(tree$tip.label) + 1
  
  # Initialise a variable to store the sum of the branch lengths back through the edges
  sum <- 0
  
  # Trace the path back to the root
  currentNodeIndex <- tipIndex
  while(currentNodeIndex != rootIndex){
    
    # Get the index of the edge to the current node
    edgeIndex <- which(tree$edge[, 2] == currentNodeIndex)
    
    # Add to the growing sum of branch lengths
    sum <- sum + tree$edge.length[edgeIndex]
    
    # Update the index of the current node
    currentNodeIndex <- tree$edge[edgeIndex, 1]
  }
  
  return(sum)
}

getTipDatesFromTipLabels <- function(tipLabels, dateFormat, dateColumn){
  
  # Intialise a vector to store dates
  tipDates <- c()
  
  # Get the date strings
  for(tipIndex in seq_along(tipLabels)){
    
    tipDates[tipIndex] <- strsplit(tipLabels[tipIndex], split="_")[[1]][dateColumn]
  }
  
  # Convert date strings to dates
  tipDates <- as.Date(tipDates, format=dateFormat)
  
  return(tipDates)
}
