# Load the necessary libraries
library(ape) # For plotting phylogenetic trees
library(grid) # Used to plot lines between plot panels

# Set the seed
set.seed(244843)

# Generate a random tree
randomTree <- rtree(n=50)

# Copy the tree, but rotate a node
treeWithRotatedNode <- rotate(randomTree, 75)

# Remove the plot margins
par(mar=c(0,0,0,0))

# Allow two columns of plots
par(mfrow=c(1,2))

# Plot the random tree and get its coordinates
plot.phylo(randomTree, show.tip.label=FALSE, edge.width=2, direction="rightwards")
tipCoordinates <- getTipCoordinates(randomTree$tip.label)

# Plot the tree with the rotated node and get its tip coordinates
plot.phylo(treeWithRotatedNode, show.tip.label=FALSE, edge.width=2, direction="leftwards")
tipCoordinatesRotated <- getTipCoordinates(treeWithRotatedNode$tip.label)

# Plot lines to link the tips
plotLinesBetweenTips(tipCoordinates, tipCoordinatesRotated, col=rgb(1,0,0, 0.5), lwd=2)

#############
# FUNCTIONS #
#############

plotLinesBetweenTips <- function(tipLocationsA, tipLocationsB, col="black", lwd=1, lty=1){
  
  # Prepare for adding lines across the plot panels - using grid package
  pushViewport(viewport())
  popViewport()
  
  # Examines the tips of A
  for(key in names(tipLocationsA)){
    
    # Prepare to add a single line
    pushViewport(viewport())
    
    # Check if the Y coordinate of the tips being compared has changed
    if(tipLocationsA[[key]][2] != tipLocationsB[[key]][2]){
      
      # Add a line when the Y coordinate has changed from the tip on the left tree to it on the right
      grid.lines(x = c(tipLocationsA[[key]][1], tipLocationsB[[key]][1]), 
                 y = c(tipLocationsA[[key]][2], tipLocationsB[[key]][2]), 
                 gp = gpar(col=col, lty=lty, lwd=lwd))
    }
    
    # Add the changes to the plot (the line)
    popViewport()
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
    tips[[as.character(tipLabels[i])]] <- c(grconvertX(lastPP$xx[i], "user", "ndc"), 
                                            grconvertY(lastPP$yy[i], "user", "ndc"))
  }
  
  return(tips)
}
