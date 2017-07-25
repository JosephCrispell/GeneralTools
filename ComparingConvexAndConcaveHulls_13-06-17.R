##################
# Load Libraries #
##################

library(alphahull)
library(igraph)

###############
# Create Data #
###############

# Create random X and Ys
sample <- data.frame(x=runif(n=500, min=0, max=1), y=runif(n=500, min=0, max=1))

# plot
plot(sample, pch=19, xlab="", ylab="", xaxt="n", yaxt="n", bty="n")

# Plot the convex hull
polygon(getConvexHull(sample), border="red")

# Plot the concave hull - note alpha has to large enough to create single hull
polygon(getConcaveHull(sample, 0.05), border=rgb(0,0,1, 0.5), lwd=3)


#############
# FUNCTIONS #
#############

getConcaveHull <- function(matrix, alpha){
  
  # Remove any duplicate points
  sample <- sample[!duplicated(paste(sample$x, sample$y)), ]
  
  # Calculate the concave hull
  concaveHullInfo <- ashape(matrix, alpha=alpha)
  
  # Convert the concave hull into a graph of edges
  concaveHullGraph <- graph.edgelist(
    cbind(as.character(concaveHullInfo$edges[, "ind1"]), 
          as.character(concaveHullInfo$edges[, "ind2"])),
    directed = FALSE)
  
  # Delete an edge to make into line
  cutGraph <- concaveHullGraph - E(concaveHullGraph)[1]
  
  # Find the ends of the line
  ends = names(which(degree(cutGraph) == 1))
  
  # Connect points via shortest path
  path = get.shortest.paths(cutGraph, ends[1], ends[2])[[1]][[1]]
  
  # Convert the path to indices
  pathIndices = as.numeric(V(concaveHullGraph)[path]$name)
  
  # Join the indices
  pathIndices = c(pathIndices, pathIndices[1])
  
  # Return the polygon
  return(concaveHullInfo$x[pathIndices, ])
}

getConvexHull <- function(matrix){
  
  # Find indices of convex points
  hullPointIndices <- chull(matrix[, 1], matrix[, 2])
  
  # Join convex hull indices
  hullPointIndices <- c(hullPointIndices, hullPointIndices[1])
  
  # Plot the polygon
  return(matrix[hullPointIndices, ])
}
