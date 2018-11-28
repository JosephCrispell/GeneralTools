#### Load libraries ####

library(contoureR)

#### Load data ####

# Read in the data
file <- "/home/josephcrispell/Desktop/Research/LiliAnalyses/MI_Elk_Data_134isolates_traits_withClades_v3.csv"
table <- read.table(file, header=TRUE, sep=",")

#### Examine spatial clustering ####

# Open the pdf
pdf("/home/josephcrispell/Desktop/Research/LiliAnalyses/ExaminingSpatialClustering_20-11-18.pdf")

# Set the plotting dimensions
par(mfrow=c(2,2))

# Define colours for the clusters 
colours <- c(rgb(1,0,0, 0.75), rgb(0,1,0, 0.75), rgb(0,0,1, 0.75), rgb(0,0,0, 0.75))
table$Colour <- colours[table$Clade]

# Plot the locations and colour by clade
plotDataWithPolygons(table, colours)

# Compare random samples of the data with real data and permuted data
differences <- calculateDistancesBetweenRandomSamples(table, nComparison=1000)
differencesPermutedClades <- calculateDistancesBetweenRandomSamples(table, nComparison=1000, permute=TRUE)
plotDifferences(differences, differencesPermutedClades)

# Kmeans clustering
runKMeansClustering(table, colours)

# Reset the plotting window dimensions
par(mfrow=c(1,1))

# Principle components analysis
runPCA(table, colours, scale=TRUE)

# Run a Logistic regression
# POINT_X:POINT_Y - adds an interaction term
results <- glm(as.factor(Clade) ~ POINT_X + POINT_Y + POINT_X:POINT_Y  + table$Year, data=table, family="binomial")
summary(results)

# Close the PDF
dev.off()

#### FUNCTIONS ####

runPCA <- function(table, colours, scale=TRUE){
  
  # Run the principle components analysis
  pcaResults <- prcomp(table[, c("Year", "POINT_X", "POINT_Y")], center = TRUE, scale. = scale)
  
  # Print summary
  print(pcaResults)
  
  # Plot the importance of each component
  plot(pcaResults, type="l", main="Variation explained by Principle Components")
  
  ## Plot the components against one another
  # 1 vs. 2
  plot(x=pcaResults$x[, "PC1"], y=pcaResults$x[, "PC2"], pch=19, col=table$Colour, cex=2,
       xlab="Principle Component 1", ylab="Principle Component 2", las=1, bty="n",
       cex.axis=0.75)
  
  # Add a legend
  legend("topleft", legend=c("Clade: ", unique(table$Clade)), text.col=c("black", colours[unique(table$Clade)]), bty="n")
  
  # 1 vs. 3
  plot(x=pcaResults$x[, "PC1"], y=pcaResults$x[, "PC3"], pch=19, col=table$Colour, cex=2,
       xlab="Principle Component 1", ylab="Principle Component 3", las=1, bty="n",
       cex.axis=0.75)
  
  # Add a legend
  legend("topleft", legend=c("Clade: ", unique(table$Clade)), text.col=c("black", colours[unique(table$Clade)]), bty="n")
  
  # 2 vs. 3
  plot(x=pcaResults$x[, "PC2"], y=pcaResults$x[, "PC3"], pch=19, col=table$Colour, cex=2,
       xlab="Principle Component 2", ylab="Principle Component 3", las=1, bty="n",
       cex.axis=0.75)
  
  # Add a legend
  legend("topleft", legend=c("Clade: ", unique(table$Clade)), text.col=c("black", colours[unique(table$Clade)]), bty="n")
}

plotDataWithPolygons <- function(table, colours){
  
  # Plot the points
  plot(table$POINT_X, table$POINT_Y, col=table$Colour, pch=19, bty="n", xlab="X", ylab="Y", 
       main="", las=1)
  
  # Add a polygon for each clade
  for(clade in unique(table$Clade)){
    
    # Get the subset of the data for the current clade
    subset <- table[table$Clade == clade, ]
    
    # Calculate the convex hull
    hull <- getConvexHull(x=subset$POINT_X, y=subset$POINT_Y)
    
    # Plot the convex hull
    polygon(x=subset[hull, "POINT_X"], y=subset[hull, "POINT_Y"], border=colours[clade])
  }
  
  # Add a legend
  legend("topright", legend=c("Clade:", unique(table$Clade)), text.col=c("black", colours[unique(table$Clade)]), bty="n")
  
  # Get the axis limits
  axisLimits <- par("usr")
  
  # Add a plot label
  mtext("A", side=3, line=1, at=axisLimits[1], cex=2.5)
}

runKMeansClustering <- function(table, colours){
  # Code mostly taken from: http://gsp.humboldt.edu/OLM/R/03_02_ClusterAnalysis.html
  
  # Create a matrix to store the coordinate data
  matrix <- as.matrix(data.frame(x=table$POINT_X, y=table$POINT_Y, clade=table$Clade))
  
  # Run a kmeans analysis for 4 clusters - the number of clades we have
  clusterInfo <- kmeans(matrix[, c("x", "y")], centers=3)
  
  # Plot the kmeans clusters and 
  plot(matrix, col=colours[matrix[, "clade"]], pch=c(21, 22, 23, 24)[clusterInfo$cluster], 
       bg=colours[matrix[, "clade"]], bty="n", las=1)
  
  # Add legend
  legend("topright", legend=c("Clade", unique(table$Clade)), 
         text.col=c("black", colours[unique(table$Clade)]), bty="n")
  legend("top", legend=c("K-means cluster", sort(unique(clusterInfo$cluster))), 
         pch=c(19, c(21, 22, 23, 24)[sort(unique(clusterInfo$cluster))]), bty="n",
         col=c("white", "black", "black", "black", "black"))
  
  # Get the axis limits
  axisLimits <- par("usr")
  
  # Add a plot label
  mtext("C", side=3, line=1, at=axisLimits[1], cex=2.5)
  
  # Calculate the variance of each column
  variance <- apply(matrix, 2, var)
  
  # Find the sum of squares for 1 cluster
  withinClusterSumOfSquares = (nrow(matrix)-1)*sum(variance)
  
  # Find the sum of squares for 2 to 15 clusters
  for (i in 2:30) {
    clusterInfo <- kmeans(matrix, centers=i)
    withinClusterSumOfSquares[i] <- sum(clusterInfo$withinss)
  }
  
  # Plot the result
  plot(1:30, withinClusterSumOfSquares, type="o", xlab="Number of Clusters", ylab="Within groups sum of squares",
       bty="n", las=1)
  
  # Get the axis limits
  axisLimits <- par("usr")
  
  # Add a plot label
  mtext("D", side=3, line=1, at=axisLimits[1], cex=2.5)
  
}

plotDifferences <- function(differences, differencesPermutedClades){
  
  # Get the Y limits
  yLim <- c(0, max(c(max(differences), max(differencesPermutedClades))))
  
  # Build an initial empty plot
  plot(x=NULL, y=NULL, xlim=c(0.75, 3.5), ylim=yLim, bty="n", xaxt="n", xlab="",
       ylab="Spatial difference (m)", las=1,
       main="")
  
  # Add an X axis
  axis(side=1, at=1:3, labels=colnames(differences))
  
  # Plot the 95% bounds of the distributions before and after permutation
  for(i in 1:3){
    
    # Get the column
    column <- colnames(differences)[i]
    
    # Plot the points from BEFORE permutating the clades
    bounds <- quantile(differences[, column], probs=c(0.025, 0.975))
    points(x=c(i-0.15, i-0.15), y=c(bounds[1], bounds[2]), type="l", col="red", lwd=3)
    points(x=c(i-0.1, i-0.2), y=c(bounds[1], bounds[1]), type="l", col="red", lwd=3)
    points(x=c(i-0.1, i-0.2), y=c(bounds[2], bounds[2]), type="l", col="red", lwd=3)
    points(x=i-0.15, y=median(differences[, column]), pch=19)
    
    # Plot the points from AFTER permutating the clades
    bounds <- quantile(differencesPermutedClades[, column], probs=c(0.025, 0.975))
    points(x=c(i+0.15, i+0.15), y=c(bounds[1], bounds[2]), type="l", col="blue", lwd=3)
    points(x=c(i+0.1, i+0.2), y=c(bounds[1], bounds[1]), type="l", col="blue", lwd=3)
    points(x=c(i+0.1, i+0.2), y=c(bounds[2], bounds[2]), type="l", col="blue", lwd=3)
    points(x=i+0.15, y=median(differencesPermutedClades[, column]), pch=19)
  }
  
  # Get the axis limits
  axisLimits <- par("usr")
  
  # Add a plot label
  mtext("B", side=3, line=1, at=axisLimits[1], cex=2.5)
  
  # Add a legend
  legend("top", legend=c("Before permutation", "After permutation"), text.col=c("red", "blue"), bty="n")
}

calculateDistancesBetweenRandomSamples <- function(table, nComparisons=100, permute=FALSE){
  
  # Permute the Clade IDs if requested
  if(permute){
    table$Clade <- sample(table$Clade)
  }
  
  # Initialise a data frame to store the distances between random samples from each clade
  differences <- data.frame(Diff_1_2=NA, Diff_1_3=NA, Diff_2_3=NA)
  
  # Compare each clade once
  for(i in unique(table$Clade)){
    for(j in unique(table$Clade)){
      
      # Skipping when i >= j means all clades will only be compared once
      if(i >= j){
        next
      }
      
      # Subset the data for Clade i and j
      subsetI <- table[table$Clade == i, ]
      subsetJ <- table[table$Clade == j, ]
      
      # Randomly pick Clade i coordinates
      randomI <- sample(seq_len(nrow(subsetI)), replace=TRUE, size=nComparisons)
      
      # Randomly pick Clade j coordinates
      randomJ <- sample(seq_len(nrow(subsetJ)), replace=TRUE, size=nComparisons)
      
      # Calculate spatial distance between random coordinates
      for(n in seq_len(nComparisons)){
        
        differences[n, paste0("Diff_", i, "_", j)] <- 
          euclideanDistance(x1=subsetI[randomI[n], "POINT_X"], y1=subsetI[randomI[n], "POINT_Y"],
                            x2=subsetJ[randomJ[n], "POINT_X"], y2=subsetJ[randomJ[n], "POINT_Y"])
      }
    }
  }
  
  return(differences)
}

euclideanDistance <- function(x1, y1, x2, y2){
  return(sqrt(sum((x1 - x2)^2 + (y1 - y2)^2)))
}


