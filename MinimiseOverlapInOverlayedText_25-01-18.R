
# Create some random points
n <- 25
coords <- data.frame(X=runif(n), Y=runif(n), Name="Test Label")

# Plot them
plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red")

# Plot text labels
overlayText(coords$X, coords$Y, coords$Name, cex=0.5)

#############
# FUNCTIONS #
#############

overlayText <- function(x, y, labels, xThresholdProp=0.05, yThresholdProp=0.05,
                        pos=NULL, cex=1){
  
  # Define your thresholds for closeness
  axisLimits <- par("usr")
  xThreshold <- xThresholdProp * (axisLimits[2] - axisLimits[1])
  yThreshold <- yThresholdProp * (axisLimits[4] - axisLimits[3])
  
  # Put coordinates and labels into dataframe
  coords <- data.frame(X=x, Y=y, Name=labels)
  
  # Get alternative points for when labels overlap
  alternatives <- getAlternativeLabelLocations(n=1000, xThreshold, yThreshold, coords,
                                               axisLimits)
  
  # Assign new locations where necessary
  newLocations <- getNewLocations(coords, alternatives, xThreshold, yThreshold)
  
  # Plot labels
  plotTextLabels(newLocations, coords, pos, cex)
}

plotTextLabels <- function(newLocations, oldLocations, pos, cex){
  
  for(i in 1:nrow(newLocations)){
    
    if(newLocations[i, "X"] != oldLocations[i, "X"] || 
       newLocations[i, "Y"] != oldLocations[i, "Y"]){
      
      text(x=newLocations[i, "X"], y=newLocations[i, "Y"],
           labels=newLocations[i, "Name"], xpd=TRUE, pos=pos, cex=cex)
      points(x=c(newLocations[i, "X"], oldLocations[i, "X"]),
             y=c(newLocations[i, "Y"], oldLocations[i, "Y"]),
             type="l", col=rgb(0,0,0, 0.5))
    }else{
      text(x=newLocations[i, "X"], y=newLocations[i, "Y"],
           labels=newLocations[i, "Name"], xpd=TRUE, pos=pos, cex=cex)
    }
  }
}

getNewLocations <- function(coords, alternatives, xThreshold, yThreshold){
  
  # Assign new locations to labels if necessary
  newLocations <- coords
  skip <- c()
  for(i in 1:nrow(coords)){
    
    # Check if current point too close to others
    for(j in 1:nrow(coords)){
      
      # Skip if same location
      if(i == j || j %in% skip){
        next
      }
      
      # Check if current locations are too close
      if(abs(coords[i, "X"] - coords[j, "X"]) < xThreshold &&
         abs(coords[i, "Y"] - coords[j, "Y"]) < yThreshold){
        
        newLocationIndex <- getNewLocation(coords[i, "X"], coords[i, "Y"], alternatives)
        newLocations[i, c(1, 2)] <- alternatives[newLocationIndex, ]
        alternatives <- alternatives[-newLocationIndex, ]
        skip[length(skip) + 1] <- i
        
        break
      }
    }
    
    # Check if ran out of new locations
    if(is.na(alternatives[1, 1]) == TRUE){
      break
    }
  }
  
  return(newLocations)
}

getNewLocation <- function(x, y, alternatives){
  
  # Calculate the distance to all alternative locations
  distances <- c()
  for(i in 1:nrow(alternatives)){
    distances[i] <- euclideanDistance(x, y, alternatives[i, "X"], alternatives[i, "Y"])
  }
  
  # Select the closest
  chosen <- which(distances == min(distances))[1]
  
  return(chosen)
}

getAlternativeLabelLocations <- function(n, xThreshold, yThreshold, coords, axisLimits){
  
  # Get a set of random points
  randomX <- runif(n, min=axisLimits[1], max=axisLimits[2])
  randomY <- runif(n, min=axisLimits[3], max=axisLimits[4])
  
  # Remove points that are too close to actual values
  remove <- c()
  for(i in 1:n){
    
    for(j in 1:nrow(coords)){
      
      if(abs(randomX[i] - coords[j, "X"]) < xThreshold &&
         abs(randomY[i] - coords[j, "Y"]) < yThreshold){
        remove[length(remove) + 1] <- i
        break
      }
    }
  }
  randomX <- randomX[-remove]
  randomY <- randomY[-remove]
  
  # Remove points that are too close to each other
  skip <- c()
  for(i in 1:length(randomX)){
    
    for(j in 1:length(randomX)){
      
      if(j %in% skip || i == j){
        next
      }
      
      if(abs(randomX[i] - randomX[j]) < xThreshold &&
         abs(randomY[i] - randomY[j]) < yThreshold){
        skip[length(skip) + 1] <- i
        break
      }
    }
  }
  randomX <- randomX[-skip]
  randomY <- randomY[-skip]
  
  alternatives <- data.frame(X=randomX, Y=randomY)
  return(alternatives)
}




