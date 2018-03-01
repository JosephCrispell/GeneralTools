##########################
# Generate random points #
##########################

# Create some random points
n <- 25
x <- runif(n)
y <- runif(n)
labels <- rep("Test Label", n)

# Plot them
plot(x=x, y=y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")

###############
# Plot labels #
###############

overlayText(x, y, labels, cex=1)


#############
# FUNCTIONS #
#############

overlayText <- function(x, y, labels, cex){
  
  ###########################################
  # Produce a list of alternative locations #
  ###########################################
  
  # Get the axis limits
  axisLimits <- par("usr")
  
  # Get the text label heights and lengths
  textHeights <- strheight(labels) * cex
  textWidths <- strwidth(labels) * cex
  
  # Define the spacer for each axis
  spacerX <- 0.01 * (axisLimits[2] - axisLimits[1])
  spacerY <- 0.01 * (axisLimits[4] - axisLimits[3])
  
  # Generate the set of points based upon the spacer
  altX <- c()
  altY <- c()
  for(i in seq(axisLimits[1], axisLimits[2], spacerX)){
    for(j in seq(axisLimits[3], axisLimits[4], spacerY)){
      
      altX[length(altX) + 1] <- i
      altY[length(altY) + 1] <- j
    }
  }
  #points(altX, altY, col=rgb(0,0,0, 0.5), pch=20, xpd=TRUE)
  
  # Remove points that are too close to actual values
  remove <- c()
  for(i in 1:length(altX)){
    
    for(j in 1:length(x)){
      
      if(abs(altX[i] - x[j]) < textWidths[j] &&
         abs(altY[i] - y[j]) < textHeights[j]){
        remove[length(remove) + 1] <- i
        break
      }
    }
  }
  #points(altX[remove], altY[remove], col=rgb(1,1,1), pch=20, xpd=TRUE)
  if(length(remove) > 0){
    altX <- altX[-remove]
    altY <- altY[-remove]
  }
  
  ##############################################################
  # Add labels to plot assigning new locations where necessary #
  ##############################################################
  
  # Plot the point label
  for(i in 1:length(x)){
    
    # Is the current point too close to others?
    if(tooClose(x, y, i, textHeights[i], textWidths[i]) == TRUE && length(altX) != 0){
      
      # Get a new location
      newLocationIndex <- chooseNewLocation(x[i], y[i], altX, altY)
      
      # Add label
      text(x=altX[newLocationIndex], y=altY[newLocationIndex],
           labels=labels[i], xpd=TRUE, cex=cex)
      
      # Add line back to previous location
      points(x=c(altX[newLocationIndex], x[i]),
             y=c(altY[newLocationIndex], y[i]),
             type="l", col=rgb(0,0,0, 0.5))
      
      # Remove new location and any locations too close to it
      output <- removeLocationAndThoseCloseToItFromAlternatives(
        altX, altY, newLocationIndex, textHeights[i], textWidths[i])
      altX <- output[["X"]]
      altY <- output[["Y"]]
      
    }else{
      text(x=x[i], y=y[i],
           labels=labels[i], xpd=TRUE, cex=cex)
    }
  }
}

removeLocationAndThoseCloseToItFromAlternatives <- function(altX, altY, index, textHeight, textWidth){
  remove <- c(index)
  for(i in 1:length(altX)){
    
    if(i == index){
      next
    }
    
    if(abs(altX[index] - altX[i]) < textWidth &&
       abs(altY[index] - altY[i]) < textHeight){
      remove[length(remove) + 1] <- i
    }
  }
  
  altX <- altX[-remove]
  altY <- altY[-remove]
  
  return(list("X" = altX, "Y" = altY))
}

chooseNewLocation <- function(x, y, altXs, altYs){
  
  # Calculate the distance from point to all alternatives
  distances <- c()
  for(i in 1:length(altXs)){
    distances[i] <- euclideanDistance(x, y, altXs[i], altYs[i])
  }
  
  return(which.min(distances))
}

tooClose <- function(x, y, index, textHeight, textWidth){
  
  result <- FALSE
  for(i in 1:length(x)){
    
    if(i == index){
      next
    }else if(abs(x[index] - x[i]) < textWidth &&
             abs(y[index] - y[i]) < textHeight){
      result <- TRUE
      break
    }
  }
 
  return(result) 
}

euclideanDistance <- function(x1, y1, x2, y2){
  return(sqrt(sum((x1 - x2)^2 + (y1 - y2)^2)))
}