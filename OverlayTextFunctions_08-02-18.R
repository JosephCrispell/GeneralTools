##########################
# Generate random points #
##########################

# Open a PDF
pdf(file="C:/Users/Joseph Crisp/Desktop/OverlayTextExample.pdf")

# Create some random points
n <- 50
coords <- data.frame(X=runif(n), Y=runif(n), Name="Test Label")

# Plot them without labels
plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")

# With basis labels
plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")
text(coords$X, coords$Y, labels=coords$Name, xpd=TRUE)

# Plot them with non-overlapping labels
plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")
overlayText(coords$X, coords$Y, coords$Name, cex=1, col.background=rgb(0,0,0, 0.75), col.label="white")

dev.off()

#############
# FUNCTIONS #
#############

overlayText <- function(xCoords, yCoords, labels, cex=1, col.label="red", col.line="black", col.background=NULL,
                        lty=1, lwd=1){
  
  ###########################################
  # Produce a list of alternative locations #
  ###########################################
  
  # Calculate the label heights and widths
  labelDimensions <- calculateLabelHeightsAndWidths(labels=labels, cex=cex, polygonAdded=is.null(col.background) == FALSE)
  textHeights <- labelDimensions[["Heights"]]
  textWidths <- labelDimensions[["Widths"]]
  
  # Generate the alternative locations
  output <- generateAlternativeLocations(xCoords=xCoords, yCoords=yCoords, textHeights=textHeights,
                                         textWidths=textWidths, cex=cex)
  altX <- output[["AltX"]]
  altY <- output[["AltY"]]
  
  ##############################################################
  # Add labels to plot assigning new locations where necessary #
  ##############################################################
  
  # Plot the point label
  for(i in 1:length(xCoords)){
    
    # Is the current point too close to others?
    if(tooClose(xCoords, yCoords, i, textHeights[i], textWidths[i]) == TRUE && length(altX) != 0){
      
      # Get a new location
      newLocationIndex <- chooseNewLocation(xCoords[i], yCoords[i], altX, altY)

      # Add line back to previous location
      addLineBackToOriginalLocation(altX=altX[newLocationIndex], altY=altY[newLocationIndex],
                                    x=xCoords[i], y=yCoords[i], label=labels[i], cex=cex, col=col.line, lty=lty, lwd=lwd)

      # Add label
      addLabel(x=altX[newLocationIndex], y=altY[newLocationIndex], label=labels[i], cex=cex, col=col.label, bg=col.background)

      # Remove new location and any locations too close to it
      output <- removeLocationAndThoseCloseToItFromAlternatives(
        altX, altY, newLocationIndex, textHeights[i], textWidths[i])
      altX <- output[["X"]]
      altY <- output[["Y"]]
      
    }else{
      
      # Add label
      addLabel(x=xCoords[i], y=yCoords[i], label=labels[i], cex=cex, col=col.label, bg=col.background)
    }
  }
}

addLineBackToOriginalLocation <- function(altX, altY, x, y, label, cex, col, lty=1, lwd=1){
  
  # Calculate the label width and height
  labelHeight <- strheight(label, cex=cex)
  labelWidth <- strwidth(label, cex=cex)
  
  # Create a set of points marking the boundaries of the label
  xMarkers <- c(seq(from=altX - (labelWidth * 0.52), to=altX + (labelWidth * 0.52), by=0.05*labelWidth), 
                altX + (labelWidth * 0.52))
  yMarkers <- c(seq(from=altY - (labelHeight * 0.6), to=altY + (labelHeight * 0.6), by=0.05*labelHeight),
                altY + (labelHeight * 0.6))
  
  # Calculate the closest pair of X and Y coordinates to the origin
  closestX <- xMarkers[which.min(abs(xMarkers - x))]
  closestY <- yMarkers[which.min(abs(yMarkers - y))]
  
  # Plot the line
  points(x=c(closestX, x), y=c(closestY, y), type="l", col=col, lty=lty, lwd=lwd)
}

calculateLabelHeightsAndWidths <- function(labels, cex, polygonAdded){
  
  # Get the text label heights and lengths
  textHeights <- strheight(labels) * cex
  textWidths <- strwidth(labels) * cex
  
  # If plotting boxes for labels make widths and heights slightly larger
  if(polygonAdded == TRUE){
    textHeights <- textHeights + (0.25 * textHeights)
    textWidths <- textWidths + (0.05 * textWidths)
  }
  
  # Create the output
  output <- list("Heights"=textHeights, "Widths"=textWidths)
  
  return(output)
}

generateAlternativeLocations <- function(xCoords, yCoords, textHeights, textWidths, cex=1){
  
  # Get the axis limits
  axisLimits <- par("usr")
  
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
    
    for(j in 1:length(xCoords)){
      
      if(abs(altX[i] - xCoords[j]) < textWidths[j] &&
         abs(altY[i] - yCoords[j]) < textHeights[j]){
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
  
  # Create an output
  output <- list("AltX"=altX, "AltY"=altY)
  
  return(output)
}

addLabel <- function(x, y, label, cex=1, col="black", bg=NULL, border=NA){
  
  # Calculate the height and width of the label
  labelHeight <- strheight(label, cex=cex)
  labelWidth <- strwidth(label, cex=cex)
  
  # Add a background polygon - if requested
  if(is.null(bg) == FALSE){
    polygon(x=c(x - (labelWidth * 0.52),
                x - (labelWidth * 0.52),
                x + (labelWidth * 0.52),
                x + (labelWidth * 0.52)),
            y=c(y - (labelHeight * 0.6),
                y + (labelHeight * 0.6),
                y + (labelHeight * 0.6),
                y - (labelHeight * 0.6)), 
            col=bg,
            border=border, xpd=TRUE)
  }
  
  
  # Add label
  text(x=x, y=y, labels=label, xpd=TRUE, cex=cex, col=col)
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