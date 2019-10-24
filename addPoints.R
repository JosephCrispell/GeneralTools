
calculatePointSize <- function(axisLimits){
  
  # Get the plotting window size in inches
  plotSizeInches <- par()$pin # width, height
  widthInches <- plotSizeInches[1]
  heightInches <- plotSizeInches[2]
  
  # Get the plotting window size in the plotting units
  widthX <- axisLimits[2] - axisLimits[1]
  heightY <- axisLimits[4] - axisLimits[3]
  
  # Calculate the size of a point in the current plot
  # Cex=1 is 1/72 inches (https://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/pdf.html)
  # Dividing by 72 is far too small - decided to choose 15?!?!
  pointWidth <- (widthX / widthInches) / 15
  pointHeight <- (heightY / heightInches) / 15
  
  return(c(pointWidth, pointHeight))
}

addPoints <- function(xCoords, yCoords, pch=19, col.line="black", col="black", bg="black",
                      lty=1, lwd=1, keepInside=TRUE, cex=1){
  
  #######################################################
  # Check that the input data are in the correct format #
  #######################################################
  
  # Are each of coordinate vectors the same length?
  if(length(xCoords) != length(yCoords)){
    stop("addPoints() The vectors containing the X and Y coodinates must be the same length.")
  }

  #######################
  # Get the axis limits #
  #######################
  
  # Get the axis limits
  axisLimits <- par("usr")
  
  ############################
  # Check for NA coordinates #
  ############################
  
  # Check if any NA coordinates present
  indicesOfNAs <- which(is.na(xCoords) | is.na(yCoords))
  if(length(indicesOfNAs) > 0){
    
    # Send warning
    warning("NA values present in coordinates provided. These are ignored.")
    
    # Check for each of the parameters that can have multiple parameters
    if(length(pch) == length(pch)){
      pch = pch[-indicesOfNAs]
    }
    if(length(col.line) == length(xCoords)){
      col.line = col.line[-indicesOfNAs]
    }
    if(length(col) == length(xCoords)){
      col = col[-indicesOfNAs]
    }
    if(length(lty) == length(xCoords)){
      lty = lty[-indicesOfNAs]
    }
    if(length(lwd) == length(xCoords)){
      lwd = lwd[-indicesOfNAs]
    }
    if(length(bg) == length(xCoords)){
      bg = bg[-indicesOfNAs]
    }
    if(length(cex) == length(xCoords)){
      cex = cex[-indicesOfNAs]
    }
    
    # Remove the NA coordinates
    xCoords <- xCoords[-indicesOfNAs]
    yCoords <- yCoords[-indicesOfNAs]
  }
  
  ############################
  # Check if axes are logged #
  ############################
  
  # Check X axis
  xAxisLogged <- FALSE
  if(par("xlog")){
    
    # Note that X axis was logged
    xAxisLogged <- TRUE
    
    # Log the X coordinates
    xCoords <- log10(xCoords)
    
    # Reset the X axis logged flag - fools points and polygon commands below
    par(xlog=FALSE)
  }
  
  # Check Y axis
  yAxisLogged <- FALSE
  if(par("ylog")){
    
    # Note that Y axis was logged
    yAxisLogged <- TRUE
    
    # Log the Y coordinates
    yCoords <- log10(yCoords)
    
    # Reset the Y axis logged flag - fools points and polygon commands below
    par(ylog=FALSE)
  }
  
  ###############################
  # Store the point information #
  ###############################
  
  # Calculate the height and width of point on current plot
  pointSize <- calculatePointSize(axisLimits)
  
  # Note the cex to be applied to each point
  if(length(cex) != length(xCoords)){
    cex <- rep(cex, ceiling(length(xCoords)/length(cex)))[1:length(xCoords)]
  }

  # Store the input coordinates and labels
  # !Note need to make addTextLabels have multiple cex values!
  pointInfo <- list("X"=xCoords, "Y"=yCoords, "N"=length(xCoords), "Heights"=pointSize[1]*cex, 
                    "Widths"=pointSize[2]*cex, "cex"=1)

  ###########################################
  # Produce a list of alternative locations #
  ###########################################
  
  # Generate the alternative locations
  alternativeLocations <- generateAlternativeLocations(axisLimits)
  
  # Calculate the distance between the actual and alternative points - rescale X axis remove axis range bias
  distances <- euclideanDistancesWithRescaledXAxis(pointInfo, alternativeLocations, axisLimits)
  
  ###############################################################
  # Create a list to store the information about plotted points #
  ###############################################################
  
  # Initialise the list to store the information about plotted labels
  plottedPointInfo <- list("X"=c(), "Y"=c(), "Height"=c(), "Width"=c(), "N"=0)
  
  ##############################################################
  # Add labels to plot assigning new locations where necessary #
  ##############################################################
  
  # Plot the point label
  for(i in seq_len(pointInfo$N)){
    
    # Set the colours for plotting - allows multiple colours and cycling through colours
    colour <- setOption(options=col, index=i)
    background <- setOption(options=bg, index=i)
    shape <- setOption(options=pch, index=i)
    size <- setOption(options=cex, index=i)

    # Set the line characteristics
    lineColour <- setOption(options=col.line, index=i)
    lineType <- setOption(options=lty, index=i)
    lineWidth <- setOption(options=lwd, index=i)
    
    # Get the information for the current point
    x <- pointInfo$X[i]
    y <- pointInfo$Y[i]
    height <- pointInfo$Heights[i]
    width <- pointInfo$Widths[i]

    # Get a new location
    newLocationIndex <- chooseNewLocation(pointInfo, i, alternativeLocations, distances,
                                          plottedPointInfo, axisLimits, keepInside)
    
    # Is the current point too close to others?
    if(alternativeLocations$N != 0 && newLocationIndex != -1 && 
       (tooClose(x, y, height, width, plottedPointInfo) || 
        outsidePlot(x, y, height, width, axisLimits))){
      
      # Get the coordinates for the chosen alternate location
      altX <- alternativeLocations$X[newLocationIndex]
      altY <- alternativeLocations$Y[newLocationIndex]
      
      # Add line back to previous location - from the outside of the circle
      points(x=c(altX, x), y=c(altY, y), type="l", col=col, lty=lty, lwd=lwd, xpd=TRUE)
      
      # Add point
      points(x=altX, y=altY, pch=shape, col=colour, bg=background, cex=size)

      # Append the plotted point information
      plottedPointInfo <- addPlottedLabel(x=altX, y=altY, height=height, width=width,
                                          plottedLabelInfo=plottedPointInfo)
      
      # Remove the alternative plotting location used
      alternativeLocations$X <- alternativeLocations$X[-newLocationIndex]
      alternativeLocations$Y <- alternativeLocations$Y[-newLocationIndex]
      alternativeLocations$N <- alternativeLocations$N - 1
      distances <- distances[, -newLocationIndex]
      
    }else{
      
      # Add point
      points(x=x, y=y, pch=shape, col=colour, bg=background, cex=size)
      
      # Append the plotted point information
      plottedPointInfo <- addPlottedLabel(x=x, y=y, height=height, width=width,
                                          plottedLabelInfo=plottedPointInfo)
    }
  }
  
  #####################################################################################
  # Return axes logged flags to original state - for if person makes any future plots #
  #####################################################################################
  
  par(xlog=xAxisLogged)
  par(ylog=yAxisLogged)
  
}