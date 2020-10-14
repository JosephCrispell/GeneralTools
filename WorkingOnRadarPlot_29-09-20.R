# Create a skills point
# Plot a radar chart
radarChart(scores=c(4,3,2,3,3,3), 
           names=c("Programming", "Statistics", "Databases", "Projects",
                   "Web", "Versioning"),
           levels=c("Unaware","Aware","Working","Practitioner", "Expert"))
radarChart(scores=c(4,4,3,4,3,5), 
           names=c("Programming", "Statistics", "Databases", "Projects",
                   "Web", "Versioning"),
           levels=c("Unaware","Aware","Working","Practitioner", "Expert"),
           col="blue", add=TRUE)
legend(x=4.5, y=5, 
       legend=c("Current", "Aim"), 
       text.col=c(rgb(1,0,0, 0.5), rgb(0,0,1, 0.5)), 
       bty="n", text.font=2, cex=1.5, xpd=TRUE)

#### FUNCTIONS ####

radarChart <- function(scores, names, levels, col="red", alpha=0.1, 
                       axisLabelPad=1.2, circles=FALSE, add=FALSE, main="",
                       margins=c(3,3,3,3), addPoints=FALSE){
  
  # Check scores and names are the same length
  if(length(scores) != length(names)){
    stop(paste0("The number of scores (", length(scores), ") provided doesn't match the number of levels provided (", length(names), ")"))
  }
  
  # Count number of levels
  nLevels <- length(levels)
  
  # Generate the axis points (for each skill)
  axesEnds <- generateEquiDistantPointsOnCircle(length(scores), radius=nLevels)
  axesInfo <- data.frame("X"=axesEnds[, 1], "Y"=axesEnds[, 2])
  
  # Get and set the plotting margins
  currentMar <- par()$mar
  par(mar=margins)
  
  # Check not adding to existing plot
  if(add == FALSE){
    
    # Create an empty plot
    plot(x=NULL, y=NULL, xlim=c(-nLevels, nLevels), ylim=c(-nLevels, nLevels),
         bty="n", asp=1, main=main, yaxt="n", xaxt="n", xlab="", ylab="")
    
    # Add in axes titles
    text(x=axesInfo$X * axisLabelPad, y=axesInfo$Y * axisLabelPad,
         labels=names, xpd=TRUE)
    
    # Add each axis line
    for(index in seq_along(scores)){
      lines(x=c(axesInfo[index, "X"] * 1/nLevels, axesInfo[index, "X"]),
            y=c(axesInfo[index, "Y"] * 1/nLevels, axesInfo[index, "Y"]),
            lwd=0.5, col=rgb(0,0,0, 0.5))
    }
    
    # Add in levels
    for(level in 1:nLevels){
      points <- generateEquiDistantPointsOnCircle(ifelse(circles, 360, length(scores)), radius=level)
      polygon(points, border=rgb(0,0,0, 0.5), col=rgb(0,0,0, 0))
    }
    
    # Add level labels
    text(x=0, y=1:nLevels, labels=levels)
  }
  
  # Add in a skills polygon
  polygon(x=scores/nLevels * axesInfo$X,
          y=scores/nLevels * axesInfo$Y,
          border=col, col=basicPlotteR::setAlpha(col, alpha))
  if(addPoints){
    points(x=scores/nLevels * axesInfo$X,
           y=scores/nLevels * axesInfo$Y,
           pch=19, col=col)
  }
  
  # Reset plotting margins
  par(mar=currentMar)
}

generateEquiDistantPointsOnCircle <- function(nPoints, radius=1, origin=c(0,0)){
  
  # Code taken from: https://stackoverflow.com/questions/5300938/calculating-the-position-of-points-in-a-circle
  
  # Define theta
  start <- 0
  end <- (2*pi) - (2*pi)/nPoints
  theta <- seq(start, end, length.out=nPoints)
  
  # Calculate the coordinates
  x <- origin[2] + (radius * sin(theta))
  y <- origin[1] + (radius * cos(theta))
  
  # Store the coordinates in a dataframe
  return(data.frame("X"=x, "Y"=y, "Theta"=theta))
}
