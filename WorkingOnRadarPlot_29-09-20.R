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
    if(circles){
      text(x=0, y=1:nLevels, labels=levels)
    }else{
      text(x=0, y=getTickPositionsOnPolygons(axesInfo, nLevels), labels=levels)
    }
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

getTickPositionsOnPolygons <- function(axesInfo, nLevels){
  
  # Get highest Y value when X is closest to zero
  maxY <- max(axesInfo$Y[which.min(abs(axesInfo$X))])
  
  # Calculate equidistant Y values for each level
  Y <- ((1:5)/nLevels) * maxY
  
  return(Y)
}

generateEquiDistantPointsOnCircle <- function(numberOfPoints, radius=1, origin=c(0,0)){
  
  # Code adapted from: https://stackoverflow.com/questions/40279052/coordinates-of-equally-distanced-n-points-on-a-circle-in-r
  
  # Generate equidistant points on circle
  points <- sapply(1:numberOfPoints, 
                   function(position, numberOfPoints){
                     
                     # Calculate the X coordinate
                     x <- cos(2 * position * pi/numberOfPoints)
                     
                     # Calculate the Y coordinate
                     y <- sin(2 * position * pi/numberOfPoints)
                     
                     return(c(x, y))
                   }, numberOfPoints)
  
  # Transpose the matrix for plotting
  points <- t(points)
  
  # Move points based on radius
  points <- points * radius
  
  # Move points based on origin
  points[, 1] <- points[, 1] + origin[1]
  points[, 2] <- points[, 2] + origin[2]
  
  return(points)
}
