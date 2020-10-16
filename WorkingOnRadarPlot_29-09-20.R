# Create a skills point
# Plot a radar chart
radarChart(scores=c(4,3.5,2,3,3,3), 
           names=c("Programming", "Statistics", "Databases", "Projects",
                   "Web", "Versioning"),
           levels=c("Unaware","Aware","Working","Practitioner", "Expert"))
radarChart(scores=c(4,4,3,4,3,4), 
           names=c("Programming", "Statistics", "Databases", "Projects",
                   "Web", "Versioning"),
           levels=c("Unaware","Aware","Working","Practitioner", "Expert"),
           polygon.col="blue", add=TRUE)
legend("topright", 
       legend=c("Current", "Aim"), 
       text.col=c(rgb(1,0,0, 0.5), rgb(0,0,1, 0.5)), 
       bty="n", text.font=2, cex=1.5, xpd=TRUE)

# How I am aiming to improve:
# - Statistics - using new-to-me spatial geospatial machine learning techniques to map HIV
# - Databases - working with team in Vanuatu to explore trade statistics data with MySQL in R
# - Projects - scrum master training
# - Versioning - exploring unit testing and GitHub Actions

#### FUNCTIONS ####

radarChart <- function(scores, names, levels, 
                       polygon.col="red", polygon.alpha=0.1, polygon.pch=19,
                       axisLabelPad=1.2, circles=FALSE, add=FALSE, main="",
                       margins=c(3,3,3,3), addPoints=FALSE, 
                       radar.col=rgb(0,0,0, 0.5), radar.lty=1, radar.lwd=0.5,
                       levels.font=1, levels.cex=1, labels.font=2, labels.cex=1){
  
  # Check scores and names are the same length
  if(length(scores) != length(names)){
    stop(paste0("The number of scores (", length(scores), ") provided doesn't match the number of names provided (", length(names), ")"))
  }
  
  # Check scores aren't outside levels
  levelsRange <- 1:length(levels)
  if(sum(scores < 1) > 0 || sum(scores > length(levels)) > 0){
    stop("The scores provided don't fall on or within the specified levels.")
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
         labels=names, xpd=TRUE, font=labels.font, cex=labels.cex)
    
    # Add each axis line
    for(index in seq_along(scores)){
      lines(x=c(axesInfo[index, "X"] * 1/nLevels, axesInfo[index, "X"]),
            y=c(axesInfo[index, "Y"] * 1/nLevels, axesInfo[index, "Y"]),
            lwd=radar.lwd, lty=radar.lty, col=radar.col)
    }
    
    # Add in levels
    for(level in 1:nLevels){
      points <- generateEquiDistantPointsOnCircle(ifelse(circles, 360, length(scores)), radius=level)
      polygon(points, border=radar.col, col=rgb(0,0,0, 0), lwd=radar.lwd, lty=radar.lty)
    }
    
    # Add level labels
    text(x=0, y=1:nLevels, labels=levels, font=levels.font, cex=levels.cex)
  }
  
  # Add in a skills polygon
  polygon(x=scores/nLevels * axesInfo$X,
          y=scores/nLevels * axesInfo$Y,
          border=polygon.col, col=basicPlotteR::setAlpha(polygon.col, polygon.alpha))
  if(addPoints){
    points(x=scores/nLevels * axesInfo$X,
           y=scores/nLevels * axesInfo$Y,
           pch=polygon.pch, col=polygon.col)
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
