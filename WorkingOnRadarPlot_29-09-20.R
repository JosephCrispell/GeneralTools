# Create a skills point
skillsPointCurrent <- data.frame("Programming"=4, 
                                  "Statistics"=3,
                                  "Databases"=2,
                                  "Projects"=3,
                                  "Web"=3,
                                  "Versioning"=3)

# Set parameters for radar chart
min <- 1
max <- 5
nLevels <- 5

# Normalise the scores to vary between zero and 1
skillsPointCurrent[2, ] <- rescale(skillsPointsCurrent[1, ], min=min, max=max, newMin=0, newMax=1)

# Generate the axis points (for each skill)
axesEnds <- generateEquiDistantPointsOnCircle(ncol(skillsPointCurrent))
skillsPointCurrent[3, ] <- axesEnds[, 1]
skillsPointCurrent[4, ] <- axesEnds[, 2]
skillsPointCurrent[5, ] <- axesEnds[, 1] * 1.15
skillsPointCurrent[6, ] <- axesEnds[, 2] * 1.15

# Note the row names
rownames(skillsPointCurrent) <- c("Score", "NormalisedScore", "X", "Y", "LabelX", "LabelY")

# Set the plotting margins
par(mar=c(0,0,4,0))

# Create an empty plot
plot(x=NULL, y=NULL, xlim=c(-1, 1), ylim=c(-1, 1), bty="n", asp=1, main="Current")

# Add in axes titles
text(x=skillsPointCurrent["LabelX", ], y=skillsPointCurrent["LabelY", ],
     labels=colnames(skillsPointCurrent))

# Add each axis line
for(column in seq_len(ncol(skillsPointCurrent))){
  lines(x=c(0, skillsPointCurrent["X", column]),
        y=c(0, skillsPointCurrent["Y", column]),
        lwd=0.5, col=rgb(0,0,0, 0.5))
}

# Add in levels
for(level in 0:nLevels){
  
  points <- generateEquiDistantPointsOnCircle(ncol(skillsPointCurrent), 
                                              radius=level * (1/nLevels))
  polygon(points, border=rgb(0,0,0, 0.5), col=rgb(0,0,0, 0))
}

# Add in a skills polygon
polygon(x=skillsPointCurrent["NormalisedScore", ] * skillsPointCurrent["X", ],
        y=skillsPointCurrent["NormalisedScore", ] * skillsPointCurrent["Y", ],
        border=NA, col=rgb(1,0,0, 0.75))


polygon(pts.circle, col=rgb(0,0,0, 0.5), border=NA)
polygon(pts.circle*0.5, col=rgb(0,0,0, 0.5), border=NA)


#### FUNCTIONS ####

rescale <- function(values, min=min(values), max=max(values), newMin, newMax){
  
  # Rescale the values to vary between new range
  rescaledValues <- (((values - min)/(max - min)) * (newMax - newMin)) + newMin
  
  return(rescaledValues)
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

findPointOnLine <- function(destination, radius, origin=c(0,0)){
  
  # Maths taken from: https://math.stackexchange.com/questions/175896/finding-a-point-along-a-line-a-certain-distance-away-from-another-point
  
  # Store the differences between destination and origin
  differences <- c(destination[1] - origin[1], destination[2] - origin[2])
  
  # Normalise those differences 
  normalised <- differences / sqrt(sum(differences^2))
  
  # Calculate the coordinates of the point on the outer circle
  point <- origin + (radius * normalised)
  
  return(point)
}