

par(mar=c(0,0,0,0))
plot(x=NULL, y=NULL, xlim=c(-1, 1), ylim=c(-1, 1), bty="n", asp=1)

points(x=0, y=0, pch=19, cex=2)

n <- 10 # number of points you want on the unit circle
#pts.circle <- t(sapply(1:n,function(r)c(cos(2*r*pi/n),sin(2*r*pi/n))))
pts.circle <- generateEquiDistantPointsOnCircle(n)
points(pts.circle, col='red', pch=19, xlab='x', ylab='y')
points(pts.circle*0.5, col='blue', pch=19, xlab='x', ylab='y')

polygon(pts.circle, col=rgb(0,0,0, 0.5), border=NA)
polygon(pts.circle*0.5, col=rgb(0,0,0, 0.5), border=NA)


#### FUNCTIONS ####

generateEquiDistantPointsOnCircle <- function(numberOfPoints, radius=1, origin=c(0,0)){
  
  # Code adapted from: https://stackoverflow.com/questions/40279052/coordinates-of-equally-distanced-n-points-on-a-circle-in-r
  
  # Generate equidistant points on circle
  points <- sapply(1:n, 
                   function(position, n){
                     
                     # Calculate the X coordinate
                     x <- cos(2 * position * pi/n)
                     
                     # Calculate the Y coordinate
                     y <- sin(2 * position * pi/n)
                     
                     return(c(x, y))
                   }, n)
  
  # Transpose the matrix for plotting
  points <- t(points)
  
  # Move points based on radius
  points <- points * radius
  
  # Move points based on origin
  points[, 1] <- points[, 1] + origin[0]
  points[, 2] <- points[, 2] + origin[1]
  
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