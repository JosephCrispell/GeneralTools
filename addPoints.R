# Load libraries
library(basicPlotteR)

# Create some random points
x <- rnorm(20)
y <- rnorm(20)

# Create an empty plot
plot(x=NULL, y=NULL, xlim=range(x), ylim=range(y), bty="n", yaxt="n", xaxt="n", xlab="", ylab="")

# Add overlapping points
#points(x, y, pch=19, col=rgb(0,0,0, 0.5), cex=5)

# Add non-overlapping points
addPoints(x, y, cex=2, col=rgb(0,0,0, 0.5), col.line="red", avoidFactor=5)


