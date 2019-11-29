##################
# Load Libraries #
##################

library(alphahull)
library(igraph)

###############
# Create Data #
###############

# Create random X and Ys
theta = runif(n <- 300, 0, 2 * pi)
r = sqrt(runif(n, 0.25^2, 0.5^2))
sample = cbind(0.5 + r * cos(theta), 0.5 + r * sin(theta))
#sample <- data.frame(x=runif(n=500, min=0, max=1), y=runif(n=500, min=0, max=1))

# plot
plot(sample, pch=19, xlab="", ylab="", xaxt="n", yaxt="n", bty="n")

# Plot the convex hull
polygon(sample[chull(sample), ], border="red")

# Plot the alpha-convex hull
alphaConvexHull <- ahull(sample, alpha=0.1)
alphaConvexHullEdges <- alphaConvexHull$arcs[, c("end1", "end2")]

for(row in 1:nrow(alphaConvexHullEdges)){
  points(sample[alphaConvexHullEdges[row, ], ], col="blue", type="l")
}

# Nice example here: https://yihui.org/en/2010/04/alphahull-an-r-package-for-alpha-convex-hull/ that I could try and replicate


#############
# FUNCTIONS #
#############

getAlphaConvexHullCoords <- function