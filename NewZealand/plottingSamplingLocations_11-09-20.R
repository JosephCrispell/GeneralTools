library(rgdal)
library(sp)
library(rgeos)

# Read in the sampling information
path <- "/home/josephcrispell/storage/Research/NewZealand/NewAnalyses_12-05-16/SamplingInformation"
samplingInfo <- read.csv(file.path(path, "combinedSamplingInformation_23-05-16.csv"))

# Read in a New Zealand shape file
path <- "/home/josephcrispell/Desktop/"
boundaries <- rgdal::readOGR(file.path(path, "gadm36_NZL_shp", "gadm36_NZL_0.shp"))
boundaries <- gSimplify(boundaries, tol=0.2, topologyPreserve=TRUE)



# Open a png file
png(file.path(path, "newzealand_sampling.png"), width=960, height=960)

# Set the background colour
par(bg=rgb(25,25,25, maxColorValue=255))

# Remove the plotting margins
par(mar=c(0,0,0,0))

# Plot the New Zealand boundaries
plot(boundaries, ylim=c(-51.54003, -33.20714), xlim=c(166.1375, 179.3339),
     col="black")

# Add the sampling locations
points(x=samplingInfo$Longitude, y=samplingInfo$Latitude, pch=19, cex=1,
       col=ifelse(samplingInfo$Host == "BOVINE", rgb(1,0,0,0.5), rgb(0,0,1, 0.5)))

# Close png file
dev.off()