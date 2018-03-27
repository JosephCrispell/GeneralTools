##########################################################################
# Script to convert UK Eastings and Northings to Latitude and Longitudes #
##########################################################################

# Load the packages
library(rgdal)

# Create variables for holding the coordinate system types
# see http://www.epsg.org/
ukGrid <- "+init=epsg:27700"
latLong <- "+init=epsg:4326"

# Read in the table containing the eastings and northings
directory <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/NewCattle_16-03-18/"
file <- paste(directory, "sampleInfo_22-03-2018.csv", sep="")
table <- read.table(file, header=TRUE, sep=",")

# Get the columns that you want
spatialInfo <- subset(table, select=c("StrainId", "Mapx", "Mapy"))

# Create a unique ID for each row
spatialInfo$ID <- 1:nrow(spatialInfo)

# Create a coordinates variable
coords <- cbind(Easting = as.numeric(as.character(spatialInfo$Mapx)),
                Northing = as.numeric(as.character(spatialInfo$Mapy)))

# Create a SpatialPointsDataFrame
spatialDF <- SpatialPointsDataFrame(coords, data=data.frame(spatialInfo$StrainId,
                                                            spatialInfo$ID),
                                    proj4string = CRS(ukGrid))

# Plot the Eastings and Northings
plot(spatialDF)

# Convert the Eastings and Northings to Latitude and Longitude
spatialDFLatLongs <- spTransform(spatialDF, CRS(latLong))

# Rename columns
colnames(spatialDFLatLongs@coords)[colnames(spatialDFLatLongs@coords) == "Easting"] <- "Latitude"
colnames(spatialDFLatLongs@coords)[colnames(spatialDFLatLongs@coords) == "Northing"] <- "Longitude"

# Add the Latitudes and Longitudes onto the original table
table$Latitude <- spatialDFLatLongs@coords[,"Latitude"]
table$Longitude <- spatialDFLatLongs@coords[,"Longitude"]

# Print the table back out
file <- paste(directory, "sampleInfo_LatLongs_22-03-2018.csv", sep="")
write.csv(table, file=file, quote=FALSE, row.names=FALSE)

plot(table$Latitude, table$Longitude)
plot(table$Mapx, table$Mapy)
