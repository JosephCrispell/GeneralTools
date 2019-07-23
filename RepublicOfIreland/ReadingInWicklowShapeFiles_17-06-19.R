#### Load libraries ####

library(sf) # Shape files
library(rgdal) # Convert X and Y to lat longs

#### Read in the shape files ####

# Set the path
path <- "J:\\WGS_Wicklow\\"

#### Read in the badger shape file ####

# Read in the shape file
shapeFile <- paste0(path, "Raw_locations__all_species_May_2019\\badgers.shp")
badgerShapeFile <- readOGR(dsn=shapeFile)

# Extract the location data - X and Y coords
badgerLocations <- data.frame("Aliquot"=as.character(badgerShapeFile$ALIQUOT),
                              "X"=badgerShapeFile$X_COORD,
                              "Y"=badgerShapeFile$Y_COORD, stringsAsFactors=FALSE)

# Convert the X and Y coords to latitude and longitudes
latLongs <- convertXYToLatLongs(x=badgerLocations$X, y=badgerLocations$Y)
badgerLocations$Lat <- latLongs$Latitude
badgerLocations$Long <- latLongs$Longitude

#### Read in the deer shape file ####

# Read in the shape file
shapeFile <- paste0(path, "Raw_locations__all_species_May_2019\\deer.shp")
deerShapeFile <- readOGR(dsn = shapeFile)

# Store the latitude and longitude
deerLocations <- data.frame("Aliquot"=as.character(deerShapeFile$ALIQUOT), 
                            "Lat"=deerShapeFile$LAT,
                            "Long"=deerShapeFile$LONG)


#### Read in the cattle land parcel polygons ####

# Read in the shape file
shapeFile <- paste0(path, "Raw_locations__all_species_May_2019\\farms_common.shp")
cattleShapeFile <- readOGR(dsn = shapeFile)

# Extract the polygon coords
polygonCoords <- getPolygonCoords(cattleShapeFile)

# Got this from: https://stackoverflow.com/questions/29803253/r-extracting-coordinates-from-spatialpolygonsdataframe
cattleShapeFile@polygons[[10]]@Polygons[[1]]@coords
cattleShapeFile@polygons[[10]]@ID
length(cattleShapeFile@polygons[[10]]@Polygons)

#### FUNCTIONS ####

getPolygonCoords <- function(spatialDataFrame){
  
  # Got some good information from: https://stackoverflow.com/questions/29803253/r-extracting-coordinates-from-spatialpolygonsdataframe
  # Also from my blog: https://josephcrispell.github.io/BlogPosts/GetPolygonsFromShapeFile_11-10-17/GetPolygonsFromShapeFile_11-10-17.html
  # WHICH NEEDS UPDATED!!!!

  
  # Note the number of polygon sets - one for each ID
  nSets <- length(cattleShapeFile@polygons)
  
  # Initialise a list to store the IDs and coordinates of each polygon
  output <- list()
  
  # Loop through all the sets of polygons
  for(setIndex in seq_len(nSets)){
    
    # Get the ID of the current polygon set
    id <- cattleShapeFile@polygons[[setIndex]]@ID
    
    # Note the number of polygons in the current set
    nPolygons <- length(cattleShapeFile@polygons[[setIndex]]@Polygons)
    
    # Initialise a list to store the coordinates of each polygon in current set
    polygons <- list()
    
    # Examine each polygon in current set
    for(polygonIndex in seq_len(nPolygons)){
      
      # Get the coordinates for the current polygon
      coords <- cattleShapeFile@polygons[[setIndex]]@Polygons[[polygonIndex]]@coords
      colnames(coords) <- c("X", "Y")
      
      # Convert them to latitude and longitude
      latLongs <- convertXYToLatLongs(x=coords[, 1], y=coords[, 2])
      
      # Store the latitude and longitudes and X and Y coords
      polygons[[polygonIndex]] <- cbind(latLongs, coords)
    }
    
    # Store the polygon coordinates
    output[[id]] <- polygons
  }
  
  return(output)
}

convertXYToLatLongs <- function(x, y, grid="+init=epsg:29903"){
  
  # Create variables for holding the coordinate system types
  # see http://www.epsg.org/
  latLong <- "+init=epsg:4326"

  # Create a coordinates variable
  coords <- cbind(Easting = x, Northing = y)
  
  # Create a SpatialPointsDataFrame
  spatialDF <- SpatialPoints(coords, proj4string = CRS(grid))

  # Convert the Eastings and Northings to Latitude and Longitude
  spatialDFLatLongs <- spTransform(spatialDF, CRS(latLong))

  # Create a table to store the converted points
  output <- data.frame(Latitude=spatialDFLatLongs@coords[,"Easting"],
                       Longitude=spatialDFLatLongs@coords[,"Northing"])
  
  # Return the lat longs
  return(output)
}
