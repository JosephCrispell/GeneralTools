##################
# Load Libraries #
##################

library(maptools) # Read shape file
library(rgeos) # Polygon centroids
library(ape) # Phylogenetic tree
library(geiger) # Phylogenetic tree tools
library(plotrix) # Draw circle

# Create a path variable
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/"

#######################################################################
# Calculate the Territory Centroids of Each Social Group in Each Year #
#######################################################################

# Note the years that territory shape files are available for
years <- c(2000:2011)
shapeFileNames <- c("territories_2000.shp", 
                    "territories_2001_autumn.shp",
                    "territories_2002_spring.shp",
                    "territories_2003.shp",
                    "territories_2004.shp",
                    "territories_2005.shp",
                    "MCPs2006 corrected.shp",
                    "MCPs_07.shp",
                    "MCP_2008.shp",
                    "MCP_2009.shp",
                    "MCPs 2010.shp",
                    "MCP_2011.shp")

# Initialise list to record each social group's location in each year
groupsCentroidsPerYear <- initialiseTableToStoreTerritoryCentroids(years)
groupsRows <- list()
index <- 0

# Initialise a list to record the territories from each year
territoriesForEachYear <- list()

# Initialise two arrays to store the min and max X and Y values for the territories
min <- c(999999999999, 99999999999)
max <- c(0, 0)

# Examine the territories in each year
for(i in 1:length(years)){
  
  year <- years[i]
  shapeFileName <- shapeFileNames[i]
  
  # Read in the shape file
  file <- paste(path, "BadgerTerritoryMarkingData/",
                "Baitmarking ", year, "/", shapeFileName, sep="")
  territories <- readShapePoly(file) # Generates SpatialPolygonsDataFrame
  
  # Extract the polygon coordinates
  territoryCoords <- getPolygonCoords(territories)
  min <- updateMin(territoryCoords, min)
  max <- updateMax(territoryCoords, max)
  
  # Get the full social group names
  territoryIDs <- getSocialGroupNames(territories@data, year)
  
  # Assign polygons to their social group names
  socialGroupTerritories <- assignTerritoriesToSocialGroupNames(territoryIDs, territoryCoords)
  territoriesForEachYear[[as.character(year)]] <- socialGroupTerritories
  
  # Calculate the territory centroids - mean X and Y
  territoryCentroids <- calculateTerritoryCentroids(territoryCoords, territoryIDs)
  
    # Give an index to all groups present
  groupsRows <- indexGroups(names(territoryCentroids), groupsRows, index)
  
  # Note the centroids of each group in the current year
  groupsCentroidsPerYear <- addGroupCentroidsFromCurrentYear(groupsCentroidsPerYear, groupsRows,
                                                             territoryCentroids, year)
  index <- nrow(groupsCentroidsPerYear)
}

# Add social group names to territory centroids table
groupsCentroidsPerYear <- addSocialGroupNames(groupsCentroidsPerYear, groupsRows)

##########################################
# Get badger isolates sampling locations #
##########################################

# Read in the newick tree
file <- paste(path, "allVCFs-IncludingPoor/vcfFiles/",
              "mlTree_Prox-10_plusRef_rmResequenced_SNPCov-0.1_28-10-16.tree", sep="")
tree <- read.tree(file=file)

# Get a list of the isolates in the clade
node <- 289
cladeTips <- tips(tree, node=node)

# Read in the isolate metadata
fileName <- paste(path, "IsolateData/", "BadgerInfo_08-04-15_LatLongs_XY.csv",
                  sep="")
metadata <- read.table(fileName, header=TRUE, stringsAsFactors=FALSE, sep=",")

# Keep only rows for the isolates
metadata <- metadata[metadata$WB_id %in% cladeTips, ]

##########################################
# Get cattle isolates sampling locations #
##########################################

# Cattle Isolates
file <- paste(path, "IsolateData/", 
              "CattleIsolateInfo_LatLongs_plusID_outbreakSize_Coverage_AddedTB1453-TB1456.csv", sep="")
cattleInfo <- read.table(file, header=TRUE, sep=",")

# Select info only for cattle isolates in clade
cattleInfo <- cattleInfo[cattleInfo$StrainId %in% cladeTips, ]

##############################################
# Plot the badger territories from all years #
##############################################

# Open a PDF
file <- paste(path, "BASTA/", "DemeAssignmentDiagram_24-05-17.pdf", sep="")
pdf(file)

# Find centre point of badger locations
badgerCentre <- findBadgerGroupsCentre(groupsCentroidsPerYear)
expand <- 7000

# Plot the badger territories from each year onto a single plot
alpha <- 0.5
plotTerritories(territoriesForEachYear, badgerCentre, expand, alpha)

# Add legend
legend(x=386000, y=200000, legend=c("BADGER", "COW"), text.col=c("red", "blue"), bty="n")

#################################
# Plot the sampled cattle herds #
#################################

for(row in 1:nrow(cattleInfo)){
  points(x=cattleInfo[row, "Mapx"], y=cattleInfo[row, "Mapy"], pch=20,
         col=rgb(0,0,1, 0.5), cex=2)
}

#########################################
# Add Badger sampling locations to plot #
#########################################

# Plot the badger isolate sampling locations
addBadgerIsolatesLocations(metadata=metadata, groupsRows=groupsRows, 
                           groupsCentroidsPerYear=groupsCentroidsPerYear,
                           col=rgb(1,0,0, 0.5), cex=2)

###############################################
# Add circles to define inner and outer demes #
###############################################

# Define inner circle radius
thresholdDistance <- 3000
outerDistance <- 7300

# Add deme defining circles
draw.circle(x=badgerCentre[1], y=badgerCentre[2], radius=thresholdDistance,
            border="black")
draw.circle(x=badgerCentre[1], y=badgerCentre[2], radius=outerDistance,
            border="black")
text(x=381938.4, y=198700, labels=paste(paste(thresholdDistance / 1000, "km")))
text(x=381938.4, y=194450, labels=paste(paste(outerDistance / 1000, "km")))

dev.off()


#############
# FUNCTIONS #
#############

findBadgerGroupsCentre <- function(groupsCentroidsPerYear){
  
  Xs <- c()
  Ys <- c()
  index <- 0
  
  for(row in 1:nrow(groupsCentroidsPerYear)){
    
    for(col in 2:ncol(groupsCentroidsPerYear)){
      
      if(is.na(groupsCentroidsPerYear[row, col]) == FALSE){
        index <- index + 1
        
        values <- as.numeric(
          strsplit(groupsCentroidsPerYear[row, col], split=":")[[1]])
        Xs[index] <- values[1]
        Ys[index] <- values[2]
      }
    }
  }
  
  means <- c(mean(Xs), mean(Ys))
  
  return(means)
}

addBadgerIsolatesLocations <- function(metadata, groupsRows, groupsCentroidsPerYear,
                                       col, cex){
  # Plot the isolate sampling locations (social group centroids)
  for(row in 1:nrow(metadata)){
    
    year <- strsplit(metadata[row, "date"], split="/")[[1]][3]
    group <- paste(strsplit(metadata[row, "Social.Group.Trapped.At"], split=" ")[[1]],
                   collapse="")
    
    territoryCentroid <- strsplit(groupsCentroidsPerYear[groupsRows[[group]], 
                                                         as.character(year)],
                                  split=":")
    
    if(length(territoryCentroid) > 0){
      territoryCentroid <- as.numeric(territoryCentroid[[1]])
      points(x=territoryCentroid[1], y=territoryCentroid[2], pch=20,
             col=col, cex=cex)
    }
  }
}

getSequenceIDsFromFastaFile <- function(fileName){
  
  # Get the lines from the FASTA file
  connection <- file(fileName, open='r')
  fileLines <- readLines(connection)
  close(connection)
  
  # Initialise a vector to store the IDs
  ids <- c()
  index <- 0
  
  # Read the fasta file line by line
  for(line in fileLines[-1]){
    
    # Skip sequence lines
    if(grepl(x=line, pattern=">") == FALSE){
      next
    }
    
    index <- index + 1
    ids[index] <- strsplit(substr(x=line, start=2, stop=nchar(line)), split="_")[[1]][1]
  }
  
  return(ids)
}

updateMax <- function(territoryCoords, max){
  
  newValues <- territoryCoords[["max"]]
  
  if(newValues[1] > max[1]){
    max[1] <- newValues[1]
  }
  
  if(newValues[2] > max[2]){
    max[2] <- newValues[2]
  }
  
  return(max)
}

updateMin <- function(territoryCoords, min){
  
  newValues <- territoryCoords[["min"]]
  
  if(newValues[1] < min[1]){
    min[1] <- newValues[1]
  }
  
  if(newValues[2] < min[2]){
    min[2] <- newValues[2]
  }
  
  return(min)
}

plotTerritories <- function(territoriesForEachYear, centre, expand,
                            alpha){
  
  par(mar=c(0,0,0,0))
  plot(x=NULL, y=NULL, yaxt="n", xaxt="n", bty="n", ylab="",
       xlim=c(centre[1] - expand, centre[1] + expand), 
       ylim=c(centre[2] - expand, centre[2] + expand), asp=1,
       xlab="")
  
  for(year in names(territoriesForEachYear)){
    
    territoryCoords <- territoriesForEachYear[[year]]
    for(id in names(territoryCoords)){
      
      polygon(territoryCoords[[id]], border=rgb(0,0,0, alpha))
    }
  }
}


addSocialGroupNames <- function(groupsCentroidsPerYear, groupsRows){
  
  for(group in names(groupsRows)){
    print(group)
    print(groupsRows[[group]])
    groupsCentroidsPerYear[groupsRows[[group]], "SocialGroup"] <- group
  }
  
  return(groupsCentroidsPerYear)
}

indexGroups <- function(groups, groupsRows, index){
  
  for(i in 1:length(groups)){
    
    if(is.null(groupsRows[[groups[i]]]) == TRUE){
      index <- index + 1
      groupsRows[[groups[i]]] <- index
    }
  }
  
  return(groupsRows)
}

addGroupCentroidsFromCurrentYear <- function(groupsCentroidsPerYear, groupsRows, territoryCentroids, year){
  
  # Update the territory centroids for the current year
  for(group in names(territoryCentroids)){
    
    coords <- territoryCentroids[[group]]
    
    groupsCentroidsPerYear[groupsRows[[group]], as.character(year)] <- paste(coords[1], ":", coords[2], sep="")
  }
  
  return(groupsCentroidsPerYear)
}

initialiseTableToStoreTerritoryCentroids <- function(years){
  groupsCentroidsPerYear <- data.frame(SocialGroup=c(NA), stringsAsFactors=FALSE)
  for(year in years){
    groupsCentroidsPerYear[, as.character(year)] <- NA
  }
  
  return(groupsCentroidsPerYear)
}

calculateTerritoryCentroids <- function(territoryCoords,
                                        territoryIDs){
  
  territoryCentroids <- list()
  polygonIDs <- names(territoryCoords)
  for(id in polygonIDs){
    
    if(id == "min" | id == "max"){
      next
    }
    
    territoryCentroids[[toupper(territoryIDs[[id]])]] <- c(
      mean(territoryCoords[[id]][, 1]),
      mean(territoryCoords[[id]][, 2]))
  }
  
  return(territoryCentroids)
}

assignTerritoriesToSocialGroupNames <- function(territoryIDs, territoryCoords){
  socialGroupTerritories <- list()
  for(id in names(territoryIDs)){
    socialGroupTerritories[[toupper(territoryIDs[[id]])]] <- territoryCoords[[id]]
  }
  
  return(socialGroupTerritories)
}

getSocialGroupNames <- function(polygonInfo, year){
  territoryIDs <- list()
  rowNames <- rownames(polygonInfo)
  
  column <- 1
  if(year == 2007){
    column <- 2
  }
  
  for(row in 1:nrow(polygonInfo)){
    socialGroupName <- as.character(polygonInfo[row, column])
    territoryIDs[[rowNames[row]]] <- removeSep(socialGroupName, sep=" ")
  }
  
  return(territoryIDs)
}

removeSep <- function(string, sep){
  parts <- strsplit(x=string, split=sep)[[1]]
  return(paste(parts, collapse=""))
}

getPolygonCoords <- function(spatialPolygonsDataFrame){
  polygonCoords <- list()
  
  polygonCoords[["min"]] <- c(99999999, 9999999)
  polygonCoords[["max"]] <- c(0, 0)
  
  for(i in 1:length(spatialPolygonsDataFrame@polygons)){
    polygonCoords[[spatialPolygonsDataFrame@polygons[[i]]@ID]] <- 
      spatialPolygonsDataFrame@polygons[[i]]@Polygons[[1]]@coords
    
    rangeX <- range(spatialPolygonsDataFrame@polygons[[i]]@Polygons[[1]]@coords[, 1])
    rangeY <- range(spatialPolygonsDataFrame@polygons[[i]]@Polygons[[1]]@coords[, 2])
    
    if(rangeX[1] < polygonCoords[["min"]][1]){
      polygonCoords[["min"]][1] <- rangeX[1]
    }
    if(rangeX[2] > polygonCoords[["max"]][1]){
      polygonCoords[["max"]][1] <- rangeX[2]
    }
    if(rangeY[1] < polygonCoords[["min"]][2]){
      polygonCoords[["min"]][2] <- rangeY[1]
    }
    if(rangeY[2] > polygonCoords[["max"]][2]){
      polygonCoords[["max"]][2] <- rangeY[2]
    }
  }
  
  return(polygonCoords)
}


