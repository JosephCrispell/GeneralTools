##################
# Load Libraries #
##################

library(maptools) # Read shape file
library(rgeos) # Polygon centroids

# Create a path variable
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/"

############################################################
# Note the centroid of the badger social group territories #
############################################################

# Taken from PlotIsolateLocationsForBASTADemeAssignment.R
badgerGroupCentroid <- c(381761.7, 200964.3)
expand <- 1600

#####################################################################
# Load Infection Category Counts for each Social Group in each Year #
#####################################################################

# Read in the infection category counts for each social group
fileName <- paste(path, "BadgerCaptureData/", 
                  "InfectionCategoryCounts_2000-2011_10-08-2017.csv", sep="")
counts <- getCountTablesFromFileLinesYears(fileName)

# Calculate the proportion infected in each social group in each year
counts <- calculatePropInfectedInEachSocialGroupInYear(counts)

# Note the column that each social group is found in
socialGroupColumns <- noteColumnsOfSocialGroups(colnames(counts[["ProportionInfected"]]))

###########################################################
# Get Number Sequenced Isolates Per Social Group Per Year #
###########################################################

# Get list of isolate IDs from fasta file
fileName <- paste(path, "vcfFiles/", 
                  "sequences_withoutHomoplasies_27-03-18.fasta", sep="")
isolateIDs <- getSequenceIDsFromFastaFile(fileName)

# Read in the isolate metadata
fileName <- paste(path, "IsolateData/", "BadgerInfo_08-04-15_LatLongs_XY_Centroids.csv", sep="")
metadata <- read.table(fileName, header=TRUE, stringsAsFactors=FALSE, sep=",")

# Keep only rows for the isolates
metadata <- metadata[metadata$WB_id %in% isolateIDs, ]

# Count number sampled animals in each social group in each year
counts <- countNSamplesPerGroupPerYear(metadata, counts)

#########################################################
# Plot Infection Proportion on Social Group Territories #
#########################################################

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
groupsCentroidsPerYear <- list()

# Open a PDF file
file <- paste(path, "BadgerCaptureData/", 
              "PropInfected_SocialGroups_2000-11_29-03-18.pdf", sep="")

# Plot on single figure
plotPropInfectedInSocialGroupsInSingleFigure(path, file, years, shapeFileNames, badgerGroupCentroid,
                                             expand, counts, socialGroupColumns)

# Plot as a giff
plotPropInfectedInSocialGroupsInGiff(path, years, shapeFileNames, badgerGroupCentroid, expand, 
                                     counts, socialGroupColumns)


#############
# FUNCTIONS #
#############

plotPropInfectedInSocialGroupsInGiff <- function(path, years, shapeFileNames,
                                                 badgerGroupCentroid, expand, 
                                                 counts, socialGroupColumns){
  # Open a png set - distuinguish by plot number
  prefix <- paste(path, "BadgerCaptureData/PropInfected_SocialGroups_2000-11_Giff/SocialGroupSampling_%02d.png", sep="")
  png(file=prefix, height=480, width=600)
  
  # Examine the territories in each year
  for(i in 1:length(years)){
    
    year <- years[i]
    shapeFileName <- shapeFileNames[i]
    
    # Read in the shape file
    file <- paste(path, "BadgerCaptureData/BadgerTerritoryMarkingData/", 
                  "Baitmarking ", year, "/", shapeFileName, sep="")
    territories <- readShapePoly(file) # Generates SpatialPolygonsDataFrame
    
    # Extract the polygon coordinates
    territoryCoords <- getPolygonCoords(territories)
    
    # Get the full social group names
    territoryIDs <- getSocialGroupNames(territories@data, year)
    
    # Assign polygons to their social group names
    socialGroupTerritories <- assignTerritoriesToSocialGroupNames(territoryIDs, territoryCoords)
    
    # Calculate the territory centroids - mean X and Y
    territoryCentroids <- calculateTerritoryCentroids(territoryCoords, territoryIDs)
    
    # Plot the group territories - colour by proportion infected
    plotTerritories(socialGroupTerritories, badgerGroupCentroid, expand, year,
                    counts, socialGroupColumns, territoryCentroids)
    
    legend("bottomright", legend=c("Prop badgers infected",
                                   "N. badgers", "N. Sampled"), 
           text.col=c("grey25", "red", "blue"), bty="n")
    
    legend("bottom", legend=paste((expand * 2) / 1000, "km"), bty="n", cex=1.5)
  }
  
  # Close the PNG file output
  dev.off()
  
  # Bind the PNG files into a Giff
  dosPath <- "C:\\Users\\Joseph Crisp\\Desktop\\UbuntuSharedFolder\\Woodchester_CattleAndBadgers\\NewAnalyses_22-03-18\\BadgerCaptureData\\"
  system(paste("magick -delay 160 ", '\"', dosPath, "PropInfected_SocialGroups_2000-11_Giff\\SocialGroupSampling_*.png\" \"",
               dosPath, "PropInfected_SocialGroups_2000-11_Giff\\SocialGroupSampling.gif\"", sep=""))
  
  # Delete the PNG files
  unlink(paste(path, "BadgerCaptureData/PropInfected_SocialGroups_2000-11_Giff/SocialGroupSampling_*.png", sep=""))

}

plotPropInfectedInSocialGroupsInSingleFigure <- function(path, pdfFile, years, shapeFileNames,
                                             badgerGroupCentroid, expand, 
                                             counts, socialGroupColumns){
  
  pdf(pdfFile, height=11.69, width=8.27)
  
  par(mfrow=c(6,2))
  
  # Examine the territories in each year
  for(i in 1:length(years)){
    
    year <- years[i]
    shapeFileName <- shapeFileNames[i]
    
    # Read in the shape file
    file <- paste(path, "BadgerCaptureData/BadgerTerritoryMarkingData/", 
                  "Baitmarking ", year, "/", shapeFileName, sep="")
    territories <- readShapePoly(file) # Generates SpatialPolygonsDataFrame
    
    # Extract the polygon coordinates
    territoryCoords <- getPolygonCoords(territories)
    
    # Get the full social group names
    territoryIDs <- getSocialGroupNames(territories@data, year)
    
    # Assign polygons to their social group names
    socialGroupTerritories <- assignTerritoriesToSocialGroupNames(territoryIDs, territoryCoords)
    
    # Calculate the territory centroids - mean X and Y
    territoryCentroids <- calculateTerritoryCentroids(territoryCoords, territoryIDs)
    
    # Plot the group territories - colour by proportion infected
    plotTerritories(socialGroupTerritories, badgerGroupCentroid, expand, year,
                    counts, socialGroupColumns, territoryCentroids)
  }
  
  dev.off()
}

countNSamplesPerGroupPerYear <- function(metadata, counts){
  sampled <- counts[["Negative"]]
  sampled[1:nrow(sampled), 2:ncol(sampled)] <- 0
  
  badgersInYearsAlreadyCounted <- list()
  for(year in sampled[1, "Years"]:sampled[nrow(sampled), "Years"]){
    badgersInYearsAlreadyCounted[[as.character(year)]] <- list()
  }
  
  for(row in 1:nrow(metadata)){
    
    year <- as.numeric(strsplit(metadata[row, "date"], split="/")[[1]][3])
    group <- paste(strsplit(metadata[row, "Social.Group.Trapped.At"], split=" ")[[1]], collapse="")

    if(group %in% colnames(sampled) == TRUE){
      
      if(is.null(badgersInYearsAlreadyCounted[[as.character(year)]][[metadata[row, "tattoo"]]]) == TRUE){
        sampled[(year - sampled[1, "Years"]) + 1, group] <- as.numeric(sampled[(year - sampled[1, "Years"]) + 1, group]) + 1
        badgersInYearsAlreadyCounted[[as.character(year)]][[metadata[row, "tattoo"]]] <- 1
      }
    }
  }
  
  counts[["Sampled"]] <- sampled
  
  return(counts)
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

divide <- function(dataframeA, dataframeB){
  
  result <- dataframeA
  
  for(row in 1:nrow(dataframeA)){
    for(col in 1:ncol(dataframeA)){
      result[row, col] <- as.numeric(dataframeA[row, col]) / as.numeric(dataframeB[row, col])
    }
  }
  
  return(result)
}

add <- function(dataframeA, dataframeB){
  
  result <- dataframeA
  
  for(row in 1:nrow(dataframeA)){
    for(col in 1:ncol(dataframeA)){
      result[row, col] <- as.numeric(dataframeA[row, col]) + as.numeric(dataframeB[row, col])
    }
  }
  
  return(result)
}

getCountTablesFromFileLinesYears <- function(fileName){
  
  connection <- file(fileName, open='r')
  fileLines <- readLines(connection)
  close(connection)
  
  counts <- list()
  counts[["Negative"]] <- data.frame(Years=rep(0, length(fileLines) - 1), stringsAsFactors=FALSE)
  counts[["Exposed"]] <- data.frame(Years=rep(0, length(fileLines) - 1), stringsAsFactors=FALSE)
  counts[["Excretor"]] <- data.frame(Years=rep(0, length(fileLines) - 1), stringsAsFactors=FALSE)
  counts[["Superexcretor"]] <- data.frame(Years=rep(0, length(fileLines) - 1), stringsAsFactors=FALSE)
  
  # Get the social group names from the header
  socialGroups <- strsplit(fileLines[1], "\t")[[1]][-1]
  
  # Read in the tables
  for(row in 2:length(fileLines)){
    cols <- strsplit(fileLines[row], "\t")[[1]]
    
    counts[["Negative"]][row-1, "Years"] <- as.numeric(cols[1])
    counts[["Exposed"]][row-1, "Years"] <- as.numeric(cols[1])
    counts[["Excretor"]][row-1, "Years"] <- as.numeric(cols[1])
    counts[["Superexcretor"]][row-1, "Years"] <- as.numeric(cols[1])
    
    for(col in 2:length(cols)){
      
      parts <- strsplit(cols[col], ":")[[1]]
      counts[["Negative"]][row-1, col] <- parts[1]
      counts[["Exposed"]][row-1, col] <- parts[2]
      counts[["Excretor"]][row-1, col] <- parts[3]
      counts[["Superexcretor"]][row-1, col] <- parts[4]
      
    }
  }
  
  # Add the column names
  colnames(counts[["Negative"]]) <- c("Years", socialGroups)
  colnames(counts[["Exposed"]]) <- c("Years", socialGroups)
  colnames(counts[["Excretor"]]) <- c("Years", socialGroups)
  colnames(counts[["Superexcretor"]]) <- c("Years", socialGroups)
  
  return(counts)
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

noteColumnsOfSocialGroups <- function(columns){
  
  socialGroupCols <- list()
  for(i in 1:length(columns)){
    if(columns[i] == "Years"){
      next
    }
    
    socialGroupCols[[columns[i]]] <- i
  }
  
  return(socialGroupCols)
}

calculatePropInfectedInEachSocialGroupInYear <- function(counts){
  
  # Subset each table to get counts for the current year
  negative <- counts[["Negative"]]
  exposed <- counts[["Exposed"]]
  excretor <- counts[["Excretor"]]
  superExcretor <- counts[["Superexcretor"]]
  
  # Calculate the proportion infected in each social group
  totalInfected <- add(add(exposed[, -1], excretor[, -1]), superExcretor[, -1])
  total <- add(totalInfected, negative[, -1])
  propInfected <- divide(totalInfected, total)
  
  # Add the dates column onto table
  propInfected$Years <- negative$Years
  total$Years <- negative$Years
  
  counts[["ProportionInfected"]] <- propInfected
  counts[["Total"]] <- total
  
  return(counts)
}

plotTerritories <- function(socialGroupTerritories, badgerCentre, expand, year,
                            counts, socialGroupColumns, territoryCentroids){
  
  par(mar=c(0,0,0,0))
  plot(x=NULL, y=NULL, yaxt="n", xaxt="n", bty="n", ylab="",
       xlim=c(badgerCentre[1] - expand, badgerCentre[1] + expand), 
       ylim=c(badgerCentre[2] - expand, badgerCentre[2] + expand), 
       asp=1, xlab="")
  legend("topright", legend=year, bty="n", cex=2)

  propInfected <- counts[["ProportionInfected"]]
  total <- counts[["Total"]]
  sampled <- counts[["Sampled"]]
  
  for(groupName in names(socialGroupTerritories)){
    
    alpha <- as.numeric(propInfected[(year - propInfected[1, "Years"]) + 1, groupName])
    n <- total[(year - propInfected[1, "Years"]) + 1, groupName]
    
    if(length(alpha) == 0 || is.nan(alpha) == TRUE){
      alpha <- 0
      n <- ""
    }
    
    nSampled <- as.numeric(sampled[(year - sampled[1, "Years"]) + 1, groupName])
    
    polygon(socialGroupTerritories[[groupName]], col=rgb(0,0,0, alpha=alpha))
    
    if(length(nSampled) != 0 && nSampled != 0){
      
      sampledLabel <- paste("(", nSampled, ")", sep="")
      
      text(x=territoryCentroids[[groupName]][1], y=territoryCentroids[[groupName]][2], 
           labels=bquote(bold(.(n)) * phantom(.(sampledLabel))), col="red", cex=1)
      text(x=territoryCentroids[[groupName]][1], y=territoryCentroids[[groupName]][2], 
           labels=bquote(bold(phantom(.(n)) * .(sampledLabel))), col="blue", cex=1)
      
    }else{
      text(x=territoryCentroids[[groupName]][1], y=territoryCentroids[[groupName]][2], 
           labels=n, font=2, col="red", cex=1)
    }
  }
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

