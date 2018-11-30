library(maptools) # Read shape file
#slibrary(rgeos) # Polygon centroids

###################################
# Read in ROI counties shape file # https://www.townlands.ie/page/download/
###################################

# Set the path
path <- "/home/josephcrispell/Desktop/Research/Paratuberculosis/"

# Read in the shape file
file <- paste(path, "CountyPolygonCoords/ROI_CountyBoundaries/counties.shp", sep="")
countyBorders <- readShapePoly(file) # Generates SpatialPolygonsDataFrame

# Get the coordinates of the counties
countyCoords <- getPolygonCoords(countyBorders)

# Get the county names associated with the coords
countyNames <- getPolygonNames(countyBorders@data, "NAME_TAG")

#########################
# Read in TP statistics #
#########################

# Read in the statistics file
file <- paste0(path, "HerdTbStatistics_2010-2018.csv")
statistics <- readTBStatisticsFile(file)

# Calculate the per county proportion test positive animals in most recent report - 2018Q3
summaryTables <- calculateSummaryStatisticsPerQuarter(statistics)

# # Read in the file
# file <- paste(path, "DAFM_TBHerdPrevalence_2015.txt", sep="")
# tbInfo <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE,
#                        check.names=FALSE)
# 
# # Combine counties that have been split
# tbInfo <- combineSplitCountyData(tbInfo)
# 
# # Get proportion herds infected per county
# countyProps <- getCountyPropHerdsInfected(tbInfo)

#####################
# Plot the counties #
#####################

plotProportionHerdsInfectedPerCounty(countyCoords, countyNames,
                                     countyProps)

#############
# FUNCTIONS #
#############

calculateSummaryStatisticsPerQuarter <- function(statistics){
  
  # Keep only the statistics of interest
  info <- statistics[grepl(statistics$Statistic, pattern="Tests on Animals|Total Yearly No. of Reactors to date|Herds Tested|Herds Restricted Since 1st of January"), ]
  
  # Examine each of the year quarters and store a summary
  output <- list()
  for(column in colnames(info)[-c(1,2)]){
    
    # Initialise a data frame to store a data summary
    dataSummary <- data.frame("County"=NA, "NumberAnimalsTested"=-1, "NumberAnimalsPositive"=-1,
                              "NumberHerdsTested"=-1, "NumberHerdsPositive"=-1, stringsAsFactors=FALSE)
    
    # Examine that TB statistics of interest
    row <- 0
    for(i in seq(from=4, to=nrow(info), by=4)){
      
      # Increment the row
      row <- row + 1
      
      # Store the information for the current county
      dataSummary[row, "County"] <- info[i, "County"]
      dataSummary[row, "NumberHerdsTested"] <- as.numeric(info[i-3, column])
      dataSummary[row, "NumberHerdsPositive"] <- as.numeric(info[i-2, column])
      dataSummary[row, "NumberAnimalsTested"] <- as.numeric(info[i-1, column])
      dataSummary[row, "NumberAnimalsPositive"] <- as.numeric(info[i, column])
    }
    
    # Calculate the proportion test positive
    dataSummary$ProportionAnimals <- dataSummary$NumberAnimalsPositive / dataSummary$NumberAnimalsTested
    dataSummary$ProportionHerds <- dataSummary$NumberHerdsPositive / dataSummary$NumberHerdsTested
    
    # Store the data summary for the current quarter
    output[[column]] <- dataSummary
  }
  
  return(output)
}

readTBStatisticsFile <- function(fileName){
  
  # Note that file was downloaded from: https://www.cso.ie/px/pxeirestat/Database/eirestat/Animal%20Disease%20Statistics/Animal%20Disease%20Statistics_statbank.asp?SP=Animal%20Disease%20Statistics&Planguage=0&ProductID=DB_DA
  # I selected all years and counties, downlaoded as csv and then removed quotes
  
  # Open a connection to a file to read (open="r")
  connection <- file(fileName, open="r")
  
  # Get all lines from file and store in vector
  fileLines <- readLines(connection)
  
  # Close file connection
  close(connection)
  
  # Intialise a dataframe to store the TB statistics
  statistics <- NULL
  county <- "NA"
  row <- 0
  
  # Loop through each of the lines in file
  for(i in 2:length(fileLines)){
    
    # Split the current line into its parts
    parts <- strsplit(fileLines[i], split=",")[[1]]
    
    # If 2nd line get the years initialise a dataframe to store the statistics
    if(i == 2){
      statistics <- data.frame(matrix(nrow=1, ncol=length(parts)), stringsAsFactors=FALSE)
      colnames(statistics) <- c("County", "Statistic", parts[-c(1,2)])
      next
    }
    
    # Check if found new county - name will present alone on new line
    if(length(parts) == 1){
      county <- parts[1]
      next
    }
    
    # Store the statistics from the current line
    row <- row + 1
    statistics[row, ] <- c(county, parts[-1])
  }
  
  # Convert the county names to upper case
  statistics$County <- toupper(statistics$County)
  
  return(statistics)
}

findMax <- function(list){
  
  max <- 0
  for(key in names(list)){
    if(list[[key]] > max){
      max <- list[[key]]
    }
  }
  
  return(max)
}

getCountyPropHerdsInfected <- function(tbInfo){
  
  countyProps <- list()
  for(row in 1:nrow(tbInfo)){
    
    countyProps[[tbInfo[row, 1]]] <- tbInfo[row, 4] / 100
  }
  
  return(countyProps)
}

combineSplitCountyData <- function(tbInfo){
  
  rowsToRemove <- c()
  for(row in 1:nrow(tbInfo)){
    
    parts <- strsplit(tbInfo[row, "RVO"], split=" ")[[1]]
    
    if(row %in% rowsToRemove){
      next
    }
    
    if(length(parts) > 1){
      
      newValues <- c()
      
      newValues[2] <- tbInfo[row, 2] + tbInfo[row + 1, 2]
      newValues[3] <- tbInfo[row, 3] + tbInfo[row + 1, 3]
      newValues[4] <- (((tbInfo[row, 3] * (tbInfo[row, 4] / 100)) +
                          (tbInfo[row + 1, 3] * (tbInfo[row + 1, 4] / 100))) /
                         newValues[3] ) * 100
      
      newValues[5] <- tbInfo[row, 5] + tbInfo[row + 1, 5]
      newValues[6] <- (((tbInfo[row, 5] * (tbInfo[row, 6] / 1000)) +
                          (tbInfo[row + 1, 5] * (tbInfo[row + 1, 6] / 1000))) /
                         newValues[5]) * 1000
      
      tbInfo[row, ] <- newValues
      tbInfo[row, 1] <- parts[1]
      rowsToRemove[length(rowsToRemove) + 1] <- row + 1
    }
  }
  tbInfo <- tbInfo[-rowsToRemove, ]
  
  return(tbInfo)
}

plotProportionHerdsInfectedPerCounty <- function(countyCoords, countyNames,
                                                 countyProps){
  
  niCounties <- c(
    "LONDONDERRY", "ANTRIM", "DOWN", "ARMAGH", "TYRONE", "FERMANAGH")
  
  # Calculate max proportion herds infected
  maxProp <- findMax(countyProps)
  
  par(mar=c(0,0,0,0))
  
  plot(x=NULL, y=NULL, yaxt="n", xaxt="n", ylab="", xlab="", bty="n",
       xlim=c(countyCoords[["min"]][1], countyCoords[["max"]][1]),
       ylim=c(countyCoords[["min"]][2], countyCoords[["max"]][2]))
  
  for(key in names(countyCoords)){
    
    if(key %in% c("min", "max")){
      next
    }
    
    if(countyNames[[key]] %in% niCounties){
      polygon(countyCoords[[key]], border=rgb(0,0,0, 1), 
              col=rgb(0,0,0, 0.75),
              lwd=2)
      next
    }
    
    prop <- countyProps[[countyNames[[key]]]]
    if(length(prop) == 0){
      print(countyNames[[key]])
      prop <- 0
    }
    
    polygon(countyCoords[[key]], border=rgb(0,0,0, 1), 
            col=rgb(1,0,0,prop / maxProp),
            lwd=2)
    
    if(prop != 0){
      text(x=mean(countyCoords[[key]][, 1]),
           y=mean(countyCoords[[key]][, 2]),
           labels=paste(round(prop * 100, digits=1), "%", sep=""),
           cex=0.75, col="blue")
    }else{
      print(countyNames[[key]])
    }
  }
}

writePolygonCoordsToFile <- function(countyCoords, countyNames, path){
  for(key in names(countyCoords)){
    
    if(key %in% c("min", "max")){
      next
    }
    
    file <- paste(path, "PolygonCoords_", countyNames[[key]], ".txt", sep="")
    table <- countyCoords[[key]]
    colnames(table) <- c("X", "Y")
    
    write.table(table, file, sep="\t",
                row.names=FALSE, quote=FALSE)
  }
}

getValues <- function(list){
  
  values <- c()
  for(key in names(list)){
    values[length(values) + 1] <- list[[key]]
  }
  
  return(values)
}

invertList <- function(list){
  
  output <- list()
  for(key in names(list)){
    output[[as.character(list[[key]])]] <- key
  }
  
  return(output)
}

countNumberSamplesPerCounty <- function(vntrInfo){
  counties <- list()
  for(row in 1:nrow(vntrInfo)){
    
    county <- strsplit(vntrInfo[row, "Herd Location"], split=" ")[[1]][1]
    if(is.null(counties[[county]]) == TRUE){
      counties[[county]] <- 1
    }else{
      counties[[county]] <- counties[[county]] + 1
    }
  }
  
  return(counties)
}

getPolygonNames <- function(polygonInfo, column){
  names <- list()
  rowNames <- rownames(polygonInfo)
  
  for(row in 1:nrow(polygonInfo)){
    names[[rowNames[row]]] <- toupper(as.character(polygonInfo[row, column]))
  }
  
  return(names)
}

removeSep <- function(string, sep){
  parts <- strsplit(x=string, split=sep)[[1]]
  return(paste(parts, collapse=""))
}

getPolygonCoords <- function(spatialPolygonsDataFrame){
  polygonCoords <- list()
  
  polygonCoords[["min"]] <- c(NA, NA)
  polygonCoords[["max"]] <- c(NA, NA)
  
  for(i in 1:length(spatialPolygonsDataFrame@polygons)){
    
    
    polygonCoords[[spatialPolygonsDataFrame@polygons[[i]]@ID]] <- 
      spatialPolygonsDataFrame@polygons[[i]]@Polygons[[1]]@coords
    
    rangeX <- range(spatialPolygonsDataFrame@polygons[[i]]@Polygons[[1]]@coords[, 1])
    rangeY <- range(spatialPolygonsDataFrame@polygons[[i]]@Polygons[[1]]@coords[, 2])
    
    if(is.na(polygonCoords[["min"]][1]) == TRUE ||
       rangeX[1] < polygonCoords[["min"]][1]){
      polygonCoords[["min"]][1] <- rangeX[1]
    }
    if(is.na(polygonCoords[["max"]][1]) == TRUE ||
       rangeX[2] > polygonCoords[["max"]][1]){
      polygonCoords[["max"]][1] <- rangeX[2]
    }
    if(is.na(polygonCoords[["min"]][2]) == TRUE ||
       rangeY[1] < polygonCoords[["min"]][2]){
      polygonCoords[["min"]][2] <- rangeY[1]
    }
    if(is.na(polygonCoords[["max"]][2]) == TRUE ||
       rangeY[2] > polygonCoords[["max"]][2]){
      polygonCoords[["max"]][2] <- rangeY[2]
    }
  }
  
  return(polygonCoords)
}
