library(maptools) # Read shape file
#slibrary(rgeos) # Polygon centroids

###################################
# Read in ROI counties shape file # https://www.townlands.ie/page/download/
###################################

# Set the path
path <- "/home/josephcrispell/storage/Research/"

# Read in the shape file
file <- paste(path, "ROI_CountyBoundaries/counties.shp", sep="")
countyBorders <- readShapePoly(file) # Generates SpatialPolygonsDataFrame

# Get the coordinates of the counties
countyCoords <- getPolygonCoords(countyBorders)

# Get the county names associated with the coords
countyNames <- getPolygonNames(countyBorders@data, "NAME_TAG")

#########################
# Read in TP statistics #
#########################

# Read in the statistics file (INCORRECT VALUE FOR WICKLOW WEST IN 2015Q2 1590/310 herds restricted
file <- paste0(path, "RepublicOfIreland/Mbovis/HerdTbStatistics_2010-2019.csv") # source: https://statbank.cso.ie/px/pxeirestat/Database/eirestat/Animal%20Disease%20Statistics/Animal%20Disease%20Statistics_statbank.asp?sp=Animal%20Disease%20Statistics&Planguage=0&ProductID=DB_DA
statistics <- readTBStatisticsFile(file)
statistics[which(statistics$County == "WICKLOW W" & statistics$Statistic == "Herds Restricted Since 1st of January (Number)"),
           "2015Q2"] <- NA

# Calculate the per county proportion test positive animals in most recent report - 2018Q3
summaryTables <- calculateSummaryStatisticsPerQuarter(statistics)

##################
# Plot summaries #
#################

# Open an output pdf
pdf(paste0(path, "RepublicOfIreland/Mbovis/HerdTbStatistics_2010-2019.pdf"))

# Plot the proportion of herds infected in the country - compare to number of herds
plotInfoPerCounty(countyCoords, countyNames, summaryTables$`2019Q3`, column="ProportionHerds", digits=1, units="%",
                  main="Proportion herds infected", cex=0.75)
plotInfoPerCounty(countyCoords, countyNames, summaryTables$`2019Q3`, column="NumberHerdsTested", units="",
                  main="Number of herds", cex=0.75)
plotInfoPerCounty(countyCoords, countyNames, summaryTables$`2019Q3`, column="NumberHerdsPositive", units="",
                  main="Number of test positive herds", cex=0.75)

# Plot the trend in the proportion of infected
yearsOfInterest <- 2010:2019
plotProportionTrends(summaryTables, yearsOfInterest)

# Plot the trends in each county
par(mar=c(5.1, 6.1, 1.1, 2.1))
plotCountyTrends(summaryTables, column="ProportionHerds")
plotCountyTrends(summaryTables, column="ProportionHerds", county="WICKLOW", county.lwd=2, county.col="red")
plotCountyTrends(summaryTables, column="ProportionHerds", county="MONAGHAN", county.lwd=2, county.col="red")

dev.off()

#############
# FUNCTIONS #
#############

plotCountyTrends <- function(summaryTables, column, county=NULL, county.lwd=1, county.col=rgb(1,0,0,1)){

  # Initialise a table to store the column values per county per quarter
  values <- data.frame(matrix(0, nrow=nrow(summaryTables[[1]]), ncol=length(summaryTables)))
  rownames(values) <- summaryTables[[1]]$County
  colnames(values) <- names(summaryTables)

  # Initialise a vector to store the quarters stored as numeric values
  quarters <- c(3, 6, 9, 12)/12
  times <- c()

  # Examine data from each quarter
  for(quarter in names(summaryTables)){

    # Store the current quarter as a time in year
    parts <- as.numeric(strsplit(quarter, split="Q")[[1]])
    times[length(times) + 1] <- parts[1] + quarters[parts[2]]

    # Store the values of the selected column
    values[, quarter] <- summaryTables[[quarter]][, column]
  }

  # Plot the trends for the selected column
  plot(x=times, y=values["STATE", ], ylim=range(values, na.rm=TRUE), las=1, bty="n", type="l", lwd=4,
       ylab="", xlab="Year", col="blue",
       cex.axis=1.5, cex.lab=2, cex.main=2)
  mtext(side=2, text="Proportion herds testing positive", cex=2, line=4)

  # Plot the trends for each county
  for(row in 2:nrow(values)){

    points(x=times, y=values[row, ], type="l",
           col=ifelse(is.null(county) == FALSE && rownames(values)[row] == county, county.col, rgb(0,0,0,0.5)),
           lwd=ifelse(is.null(county) == FALSE && rownames(values)[row] == county, county.lwd, 1))
  }

  # Add a legend
  if(is.null(county) == FALSE){
    legend("top", legend=c("Individual county trends", "Country average", county),
           text.col=c(rgb(0,0,0,0.5), "blue", county.col),
           bty="n")
  }else{
    legend("top", legend=c("Individual county trends", "Country average"), text.col=c(rgb(0,0,0,0.5), "blue"), bty="n")
  }

}

plotCounty <- function(countyCoords, countyNames, name){

  # Find the county index
  countyIndex <- NULL
  for(key in names(countyNames)){

    if(countyNames[[key]] == name){
      countyIndex <- key
      break
    }
  }

  # Create an empty plot
  plot(x=NULL, y=NULL, yaxt="n", xaxt="n", ylab="", xlab="", bty="n",
       xlim=range(countyCoords[[countyIndex]][, 1]), ylim=range(countyCoords[[countyIndex]][, 2]))

  # Plot the county polygon
  polygon(countyCoords[[countyIndex]], border=rgb(0,0,0, 1), col=rgb(0,0,1, 0), lwd=4)
}

plotCountyOutlines <- function(countyCoords, countyNames){

  # Get and set the margins
  currentMar <- par("mar")
  par(mar=c(0,0,0,0))

  # Note the names of the counties in Northern Ireland
  niCounties <- c("LONDONDERRY", "ANTRIM", "DOWN", "ARMAGH", "TYRONE", "FERMANAGH")

  # Create an empty plot
  plot(x=NULL, y=NULL, yaxt="n", xaxt="n", ylab="", xlab="", bty="n",
       xlim=c(countyCoords[["min"]][1], countyCoords[["max"]][1]),
       ylim=c(countyCoords[["min"]][2], countyCoords[["max"]][2]))

  # Examine each of the counties
  for(key in names(countyCoords)){

    # Ignore the minimum and maximum coordinates
    if(key %in% c("min", "max")){
      next
    }

    # Plot the Northern Ireland counties greyed out
    if(countyNames[[key]] %in% niCounties){
      polygon(countyCoords[[key]], border=rgb(0,0,0, 1), col=rgb(0,0,0, 0.75), lwd=2)
      next
    }

    # Plot a polygon for the current county
    polygon(countyCoords[[key]], border=rgb(0,0,0, 1), col=rgb(0,0,1, 0), lwd=2)
  }

  # Reset the margins
  par(mar=currentMar)
}

plotProportionTrends <- function(summaryTables, yearsOfInterest){

  # Initialise a vector to store the proportion of infected herds and the number of herds
  proportionHerds <- c()
  nHerds <- c()
  proportionAnimals <- c()
  nAnimals <- c()

  # Examine the summary tables available for each period
  for(key in names(summaryTables)){

    # Skip all quarters except the fourth (or third for 2018)
    if(grepl(key, pattern="Q1|Q2|Q3") && key != "2018Q3"){
      next
    }

    # Get the information for the state
    stateInfo <- summaryTables[[key]][1, ]

    # Store the proportion of herds infected and the number of herds
    proportionHerds[length(proportionHerds) + 1] <- stateInfo[1, "ProportionHerds"]
    nHerds[length(nHerds) + 1] <- stateInfo[1, "NumberHerdsTested"]

    # Store the proportion of animals infected and the number of animals
    proportionAnimals[length(proportionAnimals) + 1] <- stateInfo[1, "ProportionAnimals"]
    nAnimals[length(nAnimals) + 1] <- stateInfo[1, "NumberAnimalsTested"]
  }

  # Plot the proportion of herds infected in Ireland
  plot(x=yearsOfInterest, y=proportionHerds, ylim=c(0, 0.05), las=1, bty="n", type="o", lwd=2, pch=19, cex=2,
       ylab="Proportion", xlab="Year", main="Proportion tested herds positive")

  # Get the axis limits
  axisLimits <- par("usr")

  # Add the number of herds above each proportion
  text(x=yearsOfInterest, y=proportionHerds + 0.1*proportionHerds, labels=nHerds, xpd=TRUE)

  # Add line for European threshold
  points(x=c(yearsOfInterest[1], yearsOfInterest[length(yearsOfInterest)]), y=c(0.001, 0.001), type="l", lty=2,
         col="green", lwd=4)


  # Get and set the margins
  currentMar <- par("mar")
  par(mar=c(5.1,5.1,4.1,2.1))

  # Plot the proportion of animals infected in Ireland
  plot(x=yearsOfInterest, y=proportionAnimals, ylim=c(0, 0.003), las=1, bty="n", type="o", lwd=2, pch=19, cex=2,
       ylab="", xlab="Year", main="Proportion tested animals positive")

  # Add Y axis label
  mtext("Proportion", side=2, line=4)

  # Get the axis limits
  axisLimits <- par("usr")

  # Add the number of herds above each proportion
  text(x=yearsOfInterest, y=proportionAnimals + 0.1*proportionAnimals, labels=nAnimals, xpd=TRUE, cex=0.8)

  # Reset the margins
  par(mar=currentMar)
}

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

    # Combine Tipperary (north & south), Cork (north & south), and Wicklow (east and west)
    corkRow <- which(dataSummary$County == "CORK NORTH")
    dataSummary[corkRow, "County"] <- "CORK"
    dataSummary[corkRow, 2:5] <- dataSummary[corkRow, 2:5] + dataSummary[corkRow + 1, 2:5]
    dataSummary <- dataSummary[-(corkRow + 1), ]

    tipperaryRow <- which(dataSummary$County == "TIPPERARY NORTH")
    dataSummary[tipperaryRow, "County"] <- "TIPPERARY"
    dataSummary[tipperaryRow, 2:5] <- dataSummary[tipperaryRow, 2:5] + dataSummary[tipperaryRow + 1, 2:5]
    dataSummary <- dataSummary[-(tipperaryRow + 1), ]

    wicklowRow <- which(dataSummary$County == "WICKLOW E")
    dataSummary[wicklowRow, "County"] <- "WICKLOW"
    dataSummary[wicklowRow, 2:5] <- dataSummary[wicklowRow, 2:5] + dataSummary[wicklowRow + 1, 2:5]
    dataSummary <- dataSummary[-(wicklowRow + 1), ]

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
  for(i in 4:length(fileLines)){

    # Remove quotes
    fileLines[i] <- gsub(pattern="\"", replacement="", x=fileLines[i])

    # Split the current line into its parts
    parts <- strsplit(fileLines[i], split=",")[[1]]

    # If 4th line get the years initialise a dataframe to store the statistics
    if(i == 4){
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

plotInfoPerCounty <- function(countyCoords, countyNames, summaryTable, column, digits=1, units="%",
                              cex=1, main=""){

  # Get and set the margins
  currentMar <- par("mar")
  par(mar=c(0,0,2,0))

  # Note the names of the counties in Northern Ireland
  niCounties <- c("LONDONDERRY", "ANTRIM", "DOWN", "ARMAGH", "TYRONE", "FERMANAGH")

  # Calculate max value in column
  max <- max(summaryTable[-1, column])

  # Create an empty plot
  plot(x=NULL, y=NULL, yaxt="n", xaxt="n", ylab="", xlab="", bty="n",
       xlim=c(countyCoords[["min"]][1], countyCoords[["max"]][1]),
       ylim=c(countyCoords[["min"]][2], countyCoords[["max"]][2]),
       main=main)

  # Examine each of the counties
  for(key in names(countyCoords)){

    # Ignore the minimum and maximum coordinates
    if(key %in% c("min", "max")){
      next
    }

    # Plot the Northern Ireland counties greyed out
    if(countyNames[[key]] %in% niCounties){
      polygon(countyCoords[[key]], border=rgb(0,0,0, 1), col=rgb(0,0,0, 0.75), lwd=2)
      next
    }

    # Get the information for the current county
    countyInfo <- summaryTable[summaryTable$County == countyNames[[key]], ]

    # Plot a polygon for the current county
    polygon(countyCoords[[key]], border=rgb(0,0,0, 1),
            col=rgb(0,0,1, countyInfo[1, column] / max), lwd=2)

    # Find the mid-point of the polygon
    middleX <- mean(countyCoords[[key]][, 1])
    middleY <- mean(countyCoords[[key]][, 2])
  }

  # Add a legend
  range <- range(summaryTable[-1, column])
  legendValues <- seq(from=range[1], to=range[2], length.out=5)
  legend("bottomright", legend=signif(legendValues, digits=2), pch=22,
         pt.bg=c(rgb(0,0,1, legendValues[1]/max),
                 rgb(0,0,1, legendValues[2]/max),
                 rgb(0,0,1, legendValues[3]/max),
                 rgb(0,0,1, legendValues[4]/max),
                 rgb(0,0,1, legendValues[5]/max)),
         bty="n", pt.cex=2, col="black")

  # Reset the margins
  par(mar=currentMar)
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
