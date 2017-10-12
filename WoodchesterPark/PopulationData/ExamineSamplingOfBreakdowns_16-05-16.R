path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_13-07-17/"

##############################################
# Read in Breakdown and Location Information #
##############################################

######################
# Location Information
file <- paste(path, "CattleMovementData-Post2001/20160314_joe_cts_locations.csv", sep="")
locationInfo <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep=",", fill=TRUE)

# Remove CPHs that aren't within 15km of Woodchester Mansion
badgerCentre <- c(381761.7, 200964.3)
threshold <- 15000

locationInfo <- keepLocationsWithinThresholdDistance(locationInfo=locationInfo, 
                                                     thresholdInMetres=threshold,
                                                     badgerCentre=badgerCentre)

# Get the CPH into the right format
locationInfo$CPH <- formatCPHForLocationInfo(locationInfo$cph)

# Convert the dataframe into a list by CPH
locationList <- convertTableToList(column="CPH", table=locationInfo)


#######################
# Breakdown Information
file <- paste(path, "CattleTestData/tblccdBreakdown_22-11-16.csv", sep="")
breakdownInfo <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep=",", fill=TRUE)

# Keep information for locations within 15km of Woodchester Park
breakdownInfo <- keepRowsIfColumnValuePresentInList(column="cbCph", table=breakdownInfo,
                                                    listStructure=locationList)

# Add X, Y columns to breakdownInfo
breakdownInfo[, c("X", "Y", "DistanceToWoodchesterPark")] <- rep(0, nrow(breakdownInfo))

for(row in 1:nrow(breakdownInfo)){
  
  infoForLocation <- locationList[[as.character(breakdownInfo[row, "cbCph"])]]
  breakdownInfo[row, "X"] <- infoForLocation[1, "x"]
  breakdownInfo[row, "Y"] <- infoForLocation[1, "y"]
  breakdownInfo[row, "DistanceToWoodchesterPark"] <- infoForLocation[1, "DistanceToWoodchester"]
}

####################################
# Read in the Sampling Information #
####################################

file <- paste(path, "IsolateData/",
              "CattleIsolateInfo_LatLongs_plusID_outbreakSize_Coverage_AddedStrainIDs.csv",
              sep="")
samplingInfo <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep=",")

# Remove rows with no StrainID
samplingInfo <- samplingInfo[is.na(samplingInfo$StrainId) == FALSE, ]

##############################
# Examine Breakdown Sampling #
##############################

# Add empty column to the breakdown info table
breakdownInfo$NSequencedIsolates <- rep(0, nrow(breakdownInfo))

# Convert the dataframe into a list by breakdown ID
breakdownList <- convertTableToList(column="cbBreakId", table=breakdownInfo)

# Note the number of isolates sequenced from each breakdown
for(row in 1:nrow(samplingInfo)){
  
  # Find the breakdown information
  if(is.null(breakdownList[[samplingInfo[row, "BreakdownID"]]]) == FALSE){
    
    # Add to count for sequenced isolates from this breakdown
    info <- breakdownList[[samplingInfo[row, "BreakdownID"]]]
    info[1, "NSequencedIsolates"] <- info[1, "NSequencedIsolates"] + 1
    breakdownList[[samplingInfo[row, "BreakdownID"]]] <- info
    
  }else{
    
    # 14101000701-24/02/2006 (SL Case) Couldn't be found
    # Changed to later Breakdown: 14101000701-15/12/2006
    print(paste("Ahhhhhhhhhhhhhhhh! Can't find ", samplingInfo[row, "BreakdownID"], sep=""))
  }
}

# Convert the breakdown list into data.frame
breakdownsTable <- breakdownList[[1]]
keys <- names(breakdownList)
for(index in 1:length(keys)){
  
  breakdownsTable[index, ] <- breakdownList[[keys[index]]]
}

# Remove breakdowns that occur outside of the sampling window
dateRange <- getSampledBreakdownDateRange(samplingInfo)
breakdownsTable <- removeBreakdownsThatOccurOutsideOfSamplingWindow(breakdownsTable, dateRange)
breakdownsTable$Date <- as.Date(breakdownsTable$Date)

# Remove non-confirmed breakdowns
breakdownsTable <- breakdownsTable[breakdownsTable$cbConfFin == "Y", ]

#######################
# Produce a Figure... #
#######################

# Spatial plot - shape by herd type and size by the number of samples
expand <- 15000

file <- paste(path, "BreakdownSampling_02-10-17.pdf")
pdf(file, height=14, width=21)

par(mfrow=c(2,3))

par(mar=c(6.1,6.1,7.1,2.1)) # bottom, left, top and right

holdingColour <- rgb(1,0,0, 0.25)

### BREAKDOWNS

# Create empty plot
plot(1, type="n", yaxt="n", xaxt="n",
     xlim=c(badgerCentre[1] - expand, badgerCentre[1] + expand),
     ylim=c(badgerCentre[2] - expand, badgerCentre[2] + expand),
     xlab=paste((expand[1] * 2) / 1000, "km"), ylab=paste((expand[1] * 2) / 1000, "km"),
     main="Locations of Cattle Reactors", cex.lab=2, cex.main=3)

# Add point for Woodchester Mansion
points(x=badgerCentre[1], y=badgerCentre[2], pch=15, col="blue", cex=3)

# Add points for each cattle herd
points(x=locationInfo$x, y=locationInfo$y, 
       col=rgb(0,0,0, 0.2),
       pch=19, cex=1)
points(x=breakdownsTable$X, y=breakdownsTable$Y, 
       col=holdingColour,
       pch=17, cex=log(breakdownsTable$cbNumCattleReactors)+1)

# Add Legend
legend(x=365909, y=216377, legend=c(200, 100, 10, 1), pch=c(17, 17, 17, 17), bty="n", 
       pt.cex=log(c(200, 100, 10, 1))+1, col=rgb(0,0,0, 0.9), 
       y.intersp=2.5, x.intersp=2, cex=2)
legend("bottomleft", legend=c("COUNT", "MANSION", "HOLDING"), pch=c(17,15, 19), 
       col=c("red", "blue", "dimgrey"), text.col=c("red", "blue", "dimgrey"),
       bty="n", cex=1.5)

### CULTURE

holdingColour <- rgb(1,0,0, 0.3)

# Create empty plot
plot(1, type="n", yaxt="n", xaxt="n",
     xlim=c(badgerCentre[1] - expand, badgerCentre[1] + expand),
     ylim=c(badgerCentre[2] - expand, badgerCentre[2] + expand),
     xlab=paste((expand[1] * 2) / 1000, "km"), ylab=paste((expand[1] * 2) / 1000, "km"),
     main="Locations of Cultured Reactors", cex.lab=2, cex.main=3)

# Add point for Woodchester Mansion
points(x=badgerCentre[1], y=badgerCentre[2], pch=15, col="blue", cex=3)

# Add points for each cattle herd
points(x=locationInfo$x, y=locationInfo$y, 
       col=rgb(0,0,0, 0.2),
       pch=19, cex=1)
points(x=breakdownsTable$X, y=breakdownsTable$Y, 
       col=holdingColour,
       pch=17,
       cex=(breakdownsTable$cbNumCattleCultPos/4)+1)

# Add Legend
legend(x=365909, y=216377, legend=c(25, 10, 5, 1), pch=c(17, 17, 17, 17), bty="n", 
       pt.cex=(c(25, 10, 5, 1)/4)+1, col=rgb(0,0,0, 0.9),
       y.intersp=2.5, x.intersp=2, cex=2)
legend("bottomleft", legend=c("COUNT", "MANSION", "HOLDING"), pch=c(17, 15, 19), 
       col=c("red", "blue", "dimgrey"), text.col=c("red", "blue", "dimgrey"),
       bty="n", cex=1.5)

### SAMPLING

holdingColour <- rgb(1,0,0, 0.5)

# Create empty plot
plot(1, type="n", yaxt="n", xaxt="n",
     xlim=c(badgerCentre[1] - expand, badgerCentre[1] + expand),
     ylim=c(badgerCentre[2] - expand, badgerCentre[2] + expand),
     xlab=paste((expand[1] * 2) / 1000, "km"), ylab=paste((expand[1] * 2) / 1000, "km"),
     main="Locations of Sequenced Isolates", cex.lab=2, cex.main=3)

# Add point for Woodchester Mansion
points(x=badgerCentre[1], y=badgerCentre[2], pch=15, col="blue", cex=3)

# Add points for each cattle herd
points(x=locationInfo$x, y=locationInfo$y, 
       col=rgb(0,0,0, 0.2),
       pch=19, cex=1)
points(x=breakdownsTable$X, y=breakdownsTable$Y, 
       col=holdingColour,
       pch=17,
       cex=breakdownsTable$NSequencedIsolates*1.5)

# Add Legend
legend(x=365909, y=216377, legend=c(4, 3, 2, 1), pch=c(17, 17, 17, 17), bty="n", 
       pt.cex=c(4, 3, 2, 1)*1.5, col=rgb(0,0,0, 0.9),
       y.intersp=2.5, x.intersp=2, cex=2)
legend("bottomleft", legend=c("COUNT", "MANSION", "HOLDING"), pch=c(17, 15, 19), 
       col=c("red", "blue", "dimgrey"), text.col=c("red", "blue", "dimgrey"),
       bty="n", cex=1.5)


### By Date

countsForTimestep <- getReactorCultureAndSequencedCounts(dateRange, "3 months", breakdownsTable)

polygonColour <- rgb(0,0,0, 0)
pointsCol <- rgb(0,0,0, 0.25)
pointsSize <- 1.5
lineWidth <- 3
lineCol <- rgb(0,0,1, 1)

### Reactors
plotCountsThroughTime(column="cbNumCattleReactors", columnInCounts="NReactors",
                      countsForTimestep=countsForTimestep, log=TRUE, 
                      title="Number of Skin Test Reactors", breakdownsTable=breakdownsTable,
                      polygonColour=polygonColour, pointsCol=pointsCol,
                      pointsSize=pointsSize, lineWidth=lineWidth, lineCol=lineCol)

### Culture
plotCountsThroughTime(column="cbNumCattleCultPos", columnInCounts="NCultured",
                      countsForTimestep=countsForTimestep, log=FALSE, 
                      title="Number of Culture Positive Reactors", breakdownsTable=breakdownsTable,
                      polygonColour=polygonColour, pointsCol=pointsCol,
                      pointsSize=pointsSize, lineWidth=lineWidth, lineCol=lineCol)

### Sequenced Isolates
plotCountsThroughTime(column="NSequencedIsolates", columnInCounts="NSequenced",
                      countsForTimestep=countsForTimestep, log=FALSE, 
                      title="Number of Sequenced Isolates", breakdownsTable=breakdownsTable,
                      polygonColour=polygonColour, pointsCol=pointsCol,
                      pointsSize=pointsSize, lineWidth=lineWidth, lineCol=lineCol)

dev.off()


######################################
# Dynamic figure with space and time #
######################################

### BREAKDOWNS

prefix <- paste(path, "CattleTestData/DynamicGiff/BreakdownSampling_%02d.png", sep="")
png(file=prefix, height=480, width=960)
par(mfrow=c(1,2))
par(mar=c(6.1,6.1,7.1,2.1)) # bottom, left, top and right

breakdownsTable$Year <- format(breakdownsTable$Date, "%Y")

for(year in 1997:2012){
  
  # Get the data subset for current year
  subset <- breakdownsTable[breakdownsTable$Year == year, ]
  
  ### Breakdowns
  
  holdingColour <- rgb(1,0,0, 0.5)
  
  # Create empty plot
  plot(1, type="n", yaxt="n", xaxt="n", bty="n",
       xlim=c(badgerCentre[1] - expand, badgerCentre[1] + expand),
       ylim=c(badgerCentre[2] - expand, badgerCentre[2] + expand),
       xlab=paste((expand * 2) / 1000, "km"), ylab=paste((expand * 2) / 1000, "km"),
       main="Locations of Cattle Reactors", cex.lab=2, cex.main=2)
  
  # Add year
  text(x=badgerCentre[1] + expand, y=badgerCentre[2] + expand, labels=year, xpd=TRUE, cex=2)
  
  # Add point for Woodchester Mansion
  points(x=badgerCentre[1], y=badgerCentre[2], pch=15, col="blue", cex=3)
  
  # Add points for each cattle herd
  points(x=locationInfo$x, y=locationInfo$y, 
         col=rgb(0,0,0, 0.2),
         pch=19, cex=1)
  points(x=subset$X, y=subset$Y, 
         col=holdingColour,
         pch=17, cex=log(subset$cbNumCattleReactors)+1)
  
  # Add Legend
  legend(x=365909, y=216377, legend=c(200, 100, 10, 1), pch=c(17, 17, 17, 17), bty="n", 
         pt.cex=log(c(200, 100, 10, 1))+1, col=rgb(0,0,0, 0.9), 
         y.intersp=2.5, x.intersp=2, cex=1.05, xpd=TRUE)
  legend("bottomleft", legend=c("COUNT", "MANSION", "HOLDING"), pch=c(17,15, 19), 
         col=c("red", "blue", "dimgrey"), text.col=c("red", "blue", "dimgrey"),
         bty="n", cex=0.75)
  
  ### SAMPLING
  
  holdingColour <- rgb(1,0,0, 0.75)
  
  # Create empty plot
  plot(1, type="n", yaxt="n", xaxt="n", bty="n",
       xlim=c(badgerCentre[1] - expand, badgerCentre[1] + expand),
       ylim=c(badgerCentre[2] - expand, badgerCentre[2] + expand),
       xlab=paste((expand * 2) / 1000, "km"), ylab=paste((expand * 2) / 1000, "km"),
       main="Locations of Sequenced Isolates", cex.lab=2, cex.main=2)
  
  # Add year
  text(x=badgerCentre[1] + expand, y=badgerCentre[2] + expand, labels=year, xpd=TRUE, cex=2)
  
  # Add point for Woodchester Mansion
  points(x=badgerCentre[1], y=badgerCentre[2], pch=15, col="blue", cex=3)
  
  # Add points for each cattle herd
  points(x=locationInfo$x, y=locationInfo$y, 
         col=rgb(0,0,0, 0.2),
         pch=19, cex=1)
  points(x=subset$X, y=subset$Y, 
         col=holdingColour,
         pch=17,
         cex=subset$NSequencedIsolates*1.5)
  
  # Add Legend
  legend(x=365909, y=216377, legend=c(4, 3, 2, 1), pch=c(17, 17, 17, 17), bty="n", 
         pt.cex=c(4, 3, 2, 1)*1.5, col=rgb(0,0,0, 0.9),
         y.intersp=2.5, x.intersp=2, cex=1.05, xpd=TRUE)
  legend("bottomleft", legend=c("COUNT", "MANSION", "HOLDING"), pch=c(17, 15, 19), 
         col=c("red", "blue", "dimgrey"), text.col=c("red", "blue", "dimgrey"),
         bty="n", cex=0.75)
}

# Close the PNG file output
dev.off()

# Bind the PNG files into a Giff
dosPath <- "C:\\Users\\Joseph Crisp\\Desktop\\UbuntuSharedFolder\\Woodchester_CattleAndBadgers\\NewAnalyses_13-07-17\\CattleTestData\\"
system(paste("magick -delay 120 ", '\"', dosPath, "DynamicGiff/BreakdownSampling_*.png\" \"",
             dosPath, "DynamicGiff/BreakdownSampling.gif\"", sep=""))

# Delete the PNG files
unlink(paste(path, "CattleTestData/DynamicGiff/BreakdownSampling_*.png", sep=""))

#############
# FUNCTIONS #
#############

plotCountsThroughTime <- function(column, columnInCounts, countsForTimestep, log,
                                  breakdownsTable, title, polygonColour, pointsCol,
                                  pointsSize, lineWidth, lineCol){
  
  # Build the column names
  upperCol <- paste(columnInCounts, "Upper", sep="")
  medianCol <- paste(columnInCounts, "Median", sep="")
  lowerCol <- paste(columnInCounts, "Lower", sep="")
  
  # Find the y axis max
  yMax <- max(countsForTimestep[, upperCol], na.rm=TRUE)
  
  # Get the y values
  yValues <- countsForTimestep[, medianCol]

  # Plot the data
  
  cexAxis <- 2
  cexMain <- 3
  cexLab <- 2
  
  if(log == TRUE){
    yValues <- log(yValues)
    yMax <- log(yMax)
    plot(countsForTimestep$Date, yValues, type="l",
         pch=19, col="white", yaxt="n", cex=1, cex.axis=cexAxis, cex.lab=cexLab, cex.main=cexMain, 
         lwd=lineWidth, ylab="Frequency", xlab="Year", main=title, ylim=c(0, yMax))
    axis(side=2, at=log(c(1,10,100,200)), labels=c(1,10,100,200), las=1, cex.axis=cexAxis)
  }else{
    plot(countsForTimestep$Date, yValues, type="l",
         pch=19, col="white", cex=1, cex.axis=cexAxis, cex.lab=cexLab, cex.main=cexMain,
         lwd=lineWidth, ylab="Frequency", xlab="Year", main=title, ylim=c(0, yMax), las=1)
  }

  polygonValues <- c(countsForTimestep[, upperCol], rev(countsForTimestep[, lowerCol]))
  if(log == TRUE){
    polygonValues <- log(polygonValues)
    polygonValues[is.infinite(polygonValues)] <- -1
  }
  polygonValues <- informativelyReplaceNAs(polygonValues)
  
  polygon(y=polygonValues, 
          x=c(countsForTimestep$Date, rev(countsForTimestep$Date)),
          col=polygonColour, border=NA)
  
  yValues <- breakdownsTable[, column]
  if(log == TRUE){
    yValues <- log(yValues)
  }
  points(breakdownsTable$Date, yValues,
         pch=19, col=pointsCol, cex=pointsSize)
  
  # Re-plot moving average line
  yValues <- countsForTimestep[, medianCol]
  if(log == TRUE){
    yValues <- log(yValues)
  }
  points(countsForTimestep$Date, yValues, type="l", col=lineCol, lwd=lineWidth)
}

getReactorCultureAndSequencedCounts <- function(dateRange, timestep, breakdownsTable){
  
  # Get a range of dates from the start date - timestep to the end date
  dates <- c(seq(dateRange[1], dateRange[2], timestep), dateRange[2])
  dates <- c(seq(dateRange[1], length=2, by=paste("-", timestep, sep=""))[1], dates)
  
  table <- data.frame(Date=dates[-1], 
                      NReactorsLower=rep(0, length(dates) - 1),
                      NReactorsMedian=rep(0, length(dates) - 1),
                      NReactorsUpper=rep(0, length(dates) - 1),
                      NCulturedLower=rep(0, length(dates) - 1),
                      NCulturedMedian=rep(0, length(dates) - 1),
                      NCulturedUpper=rep(0, length(dates) - 1),
                      NSequencedLower=rep(0, length(dates) - 1),
                      NSequencedMedian=rep(0, length(dates) - 1),
                      NSequencedUpper=rep(0, length(dates) - 1))
  
  for(i in 2:length(dates)){
    
    subsetted <- breakdownsTable[breakdownsTable$Date > dates[i-1] &
                                   breakdownsTable$Date <= dates[i], ]
    
    median <- median(subsetted$cbNumCattleReactors, na.rm=TRUE)
    range <- range(subsetted$cbNumCattleReactors, na.rm=TRUE)
    table[i-1, "NReactorsLower"] <- range[1]
    table[i-1, "NReactorsMedian"] <- median
    table[i-1, "NReactorsUpper"] <- range[2]
    
    median <- median(subsetted$cbNumCattleCultPos, na.rm=TRUE)
    range <- range(subsetted$cbNumCattleCultPos, na.rm=TRUE)
    table[i-1, "NCulturedLower"] <- range[1]
    table[i-1, "NCulturedMedian"] <- median
    table[i-1, "NCulturedUpper"] <- range[2]
    
    median <- median(subsetted$NSequencedIsolates, na.rm=TRUE)
    range <- range(subsetted$NSequencedIsolates, na.rm=TRUE)
    table[i-1, "NSequencedLower"] <- range[1]
    table[i-1, "NSequencedMedian"] <- median
    table[i-1, "NSequencedUpper"] <- range[2]
  }

  # Replace Inf values
  is.na(table) <- sapply(table, is.infinite)
  
  return(table)
}

informativelyReplaceNAs <- function(vector){
  
  for(i in 1:length(vector)){
    
    if(is.na(vector[i]) == TRUE){
      
      # Does a previous and succeeding value exist?
      if(i == 1){
        vector[i] <- returnNextNonNA(vector, i)
      }else if(i == length(vector)){
        vector[i] <- vector[i-1]
      }else{
        vector[i] <- (vector[i-1] + returnNextNonNA(vector, i)) / 2
      }
    }
  }
  
  return(vector)
}

returnNextNonNA <- function(vector, index){
  
  value <- NA
  for(i in index+1:length(vector)){
    if(is.na(vector[i]) == FALSE){
      value <- vector[i]
      break
    }
  }
  
  return(value)
}

removeBreakdownsThatOccurOutsideOfSamplingWindow <- function(breakdownsTable, dateRange){
  
  rowsToKeep <- c()
  index <- 0
  
  breakdownsTable$Date <- rep("2016-11-22", nrow(breakdownsTable))
  
  for(row in 1:nrow(breakdownsTable)){
    
    # Get the breakdown date
    breakdownDate <- as.Date(strsplit(breakdownsTable[row, "cbBreakId"], split="-")[[1]][2],
                             "%d/%m/%Y")
    breakdownsTable[row, "Date"] <- as.character(breakdownDate)
    
    # Note row if date is within date range
    if(breakdownDate - dateRange[1] >= 0 && breakdownDate - dateRange[2] <= 0){
      index <- index + 1
      rowsToKeep[index] <- row
    }
  }
  
  return(breakdownsTable[rowsToKeep, ])
}

getSampledBreakdownDateRange <- function(samplingInfo){
  
  dateRange <- as.Date(c("2016-11-22", "1800-01-22"))
  for(row in 1:nrow(samplingInfo)){
    
    # Get the breakdown date
    breakdownDate <- as.Date(strsplit(samplingInfo[row, "BreakdownID"], split="-")[[1]][2],
                             "%d/%m/%Y")
    
    # Is the breakdown date less than minDate?
    if(breakdownDate - dateRange[1] < 0){
      dateRange[1] <- breakdownDate
    }
    
    # Is the breakdown date after the maxDate?
    if(breakdownDate - dateRange[2] > 0){
      dateRange[2] <- breakdownDate
    }
  }
  
  return(dateRange)
}

keepRowsIfColumnValuePresentInList <- function(column, table, listStructure){
  
  # Initialise an array to store the rows to keep
  rowsToKeep <- c()
  index <- 0
  
  # Examine each row of table
  for(row in 1:nrow(table)){
    
    # Keep row if column value is in list
    if(is.null(listStructure[[as.character(table[row, column])]]) == FALSE){
      index <- index + 1
      rowsToKeep[index] <- row
    }
  }
  
  return(table[rowsToKeep, ])
}

convertTableToList <- function(column, table){
  
  tableAsList <- list()
  for(row in 1:nrow(table)){
    
    tableAsList[[table[row, column]]] <- table[row, ]
  }
  
  return(tableAsList)
}

formatCPHForLocationInfo <- function(cphs){
  
  formattedCPHs <- c()
  for(i in 1:length(cphs)){
    
    formattedCPHs[i] <- paste(strsplit(
      strsplit(cphs[i], split="-")[[1]][1], split="/")[[1]], collapse="")
  }
  
  return(formattedCPHs)
}

keepLocationsWithinThresholdDistance <- function(locationInfo, thresholdInMetres, 
                                                 badgerCentre){
  
  locationInfo$DistanceToWoodchester <- rep(0, nrow(locationInfo))
  
  # Initialise an array to record the rows to keep
  rowsToKeep <- c()
  index <- 0
  
  # Examine each location
  for(row in 1:nrow(locationInfo)){
    
    # Skip if no location information available
    if(is.na(locationInfo[row, "x"]) == TRUE || is.na(locationInfo[row, "y"]) == TRUE){
      next
    }
    
    # Calculate distance to Woodchester Mansion
    distance <- euclideanDistance(x1=badgerCentre[1], y1=badgerCentre[2], 
                                  x2=locationInfo[row, "x"], y2=locationInfo[row, "y"])
    locationInfo[row, "DistanceToWoodchester"] <- distance
    
    # Keep the row if distance is <=threshold
    if(distance <= thresholdInMetres){
      index <- index + 1
      rowsToKeep[index] <- row
    }
    
    if(row %% 10000 == 0){
      print(paste("Finished reading row ", row, " of ", nrow(locationInfo), sep=""))
    }
  }
  
  return(locationInfo[rowsToKeep, ])
}

euclideanDistance <- function(x1, y1, x2, y2){
  return(sqrt(sum((x1 - x2)^2 + (y1 - y2)^2)))
}
