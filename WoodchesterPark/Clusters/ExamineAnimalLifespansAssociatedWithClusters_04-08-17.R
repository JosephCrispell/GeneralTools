###################################
# Load and process the input data #
###################################

# Create a path variable
path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/InterSpeciesClusters/"

# Load the life history summaries
file <- paste(path, "sampledAnimalsLifeHistories_22-11-2018.txt", sep="")
table <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep="\t",
                    colClasses = "character")

# Create a table for each cluster
clusterTables <- splitInputTableIntoInfoForClusters(table)

#############################
# Plot the Animal Lifespans #
#############################

# Time Period
xLim <- as.Date(c("1988-01-01", "2015-12-31")) # year-month-day

# Create a condensed version of the above plot
file <- paste(path, "sampledAnimalLifespansInClusters-CONDENSED_22-11-18.pdf", sep="")
pdf(file, height=11, width=8)

# Get an array of all the clusters
clusters <- sort(names(clusterTables))[-length(names(clusterTables))] # Remove last cluster (all isolates)

for(cluster in clusters){
  
  plotAnimalLifespansForClusterSummariseUnSampled(clusterTables, cluster,
                                                  xLim)
}

dev.off()

#########################
# Functions - Lifespans #
#########################

roundToNearestMultiple <- function(value, multiple){
  
  output <- value
  while(output %% multiple != 0){
    output <- output + 1
  }
  
  return(output)
}

plotAnimalLifespansForClusterSummariseUnSampled <- function(clusterTables, cluster,
                                                            xLim){

  # Note the dates within range of interest
  datesInRange <- seq(xLim[1], xLim[2], 1)
  indexedDates <- indexArrayOfDates(datesInRange)
  
  # Get the life history information for the current cluster
  clusterTable <- clusterTables[[cluster]]
  
  # Get only sampled animal life histories
  sampled <- orderTableByStartDate(clusterTable[is.na(clusterTable$Isolates) == FALSE, ])

  # Get not sampled and split by their testing
  unSampledAnimalTables <- getUnSampledAndSplitByTestHistory(clusterTable)
  unSampledCattleReactors <- orderTableByStartDate(unSampledAnimalTables[["Cattle - reactors"]])
  unSampledCattleNegatives <- orderTableByStartDate(unSampledAnimalTables[["Cattle - negatives"]])
  unSampledBadgerReactors <- orderTableByStartDate(unSampledAnimalTables[["Badger - reactors"]])
  unSampledBadgerNegatives <- orderTableByStartDate(unSampledAnimalTables[["Badger - negatives"]])

  # Summarise the unsampled cattle that encountered the sampled cattle
  unSampledCattle <- combineTwoTables(unSampledCattleReactors, unSampledCattleNegatives)
  cattleCounts <- countNegativeAndTestReactingCattleDuringPeriod(unSampledCattle,
                                                                 datesInRange,
                                                                 indexedDates)
  
  # Summarise the unsampled badgers that encountered the sampled badgers
  unSampledBadgers <- combineTwoTables(unSampledBadgerReactors, unSampledBadgerNegatives)
  badgerCounts <- countNegativeAndPositiveBadgersDuringPeriod(unSampledBadgers, datesInRange,
                                                              indexedDates)
  
  # Check if in-contact animals are available
  unSampledBadgersPresent <- 
    nrow(unSampledBadgers) > 0 && is.na(unSampledBadgers[1, 1]) == FALSE
  unSampledCattlePresent <- 
    nrow(unSampledCattle) > 0 && is.na(unSampledCattle[1, 1]) == FALSE
  
  if(unSampledBadgersPresent == TRUE && unSampledCattlePresent == TRUE){
    
    # Set the layout of the plots within the plotting window
    layout(matrix(c(1, 2, 3, 3), nrow=4, ncol=1, byrow=TRUE))
    if((as.numeric(cluster) + 1) == 4){
      layout(matrix(c(1, 2, 3, 3, 3, 3), nrow=6, ncol=1, byrow=TRUE))
    }
    
    # Plot the badger counts
    plotUnSampledBadgerCounts(badgerCounts, datesInRange, cluster, TRUE)
    
    # Plot the cattle counts
    plotUnSampledCattleCounts(cattleCounts, datesInRange, cluster, FALSE)
    
    # Plot the lifespans of the sampled animals
    plotSampledAnimalLifespans(sampled, datesInRange, cluster, FALSE)
    
    if((as.numeric(cluster) + 1) %in% c(2, 4)){

      par(mfrow=c(1,1))
      # Plot the lifespans of the sampled animals
      plotSampledAnimalLifespans(sampled, datesInRange, cluster, TRUE)
    }

  }else if(unSampledBadgersPresent == TRUE && unSampledCattlePresent == FALSE){
    
    par(mfrow=c(2,1))
    
    # Plot the badger counts
    plotUnSampledBadgerCounts(badgerCounts, datesInRange, cluster, TRUE)

    # Plot the lifespans of the sampled animals
    plotSampledAnimalLifespans(sampled, datesInRange, cluster, FALSE)
    
  }else if(unSampledBadgersPresent == FALSE && unSampledCattlePresent == TRUE){
    
    par(mfrow=c(2,1))

    # Plot the cattle counts
    plotUnSampledCattleCounts(cattleCounts, datesInRange, cluster, TRUE)
    
    # Plot the lifespans of the sampled animals
    plotSampledAnimalLifespans(sampled, datesInRange, cluster, FALSE)
  }else{
    
    par(mfrow=c(1,1))

    # Plot the lifespans of the sampled animals
    plotSampledAnimalLifespans(sampled, datesInRange, cluster, TRUE)
  }
}

plotSampledAnimalLifespans <- function(sampled, datesInRange, cluster,
                                       addCluster){
  
  # Set the plotting margins
  par(mar = c(4.1,4.1,1.2,4.1)) # bottom left top right
  if(addCluster == TRUE){
    par(mar = c(4.1,4.1,3.1,4.1)) # bottom left top right
  }
  
  # Note number of sampled animals
  nAnimals <- nrow(sampled)
  
  # Create an empty plot
  plot(x = xLim, y=rep(1, 2), type="n",
       ylim=c(1, nAnimals),
       yaxt="n", ylab="", xlab="Year",
       main="")
  
  # Plot title
  mtext("Lifespans of Sampled Animals", side=3)
  
  # Add the animal lifespans
  addPointsForAnimals(sampled, cluster)
  
  # Add Species labels on Y axis
  axis(side=2, at=seq(1, nrow(sampled)), labels=FALSE)
  xLabPosition <- xLim[1] - (0.05 * (xLim[2] - xLim[1]))
  text(labels=sampled$Species, 
       col=ifelse(sampled$Species == "BADGER", rgb(0,0,0, 0.5), "black"),
       x=rep(xLabPosition,length(sampled$Species)),
       y=1:length(sampled$Species),
       srt = 0, pos = 2, xpd = TRUE, cex=0.7)
  
  # Add lines to highlight years
  addLinesForYears(datesInRange)
  
  # Add Legend
  if(cluster == 3){
    legend("topleft", legend=c("Isolate Obtained", "Detection/Breakdown",
                               "Positive Test"),
           pch=c(18, 16, 16, 16), 
           col=c("blue", "cyan", "red"),
           text.col=c("blue", "cyan", "red"),
           bty="n", cex=1)
  }else{
    legend("topleft", legend=c("Isolate Obtained", "Detection/Breakdown",
                               "Positive Test", "Inconclusive Test"),
           pch=c(18, 16, 16, 16, 16), 
           col=c("blue", "cyan", "red", "orange"),
           text.col=c("blue", "cyan", "red", "orange"),
           bty="n", cex=1)
  }
  
  
  # Add the cluster number to plot
  if(addCluster == TRUE){
    mtext(paste("Clade: ", (as.numeric(cluster) + 1)), side=3, 
          at=as.Date(paste(format(datesInRange[length(datesInRange)], "%Y"), "-01-01", sep="")),
          line=1.5, font=2)
  }
}

addLinesForYears <- function(datesInRange){
  years <- seq(from=format(datesInRange[1], "%Y"),
               to=format(datesInRange[length(datesInRange)], "%Y"), by=5)
  for(i in 1:length(years)){
    
    year <- as.Date(paste(years[i],"01", "01", sep="-"))
    
    abline(v=year, col=rgb(0,0,0, 0.3), lty=2)
  }
}

plotUnSampledBadgerCounts <- function(badgerCounts, datesInRange, cluster,
                                      addCluster){
  
  # Set the plot margins
  par(mar = c(0.1,4.1,3.1,4.1)) # bottom left top right
  
  # Get the counts
  negativeCounts <- badgerCounts[["Negative"]]
  positiveCounts <- badgerCounts[["Positive"]]
  
  # Set the y axis limits
  yLim <- c(0,max(c(negativeCounts + positiveCounts)))
  
  # Plot the number of animals reacting to the testing
  plot(x = datesInRange, y=negativeCounts, type="n",
       ylim=yLim, ylab="Number Badgers Present", main="",
       xaxt="n", xlab="", las=1)
  
  # Add the points for present and positive badgers
  points(x=datesInRange, y=negativeCounts + positiveCounts, type="l", col=rgb(0,0,0, alpha=0.5))
  points(x=datesInRange, y=positiveCounts, type="l", col=rgb(1,0,0, alpha=0.75))
  
  # Plot title
  mtext("Number of Badgers Encountered by Sampled Badgers",
        side=3)
  
  # Add lines to highlight years
  addLinesForYears(datesInRange)
  
  # Add a legend
  legend("topleft", legend=c("Total", "Positive"), 
         text.col=c(rgb(0,0,0, alpha=0.5), rgb(1,0,0, alpha=1)), bty="n", cex=1)

  # Add the cluster number to plot
  if(addCluster == TRUE){
    mtext(paste("Clade: ", (as.numeric(cluster) + 1)), side=3, 
          at=as.Date(paste(format(datesInRange[length(datesInRange)], "%Y"), "-01-01", sep="")),
          line=1.5, font=2)
  }
}

plotUnSampledCattleCounts <- function(cattleCounts, datesInRange, cluster,
                                      addCluster){
  
  # Set the plot margins
  par(mar = c(0.1,4.1,1.2,4.1)) # bottom left top right

  # Get the counts
  negativeCounts <- cattleCounts[["Negative"]]
  inconclusiveCounts <- cattleCounts[["Inconclusive"]]
  positiveCounts <- cattleCounts[["Positive"]]
  
  # Set the y axis limits
  yLim <- c(0,max(c(inconclusiveCounts, positiveCounts)))
  
  # Plot the number of animals reacting to the testing
  plot(x = datesInRange, y=inconclusiveCounts, type="n",
       ylim=yLim, ylab="Number Animals Reacting to SICCT", main="",
       xaxt="n", xlab="", las=1)

  # Add the lines to indicate number of animals encountered that are present
  present <- negativeCounts + inconclusiveCounts + positiveCounts
  rightAxisLimit <- roundToNearestMultiple(max(present), length(seq(yLim[1], yLim[2], yLim[2]/5)) - 1)
  presentScaled <- (((present - min(present)) * (yLim[2] - yLim[1])) / 
                      (rightAxisLimit - min(present))) + yLim[1]
  points(x=datesInRange, y=presentScaled, type="l", col=rgb(0,0,0, alpha=0.5))
  
  # Add axis on right
  at <- seq(yLim[1], yLim[2], yLim[2]/5)
  labels <- seq(min(present), rightAxisLimit, 
                (rightAxisLimit - min(present))/(length(at) - 1))
  axis(side=4, at=at, las=1,labels=round(labels, digits=0))
  mtext(side=4, text="Number Animals Present", line=3, cex=0.75)
  
  legend("topright", legend="Total", 
         text.col=rgb(0,0,0, 0.5), bty="n", cex=1)
  
  # Add the points for test reactors
  points(x=datesInRange, y=inconclusiveCounts, type="l", col=rgb(0,0,1, alpha=0.75))
  points(x=datesInRange, y=positiveCounts, type="l", col=rgb(1,0,0, alpha=0.75))
  
  # Plot title
  mtext("Number of Cattle Encountered by Sampled Cattle",
        side=3)
  
  # Add lines to highlight years
  addLinesForYears(datesInRange)
  
  # Add a legend
  legend("topleft", legend=c("Inconclusive", "Positive"), 
         text.col=c("blue", "red"), bty="n", cex=1)
  
  # Add the cluster number to plot
  if(addCluster == TRUE){
    mtext(paste("Clade: ", (as.numeric(cluster) + 1)), side=3, 
          at=as.Date(paste(format(datesInRange[length(datesInRange)], "%Y"), "-01-01", sep="")),
          line=1.5, font=2)
  }
  
}

countNegativeAndPositiveBadgersDuringPeriod <- function(unSampled, datesInRange, indexedDates){
  
  negativeCounts <- rep(0, length(datesInRange))
  positiveCounts <- rep(0, length(datesInRange))
  
  for(row in 1:nrow(unSampled)){
    
    # Indicate progress
    if(row %% 100 == 0){
      cat(paste("\rReading row ", row, " of ", nrow(unSampled), sep=""))
    }
    
    # Skip animals for which there are no movements
    if(is.na(unSampled[row, "MovementDates"]) == TRUE){
      next
    }
    
    # Get the animal movement dates
    movementDates <- as.Date(strsplit(
      as.character(unSampled[row, "MovementDates"]), split=",")[[1]], "%d-%m-%Y")
    
    # Remove movement dates outside of period of interest
    movementDates <- movementDates[movementDates >= datesInRange[1]]
    movementDates <- movementDates[movementDates <= datesInRange[length(datesInRange)]]
    
    # Note all days that current animal alive for
    daysPresent <- seq(movementDates[1], movementDates[length(movementDates)], 1)
    for(i in 1:length(daysPresent)){
      negativeCounts[indexedDates[[as.character(daysPresent[i])]]] <- 
        negativeCounts[indexedDates[[as.character(daysPresent[i])]]] + 1
    }
    
    # Skip animals that don't have a detection date
    if(is.na(unSampled[row, "DetectionDate"]) == TRUE){
      next
    }
    
    # Examine the current cow's test history
    detectionDate <- as.Date(unSampled[row, "DetectionDate"], "%d-%m-%Y")
    
    if(is.null(indexedDates[[as.character(detectionDate)]]) == FALSE){
      
      startOfPeriodIndex <- indexedDates[[as.character(daysPresent[1])]]
      endOfPeriodIndex <- indexedDates[[as.character(daysPresent[length(daysPresent)])]]
      detectionDateIndex <- indexedDates[[as.character(detectionDate)]]
      
      positiveCounts[detectionDateIndex:endOfPeriodIndex] <- 
        positiveCounts[detectionDateIndex:endOfPeriodIndex] + 1
      negativeCounts[detectionDateIndex:endOfPeriodIndex] <- 
        negativeCounts[detectionDateIndex:endOfPeriodIndex] - 1
    
    # Check if badger infected before period of interest
    }else if(detectionDate < datesInRange[1]){
      
      startOfPeriodIndex <- indexedDates[[as.character(daysPresent[1])]]
      endOfPeriodIndex <- indexedDates[[as.character(daysPresent[length(daysPresent)])]]

      positiveCounts[startOfPeriodIndex:endOfPeriodIndex] <- 
        positiveCounts[startOfPeriodIndex:endOfPeriodIndex] + 1
      negativeCounts[startOfPeriodIndex:endOfPeriodIndex] <- 
        negativeCounts[startOfPeriodIndex:endOfPeriodIndex] - 1
    }
  }
  
  output <- list(
    "Negative" = negativeCounts,
    "Positive" = positiveCounts
  )
  
  return(output)
}

countNegativeAndTestReactingCattleDuringPeriod <- function(unSampled, datesInRange, indexedDates){
  
  negativeCounts <- rep(0, length(datesInRange))
  inconclusiveCounts <- rep(0, length(datesInRange))
  positiveCounts <- rep(0, length(datesInRange))
  
  for(row in 1:nrow(unSampled)){
    
    # Indicate progress
    if(row %% 100 == 0){
      cat(paste("\rReading row ", row, " of ", nrow(unSampled), sep=""))
    }
    
    # Skip animals for which there are no movements
    if(is.na(unSampled[row, "MovementDates"]) == TRUE){
      next
    }
    
    # Calculate the start and end of the current animals movements
    movementDates <- as.Date(strsplit(
      as.character(unSampled[row, "MovementDates"]), split=",")[[1]], "%d-%m-%Y")
    
    # Note all days that current animal alive for
    daysPresent <- seq(movementDates[1], movementDates[length(movementDates)], 1)
    for(i in 1:length(daysPresent)){
      negativeCounts[indexedDates[[as.character(daysPresent[i])]]] <- 
        negativeCounts[indexedDates[[as.character(daysPresent[i])]]] + 1
    }
    
    # Skip animals that don't have a test history - if looking at cattle
    if(is.na(unSampled[row, "CattleTestDates"]) == TRUE){
      next
    }
    
    # Examine the current cow's test history
    testDates <- as.Date(strsplit(
      as.character(unSampled[row, "CattleTestDates"]), split=",")[[1]], "%d-%m-%Y")
    testResults <- strsplit(as.character(unSampled[row, "CattleTestResults"]), split=",")[[1]]
    
    for(i in 1:length(testDates)){
      
      # Was the current test an Inconclusive?
      if(testResults[i] == "IR"){
        
        if(i < length(testDates)){
          inconclusiveDateIndex <- indexedDates[[as.character(testDates[i])]]
          endOfPeriodIndex <- indexedDates[[as.character(testDates[i+1])]]
          
          inconclusiveCounts[inconclusiveDateIndex:endOfPeriodIndex] <- 
            inconclusiveCounts[inconclusiveDateIndex:endOfPeriodIndex] + 1
          
          negativeCounts[inconclusiveDateIndex:endOfPeriodIndex] <- 
            negativeCounts[inconclusiveDateIndex:endOfPeriodIndex] - 1
        }else{
          inconclusiveCounts[indexedDates[[as.character(testDates[i])]]] <- 
            inconclusiveCounts[indexedDates[[as.character(testDates[i])]]] + 1
          
          negativeCounts[indexedDates[[as.character(testDates[i])]]] <- 
            negativeCounts[indexedDates[[as.character(testDates[i])]]] - 1
        }

      }else{
        positiveCounts[indexedDates[[as.character(testDates[i])]]] <- 
          positiveCounts[indexedDates[[as.character(testDates[i])]]] + 1
        
        negativeCounts[indexedDates[[as.character(testDates[i])]]] <- 
          negativeCounts[indexedDates[[as.character(testDates[i])]]] - 1
      }
    }
  }
  
  output <- list(
    "Negative" = negativeCounts,
    "Inconclusive" = inconclusiveCounts,
    "Positive" = positiveCounts
  )
  
  return(output)
}

addPointsForAnimals <- function(table, cluster){
  
  # Add lifespan information from each animal into current plotting window
  for(row in 1:nrow(table)){
    
    ## Plot start and end of animal lifespan - distinguish cattle and badgers
    movementDates <- as.Date(strsplit(
      as.character(table[row, "MovementDates"]), split=",")[[1]], "%d-%m-%Y")
    points(x=c(movementDates[1], movementDates[length(movementDates)]), y=c(row,row),
           type="l", lwd=4, 
           col=ifelse(table[row, "Species"] == "BADGER", "grey", "black"))

    ## Cattle Testing Reactors
    if(is.na(table[row, "CattleTestResults"]) == FALSE){
      
      testDates <- as.Date(strsplit(
        as.character(table[row, "CattleTestDates"]), split=",")[[1]], "%d-%m-%Y")
      testResults <- strsplit(as.character(table[row, "CattleTestResults"]), split=",")[[1]]
      
      # Examine each test
      for(i in 1:length(testDates)){
        
        # Inconclusive
        if(testResults[i] == "IR"){
          colour <- "orange"
          points(x=testDates[i], y=row, pch=16, col=colour)
          # Test Positive
        }else if(testResults[i] == "R" || testResults[i] == "SL"){
          colour <- "red"
          points(x=testDates[i], y=row, pch=16, col=colour)
        }
      }
    }
    
    ## Infection Detection
    if(is.na(table[row, "DetectionDate"]) == FALSE){
      
      points(x=as.Date(table[row, "DetectionDate"], "%d-%m-%Y"), y=row, pch=16,
             col="cyan")
    }
    
    ## Sampling Events - note those from different clusters
    if(is.na(table[row, "SamplingDates"]) == FALSE){
      
      # Get Sampling information
      samplingDates <-  testDates <- as.Date(strsplit(
        as.character(table[row, "SamplingDates"]), split=",")[[1]], "%d-%m-%Y")
      clusters <- strsplit(as.character(table[row, "Clusters"]), split=",")[[1]]

      # Add points of sampled isolates in current cluster
      for(i in 1:length(clusters)){
        
        if(clusters[i] == cluster){
          
          # Add point for isolate in cluster
          points(x=samplingDates[i], y=row, pch=18, col="blue")

        }else{
          
          # Add point for isolate that isn't in cluster
          points(x=samplingDates[i], y=row, pch=5, col="blue")
          
        }
      }
    }
  }
}

#######################
# Functions - General #
#######################

getUnSampledAndSplitByTestHistory <- function(clusterTable){
  
  # Get only the sampled animals
  unSampled <- clusterTable[is.na(clusterTable$Isolates) == TRUE, ]
  unSampledCattle <- unSampled[unSampled$Species == "COW", ]
  unSampledBadgers <- unSampled[unSampled$Species == "BADGER", ]
  
  # Split the unsampled animals by whether they were tested
  splitTables <- list()
  splitTables[["Cattle - reactors"]] <- unSampledCattle[is.na(unSampledCattle$CattleTestResults) == FALSE, ]
  splitTables[["Cattle - negatives"]] <- unSampledCattle[is.na(unSampledCattle$CattleTestResults) == TRUE, ]
  splitTables[["Badger - reactors"]] <- unSampledBadgers[is.na(unSampledBadgers$DetectionDate) == FALSE, ]
  splitTables[["Badger - negatives"]] <- unSampledBadgers[is.na(unSampledBadgers$DetectionDate) == TRUE, ]
  
  return(splitTables)
}

indexArrayOfDates <- function(array){
  output <- list()
  for(i in 1:length(array)){
    output[[as.character(array[i])]] <- i
  }
  
  return(output)
}

removeMovementsWhereLocationAreUnknownOrIfToSlaughter <- function(table, premisesToIgnore){
  
  # Examine each row of the table
  for(row in 1:nrow(table)){
    
    # Skip animals with no movement data
    if(is.na(table[row, "MovementDates"]) == TRUE){
      next
    }
    
    # Get the movement dates and locations
    movementDates <- as.Date(strsplit(
      as.character(table[row, "MovementDates"]), split=",")[[1]], "%d-%m-%Y")
    Xs <- as.numeric(strsplit(table[row, "Xs"], ",")[[1]])
    Ys <- as.numeric(strsplit(table[row, "Ys"], ",")[[1]])
    
    # Remove movements to premises we want to ignore - CATTLE
    if(table[row, "Species"] == "COW"){
      premisesTypes <- strsplit(table[row, "PremisesTypes"], ",")[[1]]
      
      keep <- c()
      for(i in 1:length(premisesTypes)){
        
        keep[i] <- TRUE
        if(is.null(premisesToIgnore[[premisesTypes[i]]]) == FALSE){
          keep[i] <- FALSE
        }
      }
      movementDates <- movementDates[keep]
      Xs <- Xs[keep]
      premisesTypes <- premisesTypes[keep]
      Ys <- Ys[keep]
    }else{
      premisesTypes <- rep("NA", length(movementDates))
    }
    
    # Remove info associated with NA locations
    movementDates <- movementDates[is.na(Ys) == FALSE]
    Xs <- Xs[is.na(Ys) == FALSE]
    premisesTypes <- premisesTypes[is.na(Ys) == FALSE]
    Ys <- Ys[is.na(Ys) == FALSE]
    
    # Order the movements by date
    order <- order(movementDates)
    movementDates <- movementDates[order]
    Xs <- Xs[order]
    Ys <- Ys[order]
    premisesTypes <- premisesTypes[order]
    
    # Convert Movement dates back into previous format
    movementDates <- as.character(movementDates, format="%d-%m-%Y")
    
    # Put these data back into table - make sure there are movements left
    if(length(movementDates) > 0){
      table[row, "MovementDates"] <- paste(movementDates, collapse=",")
      table[row, "Xs"] <- paste(Xs, collapse=",")
      table[row, "Ys"] <- paste(Ys, collapse=",")
      table[row, "PremisesTypes"] <- paste(premisesTypes, collapse=",")
    }else{
      table[row, "MovementDates"] <- NA
      table[row, "Xs"] <- NA
      table[row, "Ys"] <- NA
      table[row, "PremisesTypes"] <- NA
    }
  }
  
  return(table)
}

combineTwoTables <- function(tableA, tableB){
  
  # Calculate the number of rows in the combined table
  nRow <- nrow(tableA) + nrow(tableB)
  
  # Create a new empty table
  combinedTable <- as.data.frame(matrix(ncol=ncol(tableA)))
  colnames(combinedTable) <- colnames(tableA)
  
  # Fill the table
  for(row in 1:nrow(tableA)){
    combinedTable[row, ] <- tableA[row, ]
  }
  from <- nrow(tableA)
  if(nrow(tableB) > 0){
    for(row in 1:nrow(tableB)){
      combinedTable[row + from, ] <- tableB[row, ]
    }
    from <- from + nrow(tableB)
  }
  
  return(combinedTable)
}

orderTableByDetectionDate <- function(table){
  
  # Check that there are rows in the table
  if(nrow(table) > 0){
    # Initialise an array to store the start dates of each animals observed lifespan
    detectionDates <- c()
    
    # Make sure the MovementDates column is a character
    table[, "DetectionDate"] <- as.character(table[, "DetectionDate"])
    
    # Examine each of the animals
    for(row in 1:nrow(table)){
      
      # Get the movement dates for the current animal - check that they're available
      if(is.na(table[row , "DetectionDate"]) == FALSE){
        detectionDate <- as.Date(as.character(table[row, "DetectionDate"]), "%d-%m-%Y")
        
        # Store the start date
        detectionDates[row] <- detectionDate
      }else{
        detectionDates[row] <- NA
      }
    }
    
    # Order the table
    orderedTable <- table[order(detectionDates), ]
  }else{
    orderedTable <- table
  }
  
  # Return the table ordered by the start dates
  return(orderedTable)
}

orderTableByStartDate <- function(table){
  
  # Check that there are rows in the table
  if(nrow(table) > 0){
    # Initialise an array to store the start dates of each animals observed lifespan
    startDates <- c()
    
    # Make sure the MovementDates column is a character
    table[, "MovementDates"] <- as.character(table[, "MovementDates"])
    
    # Examine each of the animals
    for(row in 1:nrow(table)){
      
      # Get the movement dates for the current animal - check that they're available
      if(is.na(table[row , "MovementDates"]) == FALSE){
        movementDates <- as.Date(strsplit(
          as.character(table[row, "MovementDates"]), split=",")[[1]], "%d-%m-%Y")
        
        # Store the start date
        startDates[row] <- movementDates[1]
      }else if(is.na(table[row , "DetectionDate"]) == FALSE){
        
        # Store detection date if no life history information available
        startDates[row] <- as.Date(as.character(table[row, "DetectionDate"]),
                                   "%d-%m-%Y")
      }else{
        startDates[row] <- NA
      }
    }
    
    # Order the table
    orderedTable <- table[order(startDates), ]
  }else{
    orderedTable <- table
  }
  
  # Return the table ordered by the start dates
  return(orderedTable)
}

splitInputTableIntoInfoForClusters <- function(table){
  clusters <- list()
  
  for(row in 1:nrow(table)){
    
    clustersIsolatesFoundIn <- unique(strsplit(table[row, "Clusters"], split=",")[[1]])
    
    # Add the information for the current animal into the clusters its isolates were found in
    for(clusterID in clustersIsolatesFoundIn){
      
      # If haven't already encountered this cluster
      if(is.null(clusters[[clusterID]]) == TRUE){
        
        clusters[[clusterID]] <- matrix(ncol=ncol(table))
        colnames(clusters[[clusterID]]) <- colnames(table)
        
        clusters[[clusterID]][1, ] <- table[row, ]
      }else{
        
        clusters[[clusterID]] <- rbind(clusters[[clusterID]], table[row, ])
      }
    }
  }
  
  # Convert each of the tables to dataframes
  for(clusterID in names(clusters)){
    
    clusters[[clusterID]] <- as.data.frame(clusters[[clusterID]], 
                                           row.names = 1:nrow(clusters[[clusterID]]))
  }
  
  return(clusters)
}
