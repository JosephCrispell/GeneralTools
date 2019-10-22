###################################
# Load and process the input data #
###################################

# Create a path variable
path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/InterSpeciesClusters/"

# Load the life history summaries
file <- paste(path, "sampledAnimalsLifeHistories_05-04-2018.txt", sep="")
table <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep="\t",
                    colClasses = "character")

# Remove movements where location is unknown or if movement to premises we want to ignore
premisesToIgnore <- list("SR"=1, "SW"=1, "EX"=1, "CC"=1)
table <- suppressWarnings(
  removeMovementsWhereLocationAreUnknownOrIfToSlaughter(table, premisesToIgnore))

# Find unique territory centroids - only keep one representative within threshold distance
territoryCentroids <- findUniqueTerritoryCentroids(table, threshold=150)

# Create a table for each cluster
clusterTables <- splitInputTableIntoInfoForClusters(table)

###############################################
# Plot the animal movements in space - Static #
###############################################

# Set margin size and plotting window
par(mar=c(2.5,2.5,3,2.1)) # bottom, left, top and right

# Get an array of all the clusters
clusters <- sort(names(clusterTables))[-length(names(clusterTables))] # Remove last cluster - sampled animals

# Note the centre of the badger territories
badgerCentre <- c(381761.7, 200964.3)
expand <- c(15000, 5000)

# Open a pdf
file <- paste(path, "sampledAnimalMovementsInClusters_06-04-18.pdf", sep="")
pdf(file, height=7, width=14)

# Set plotting window dimensions
par(mfrow=c(1,2))

# Plot the animal movements on a static figure
for(cluster in clusters){
  
  plotMovementsOfAnimalsInCluster(cluster, clusterTables, expand)
}

dev.off()

###########################################
# Plot animal movements in space and time #
###########################################

# Plot the animal movements for each year
for(cluster in clusters){
  
  plotMovementsOfAnimalsInClusterInEachYear(cluster, clusterTables, expand, alpha=0.75)
}


#######################
# Functions - Dynamic #
#######################

plotMovementsOfAnimalsInClusterInEachYear <- function(cluster, clusterTables, expand,
                                                      alpha){
  
  
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
  
  # Create one large stacked table with each ordered subset
  clusterTable <- combineTablesIntoOne(unSampledCattleNegatives, 
                                       unSampledBadgerNegatives,
                                       unSampledCattleReactors,
                                       unSampledBadgerReactors, sampled)
  
  prefix <- paste(path, "MovementGiff/Movements_cluster-", cluster, "_%02d.png", sep="")
  png(file=prefix, height=480, width=960)
  
  par(mfrow=c(1,2))
  par(mar=c(2.5,2.5,3,2.1)) # bottom, left, top and right
  
  # Set year
  for(year in 1988:2015){
    
    # Create empty plot
    createEmptyPlot(badgerCentre, expand[1], cluster)
    mtext(side=3, text=year, line=0.5, at=badgerCentre[1] - expand[1], cex=2)
    
    # Plot the movements of the animals in current year
    plotMovementsOfAnimalsInClusterInCurrentYear(year, clusterTable, alpha)
    addLegend()
    
    # Create empty plot
    createEmptyPlot(badgerCentre, expand[2], cluster)

    # Plot the movements of the animals in current year
    plotMovementsOfAnimalsInClusterInCurrentYear(year, clusterTable, alpha)
  }
  
  # Close the PNG file output
  dev.off()
  
  # Bind the PNG files into a Giff
  dosPath <- "C:\\Users\\Joseph Crisp\\Desktop\\UbuntuSharedFolder\\Woodchester_CattleAndBadgers\\NewAnalyses_22-03-18\\InterSpeciesClusters\\"
  system(paste("magick -delay 80 ", '\"', dosPath, "MovementGiff\\Movements_cluster-",
               cluster, "_*.png\" \"", dosPath, "MovementGiff\\Movements_cluster-", cluster,
               ".gif\"", sep=""))
  
  # Delete the PNG files
  unlink(paste(path, "MovementGiff/Movements_cluster-", cluster, "_*.png", sep=""))
}

plotMovementsOfAnimalsInClusterInCurrentYear <- function(year, clusterTable, alpha){
  
  # Examine each animal in the table - add movements where appropriate
  for(row in 1:nrow(clusterTable)){
    
    # Skip animals for which there are no movements
    if(is.na(clusterTable[row, "MovementDates"]) == TRUE){
      next
    }
    
    # Plot the current animals movements for current year
    plotMovementsOfAnimalInCurrentYear(year, row, clusterTable, alpha)
  }
}

plotMovementsOfAnimalInCurrentYear <- function(year, row, clusterTable, alpha){
  
  # Get the movement information for the current animal
  movementDates <- as.Date(strsplit(
    as.character(clusterTable[row, "MovementDates"]), split=",")[[1]], "%d-%m-%Y")
  Xs <- as.numeric(strsplit(as.character(clusterTable[row, "Xs"]), split=",")[[1]])
  Ys <- as.numeric(strsplit(as.character(clusterTable[row, "Ys"]), split=",")[[1]])
  
  # Examine testing history
  detectionHistory <- noteDetectionDatesAndStatuses(row, clusterTable)
  
  # Create variables to store previous location information
  previousPchInfo <- list()
  previousMovementYear <- 1900
  
  # Create a variable to note whether point was added
  locationAddedInCurrentYear <- FALSE
  
  # Examine each movement
  for(i in 1:length(movementDates)){
    
    # Get current status shape info
    pchInfo <- returnCurrentPchRelatingToStatus(date=movementDates[i],
                                                detectionHistory=detectionHistory,
                                                unSampled=is.na(clusterTable[row, "Isolates"]),
                                                species=clusterTable[row, "Species"])
    
    # Get year of current movement
    movementYear <- as.numeric(format(movementDates[i],'%Y'))
    
    # Add point - if location or status has changed and movement date in current year
    if(movementYear == year && (
      i == 1 || previousPchInfo[["col"]] != pchInfo[["col"]] || 
      checkIfLocationChanged(i, Xs, Ys) == TRUE ||
      previousMovementYear != movementYear)){
      
      # Record that adding point
      locationAddedInCurrentYear <- TRUE
      
      # Add point
      points(x=Xs[i], y=Ys[i], pch=pchInfo[["pch"]], 
             col=setAlpha(pchInfo[["col"]], alpha),
             bg=setAlpha(pchInfo[["bg"]], alpha), cex=1)
      
      # Add arrow from previous movement location if necessary
      addArrow(movementIndex=i, Xs=Xs, Ys=Ys, previousPchInfo=previousPchInfo)
      
      # Add previous point if movement date is in current year and location/status changed
      if(i != 1 && previousMovementYear != year && 
         (previousPchInfo[["col"]] != pchInfo[["col"]] || checkIfLocationChanged(i, Xs, Ys) == TRUE)){
        points(x=Xs[i-1], y=Ys[i-1], pch=previousPchInfo[["pch"]], 
               col=setAlpha(previousPchInfo[["col"]], alpha),
               bg=setAlpha(previousPchInfo[["bg"]], alpha), cex=1)
      }
    }
    
    # Store current pch info as previous
    previousPchInfo <- pchInfo
    previousMovementYear <- movementYear
  }
  
  # If no location added for current year - were there movements that spanned current year?
  if(locationAddedInCurrentYear == FALSE && length(movementDates) > 1){
    
    # Create variable to note whether found a movement in previous years
    foundMovementPriorToCurrentYear <- FALSE
    
    # Find last movement before current year
    for(i in 1:length(movementDates)){
      
      # Get year of current movement
      movementYear <- as.numeric(format(movementDates[i],'%Y'))
      
      if(movementYear < year){
        foundMovementPriorToCurrentYear == TRUE
      }else if(foundMovementPriorToCurrentYear == TRUE && movementYear > year){
        
        ## Add point from previous movement
        
        # Get current status shape info
        previousPchInfo <- returnCurrentPchRelatingToStatus(date=movementDates[i-1],
                                                            detectionHistory=detectionHistory,
                                                            unSampled=is.na(clusterTable[row, "Isolates"]),
                                                            species=clusterTable[row, "Species"])
        # Add point
        points(x=Xs[i-1], y=Ys[i-1], pch=previousPchInfo[["pch"]], 
               col=setAlpha(previousPchInfo[["col"]], alpha),
               bg=setAlpha(previousPchInfo[["bg"]], alpha), cex=1)
      }
    }
  }
}

######################
# Functions - Static #
######################

addLegend <- function(){
  
  # Statuses:
  #   Cow <- 24
  #   Badger <- 21
  #   Sampled <- bg="blue"
  #   Unsampled <- bg="grey"
  #   Inconclusive <- col="orange"
  #   Reactor <- col="red"
  
  legend("topleft", 
         legend=c("Badger", "Cow", "Sampled", "UnSampled", "Inconclusive", "Positive"),
         pch=c(21, 24, 21, 21, 21, 21),
         text.col=c("black", "black", "blue", "grey", "orange", "red"),
         col=c("black", "black", "blue", "grey", "orange", "red"),
         pt.bg=c("white", "white", "blue", "grey", "white", "white"),
         bty="o", cex=0.75, bg="white", box.col="white")
}

plotMovementsOfAnimalsInCluster <- function(cluster, clusterTables, expand){
  
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
  
  # Create one large stacked table with each ordered subset
  clusterTable <- combineTablesIntoOne(unSampledCattleNegatives, 
                                       unSampledBadgerNegatives,
                                       unSampledCattleReactors,
                                       unSampledBadgerReactors, sampled)
  
  ### Large spatial range
  
  # Create empty plot
  createEmptyPlot(badgerCentre, expand[1], cluster)
  
  # Add badger social group territory centroids
  # points(x=territoryCentroids$X, y=territoryCentroids$Y, pch=1,
  #        col=rgb(0,0,0, 0.25))
  
  # Examine each animal in the table - add movements where appropriate
  for(row in 1:nrow(clusterTable)){
    
    # Skip animals for which there are no movements
    if(is.na(clusterTable[row, "MovementDates"]) == TRUE){
      next
    }

    # Plot the animals movements
    plotMovementsOfAnimal(row, clusterTable, alpha=0.5)
  }
  addLegend()

  ### Small spatial range
  
  # Create empty plot
  createEmptyPlot(badgerCentre, expand[2], cluster)
  
  # Add badger social group territory centroids
  # points(x=territoryCentroids$X, y=territoryCentroids$Y, pch=1,
  #        col=rgb(0,0,0, 0.25))
  
  # Examine each animal in the table - add movements where appropriate
  for(row in 1:nrow(clusterTable)){
    
    # Skip animals for which there are no movements
    if(is.na(clusterTable[row, "MovementDates"]) == TRUE){
      next
    }
    
    # Plot the animals movements
    plotMovementsOfAnimal(row, clusterTable, alpha=0.5)
  }
}

plotMovementsOfAnimal <- function(row, clusterTable, alpha){
  
  # Get the movement information for the current animal
  movementDates <- as.Date(strsplit(
    as.character(clusterTable[row, "MovementDates"]), split=",")[[1]], "%d-%m-%Y")
  Xs <- as.numeric(strsplit(as.character(clusterTable[row, "Xs"]), split=",")[[1]])
  Ys <- as.numeric(strsplit(as.character(clusterTable[row, "Ys"]), split=",")[[1]])
  
  # Examine testing history
  detectionHistory <- noteDetectionDatesAndStatuses(row, clusterTable)
  
  # Plot movements
  previousPchInfo <- list()
  for(i in 1:length(movementDates)){
    
    # Get current status shape info
    pchInfo <- returnCurrentPchRelatingToStatus(date=movementDates[i],
                                                detectionHistory=detectionHistory,
                                                unSampled=is.na(clusterTable[row, "Isolates"]),
                                                species=clusterTable[row, "Species"])
    # Add point - if location or status has changed
    if(i == 1 || previousPchInfo[["col"]] != pchInfo[["col"]] || 
       checkIfLocationChanged(i, Xs, Ys) == TRUE){

      points(x=Xs[i], y=Ys[i], pch=pchInfo[["pch"]], 
             col=setAlpha(pchInfo[["col"]], alpha),
             bg=setAlpha(pchInfo[["bg"]], alpha), cex=1)
      
      # Add arrow from previous movement location if necessary
      addArrow(movementIndex=i, Xs=Xs, Ys=Ys, previousPchInfo=previousPchInfo)
    }
    
    # Store current pch info as previous
    previousPchInfo <- pchInfo
  }
}

checkIfLocationChanged <- function(movementIndex, Xs, Ys){
  result <- FALSE
  
  # Compare current and previous position
  previous <- paste(Xs[movementIndex-1], Ys[movementIndex-1], sep=":")
  current <- paste(Xs[movementIndex], Ys[movementIndex], sep=":")
  if(previous != current){
    result <- TRUE
  }
  
  return(result)
}

addArrow <- function(movementIndex, Xs, Ys, previousPchInfo){
  
  # Check if movement Index isn't 1
  # Check if X and Y coordinates changed from previous to current position
  if(movementIndex != 1 && checkIfLocationChanged(movementIndex, Xs, Ys) == TRUE){

    # Add Arrow
    arrows(x0=Xs[movementIndex-1], y0=Ys[movementIndex-1],
             x1=Xs[movementIndex], y1=Ys[movementIndex],
             angle=20, length=0.1, col=previousPchInfo[["col"]])
  }
}

returnCurrentPchRelatingToStatus <- function(date, detectionHistory, unSampled, species){
  
  # Statuses:
  #   Cow <- 24
  #   Badger <- 21
  #   Sampled <- bg="blue"
  #   Unsampled <- bg="grey"
  #   Inconclusive <- col="orange"
  #   Reactor <- col="red"
  
  # Initialise a list to store the point information
  pchInfo <- list(
    "pch" = 24,
    "bg" = "grey",
    "col" = "grey"
  )
  
  # Check if badger
  if(species == "BADGER"){
    pchInfo[["pch"]] <- 21
  }
  
  # Check if sampled
  if(unSampled == FALSE){
    pchInfo[["bg"]] <- "blue"
    pchInfo[["col"]] <- "blue"
  }
  
  # Check status
  status <- returnStatus(date, detectionHistory)
  if(status == "IR"){
    pchInfo[["col"]] <- "orange"
  }else if(status == "R"){
    pchInfo[["col"]] <- "red"
  }
  
  return(pchInfo)
}

returnStatus <- function(date, detectionHistory){
  
  status <- "S" # Susceptible
  
  if(detectionHistory[["Reactor"]] == TRUE){
    
    # Get detection dates
    dates <- detectionHistory[["Dates"]]
    statuses <- detectionHistory[["Statuses"]]
    
    # Check if current date is after any of above dates
    for(i in 1:length(dates)){
      if(date > dates[i]){
        status <- statuses[i]
      }
    }
  }
  
  return(status)
}

noteDetectionDatesAndStatuses <- function(row, clusterTable){
  
  # Badgers have single infected status <- "R" (Infected)
  # Cattle have two states <- "IR" (Inconclusive) and "R" (Infected)
  
  # Initialise an array to store status dates and statuses
  detectionDates <- c()
  detectionStatuses <- c()
  
  # Check if badger
  if(clusterTable[row, "Species"] == "BADGER"){
    
    # Where present store detection date
    if(is.na(clusterTable[row, "DetectionDate"]) == FALSE){
      detectionDates[length(detectionDates) + 1] <-
        clusterTable[row, "DetectionDate"]
      detectionStatuses[length(detectionStatuses) + 1] <- "R"
    }

  }else{
    
    # Examine cattle with testing histories
    if(is.na(clusterTable[row, "CattleTestDates"]) == FALSE){
      
      # Get test history
      testDates <- strsplit(
        as.character(clusterTable[row, "CattleTestDates"]), split=",")[[1]]
      testResults <- strsplit(as.character(clusterTable[row, "CattleTestResults"]),
                              split=",")[[1]]
      
      # Examine test history
      for(i in 1:length(testDates)){
        
        # Skip non-reactors/inconclusives
        if(testResults[i] %in% c("SL", "R", "IR") == FALSE){
          next
        }
        
        # Store detection date and status
        if(testResults[i] == "IR"){
          detectionDates[length(detectionDates) + 1] <- testDates[i]
          detectionStatuses[length(detectionStatuses) + 1] <- "IR"
        }else{
          detectionDates[length(detectionDates) + 1] <- testDates[i]
          detectionStatuses[length(detectionStatuses) + 1] <- "R"
        }
      }
    }
  }
  
  # Return detection history
  output <- list()
  if(length(detectionDates) > 0){
    output[["Dates"]] <- as.Date(detectionDates, "%d-%m-%Y")
    output[["Statuses"]] <- detectionStatuses
    output[["Reactor"]] <- TRUE
  }else{
    output[["Reactor"]] <- FALSE
  }
  
  return(output)
}

euclideanDistance <- function(x1, y1, x2, y2){
  
  return(sqrt(sum((x1 - x2)^2 + (y1 - y2)^2)))
}

findUniqueTerritoryCentroids <- function(table, threshold){
  
  # Create an array of the territory centroid coordinates
  territoryXs <- c()
  territoryYs <- c()
  
  # Initialise a list to note the centroids already examined
  checked <- list()
  
  # Examine each life history
  for(row in 1:nrow(table)){
    
    # Skip cattle
    if(table[row, "Species"] == "COW"){
      next
    }
    
    # Get X and Y coordinates
    Xs <- as.numeric(strsplit(table[row, "Xs"], split=",")[[1]])
    Ys <- as.numeric(strsplit(table[row, "Ys"], split=",")[[1]])
    
    # Examine the coordinates
    for(i in 1:length(Xs)){
      
      # Skip NAs
      if(is.na(Xs[i]) == TRUE || is.na(Ys[i]) == TRUE){
        next
      }
      
      # Create key
      key <- paste(Xs[i], Ys[i], sep=":")
      
      # Skip keys we have already checked
      if(is.null(checked[[key]]) == FALSE){
        next
      }
      
      # Note that checking key
      checked[[key]] <- 1
      
      # Skip coordinates that are close to those already added
      close <- FALSE
      if(length(territoryXs) > 0){
        for(pos in 1:length(territoryXs)){
          
          distance <- euclideanDistance(Xs[i], Ys[i],
                                        territoryXs[pos], territoryYs[pos])
          if(distance < threshold){
            close <- TRUE
            break
          }
        }
      }
      
      # Add coordinates that aren't close to those already added
      if(close == FALSE){
        territoryXs[length(territoryXs) + 1] <- Xs[i]
        territoryYs[length(territoryYs) + 1] <- Ys[i]
      }
    }
  }
  
  # Create table of Xs and Ys of unique territory centroids
  output <- data.frame(X=territoryXs, Y=territoryYs, stringsAsFactors=FALSE)
  
  return(output)
}

createEmptyPlot <- function(badgerCentre, expand, cluster){
  # Create empty plot
  plot(1, type="n", yaxt="n", xaxt="n",
       xlim=c(badgerCentre[1] - expand, badgerCentre[1] + expand),
       ylim=c(badgerCentre[2] - expand, badgerCentre[2] + expand),
       xlab="", ylab="", bty="n",
       main=paste("Cluster ", cluster, sep=""))
  
  # Add axis labels
  scale <- paste((expand * 2) / 1000, "km")
  mtext(text=scale, side=1, line=1)
  mtext(text=scale, side=2, line=1)
}

#######################
# Functions - General #
#######################

setAlpha <- function(colour, alpha){
  
  rgbValues <- col2rgb(colour)
  
  rgbColour <- rgb(rgbValues["red", 1], rgbValues["green", 1], rgbValues["blue", 1],
                   alpha=alpha*255, maxColorValue=255)
  
  return(rgbColour)
}

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

combineTablesIntoOne <- function(a,b,c,d,e){
  
  # Calculate the number of rows in the combined table
  nRow <- nrow(a) + nrow(b) + nrow(c) + nrow(d) + nrow(e)
  
  # Create a new empty table
  combinedTable <- as.data.frame(matrix(ncol=ncol(e)))
  colnames(combinedTable) <- colnames(e)
  
  # Fill the table
  from <- 0
  table <- a
  if(nrow(table) > 0){
    for(row in 1:nrow(table)){
      combinedTable[row + from, ] <- table[row, ]
    }
    from <- from + nrow(table)
  }
  table <- b
  if(nrow(table) > 0){
    for(row in 1:nrow(table)){
      combinedTable[row + from, ] <- table[row, ]
    }
    from <- from + nrow(table)
  }
  table <- c
  if(nrow(table) > 0){
    for(row in 1:nrow(table)){
      combinedTable[row + from, ] <- table[row, ]
    }
    from <- from + nrow(table)
  }
  table <- d
  if(nrow(table) > 0){
    for(row in 1:nrow(table)){
      combinedTable[row + from, ] <- table[row, ]
    }
    from <- from + nrow(table)
  }
  table <- e
  if(nrow(table) > 0){
    for(row in 1:nrow(table)){
      combinedTable[row + from, ] <- table[row, ]
    }
    from <- from + nrow(table)
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
