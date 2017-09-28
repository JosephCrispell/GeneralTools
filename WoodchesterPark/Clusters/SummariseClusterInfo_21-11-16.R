###################################
# Load and process the input data #
###################################

# Create a path variable
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/"

## Load the life history summaries
file <- paste(path, "InterSpeciesClusters/sampledAnimalsLifeHistories_17-11-2016.txt", sep="")
table <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep="\t",
                    colClasses = "character")

# Remove test types not interested in
testTypesToIgnore <- list("CT"=1, "DC"=1)
table <- removeTestDatesForCertainTestTypes(table, testTypesToIgnore)


## Add locations of breakdowns

# Cattle sampling info
file <- paste(path, "IsolateData/",
              "CattleIsolateInfo_LatLongs_plusID_outbreakSize_Coverage_AddedTB1453-TB1456.csv",
              sep="")
cattleSamplingInfo <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep=",",
                                 colClasses = "character")
# CTS Location data
file <- paste(path, "CattleMovementData/20160314_joe_cts_locations.csv", sep="")
locationInfo <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep=",", fill=TRUE)

# Add the breakdown coordinates
mansionX <- 380909
mansionY <- 201377
threshold <- 15000
table <- addBreakdownCoordinates(table, cattleSamplingInfo, locationInfo, mansionX, 
                                 mansionY, threshold)

# Create a table for each cluster
clusterTables <- splitInputTableIntoInfoForClusters(table)

#########################################
# Create Summary Table For Each Cluster #
#########################################

# Note the cattle premises types to ignore
premisesToIgnore <- list("SR"=1, "SW"=1, "EX"=1)

# Get an array of all the clusters
clusters <- sort(names(clusterTables))

# Initialise a table to store the cluster summaries
summary <- data.frame(
  # Cluster
  "Cluster" = clusters,
  
  # Counts
  "NumberSampledBadgers" = rep(0, length(clusters)),
  "NumberSampledCattle" = rep(0, length(clusters)),
  "NumberUnSampledDetectedBadgers" = rep(0, length(clusters)),
  "NumberUnSampledCattleReactors" = rep(0, length(clusters)),
  "NumberNegativeBadgers" = rep(0, length(clusters)),
  "NumberNegativeCattle" = rep(0, length(clusters)),
  
  # Earliest Detection
  "EarliestDetectionDateForSampledBadgers" = rep(as.Date("2016-11-21"), length(clusters)),
  "EarliestDetectionDateForSampledCattle" = rep(as.Date("2016-11-21"), length(clusters)),
  "EarliestDetectionDateForUnSampledBadgers" = rep(as.Date("2016-11-21"), length(clusters)),
  "EarliestDetectionDateForUnSampledCattle" = rep(as.Date("2016-11-21"), length(clusters)),
  
  # Min MRCA Distance
  "MinDistanceOfBadgersToMRCA" = rep(0, length(clusters)),
  "MinDistanceOfCattleToMRCA" = rep(0, length(clusters)),
  
  # Mean Spatial Distance
  "MeanSpatialDistanceOfSampledHerds" = rep(0, length(clusters)),
  
  # Degree of Herds Sampled Cattle Lived on
  "MeanNumberSampledCattleMovementsToFromSampledHerds" = rep(0, length(clusters)),
  "MeanNumberReactorCattleMovementsToFromSampledHerds" = rep(0, length(clusters)),
  "MeanNumberNegativeCattleMovementsToFromSampledHerds" = rep(0, length(clusters))
)

for(index in 1:length(clusters)){
  cluster <- clusters[index]

  clusterTable <- clusterTables[[cluster]]
  
  rowValues <- summariseCluster(clusterTable, premisesToIgnore, mansionX, mansionY)
  
  summary[index, 2:ncol(summary)] <- rowValues
}

# Convert the distances to MRCA to SNPs
fastaLength <- 9464
summary[, 12] <- round(summary[, 12] * fastaLength, digits=2)
summary[, 13] <- round(summary[, 13] * fastaLength, digits=2)

# Round the columns of the summary table
summary[, 14] <- round(summary[, 14], digits=2)
summary[, 15] <- round(summary[, 15], digits=2)
summary[, 16] <- round(summary[, 16], digits=2)
summary[, 17] <- round(summary[, 17], digits=2)

# Flip the table
output <- as.data.frame(t(summary))
colnames(output) <- paste("Cluster-", clusters, sep="")
output <- output[-1,]

# Define new rownames
rownames(output) <- c(
  "Number of badgers sampled",
  "Number of cattle sampled",
  "Number of in-contact badgers that tested positive",
  "Number of in-contact cattle that tested positive",
  "Number of in-contact badgers that NEVER tested positive",
  "Number of in-contact cattle that NEVER tested positive",
  "Earliest date that a sampled badger tested positive",
  "Earliest date that a sampled cow tested positive",
  "Earliest date that an in-contact badger tested positive",
  "Earliest date that an in-contact cow tested positive",
  "Minimum patristic distance (SNPs) of the sampled badgers to the MRCA of cluster",
  "Minimum patristic distance (SNPs) of the sampled cattle to the MRCA of cluster",
  "Mean spatial distance (KM) from the sampled herds to Woodchester Park",
  "Mean number of movements of sampled cattle to or from the sampled herds",
  "Mean number of movements of in-contact animals that tested positive to or from the sampled herds",
  "Mean number of movements of in-contact animals that NEVER tested positive to or from the sampled herds"
)

# Print the table out to file
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/InterSpeciesClusters/"
file <- paste(path, "clustersSummaryTable_29-11-16.csv", sep="")
write.table(output, file=file, row.names=TRUE, quote=FALSE, sep=",")


#############
# FUNCTIONS #
#############

###########
# Summary #
###########

addBreakdownCoordinates <- function(table, cattleSamplingInfo, locationInfo, mansionX, 
                                    mansionY, threshold){
  
  cattleSamplingList <- convertTableToList("StrainId", cattleSamplingInfo)
  
  locationInfo <- keepLocationsWithinThresholdDistance(locationInfo=locationInfo, 
                                                       thresholdInMetres=threshold,
                                                       mansionX=mansionX,
                                                       mansionY=mansionY)
  
  # Get the CPH into the right format
  locationInfo$CPH <- formatCPHForLocationInfo(locationInfo$cph)
  locationList <- convertTableToList("CPH", locationInfo)
  
  
  # Get the coordinates of each breakdown location for each sampled cow
  table$BreakdownX <- rep(NA, nrow(table))
  table$BreakdownY <- rep(NA, nrow(table))
  
  for(row in 1:nrow(table)){
    
    # Skip badgers and unsampled
    if(is.na(table[row, "Isolates"]) == TRUE || table[row, "Species"] == "BADGER"){
      next
    }
    
    # Get the Isolate ID
    isolateId <- table[row, "Isolates"]
    
    # Get the CPH from the breakdown ID
    samplingInfo <- cattleSamplingList[[isolateId]]
    breakdownCPH <- strsplit(samplingInfo[1, "BreakdownID"], split="-")[[1]][1]
    breakdownCPH <- substr(breakdownCPH, start=1, stop=nchar(breakdownCPH) - 2)
    
    # Get the X and Y coordinates of location
    breakdownLocationInfo <- locationList[[breakdownCPH]]
    table[row, "BreakdownX"] <- paste(breakdownLocationInfo[1, "x"], ".0", sep="")
    table[row, "BreakdownY"] <- paste(breakdownLocationInfo[1, "y"], ".0", sep="")
  }
  
  return(table)
}

summariseCluster <- function(clusterTable, premisesToIgnore, mansionX, mansionY){
  
  # Initialise an array to keep counts
  counts <- rep(0, 6) # SampledBadgers, SampledCattle, DetectedBadgers, ReactorCattle, UnDetectedBadgers, NegativeCattle
  
  # Initialise a variable to note which species infection was detected in first and when
  earliestDetectionDates <- rep(as.Date("2016-11-21"), 4) # SampledBadgers, SampledCattle, unSampledBadgers, unSampledCattle
  
  # Initialise a variable to note which species was closer to the root
  minMRCADistance <- c(999999, 999999) # Badgers, Cattle
  
  # Initialise variables to record the herds involved
  sampledHerds <- list()
  
  # Examine each animal associated with the cluster
  for(row in 1:nrow(clusterTable)){

    # Record information for badger
    if(clusterTable[row, "Species"] == "BADGER"){
      
      # Was this animal sampled?
      if(is.na(clusterTable[row, "Isolates"]) == FALSE){
        
        # Add Sampled badger to count
        counts[1] <- counts[1] + 1
        
        # Was the infection in sampled badger detected earliest?
        detectionDate <- as.Date(clusterTable[row, "DetectionDate"], "%d-%m-%Y")
        if(detectionDate - earliestDetectionDates[1] < 0){
          earliestDetectionDates[1] <- detectionDate
        }
        
        # Is the distance from this badger's isolate to the MRCA the minimum?
        distancesToMRCA <- as.numeric(strsplit(clusterTable[row, "DistancesToMRCA"], split=",")[[1]])
        for(i in 1:length(distancesToMRCA)){
          if(distancesToMRCA[i] < minMRCADistance[1]){
            minMRCADistance[1] <- distancesToMRCA[i]
          }
        }

      }else{
        
        # Was infection ever detected in this badger?
        if(is.na(clusterTable[row, "DetectionDate"]) == FALSE){
          
          # Add Detected Badger to count
          counts[3] <- counts[3] + 1
          
          # Was the infection in this unsampled badger detected earliest?
          detectionDate <- as.Date(clusterTable[row, "DetectionDate"], "%d-%m-%Y")
          if(detectionDate - earliestDetectionDates[3] < 0){
            earliestDetectionDates[3] <- detectionDate
          }
          
        }else{
          
          # Add Negative Badger to count
          counts[5] <- counts[5] + 1
        }
      }
      
      # Record information for cow  
    }else{
      
      # Was this animal sampled?
      if(is.na(clusterTable[row, "Isolates"]) == FALSE){
        
        # Add sampled cow to count
        counts[2] <- counts[2] + 1
        
        # Was the infection in this sampled cow detected ealiest?
        testResults <- strsplit(clusterTable[row, "CattleTestResults"], split=",")[[1]]
        testDates <- as.Date(strsplit(
          as.character(clusterTable[row, "CattleTestDates"]), split=",")[[1]], "%d-%m-%Y")
        for(i in 1:length(testResults)){
          if(testResults[i] == "R" || testResults[i] == "SL"){
            
            # Note whether this unsampled cow's disease was detected earliest
            if(testDates[i] - earliestDetectionDates[2] < 0){
              earliestDetectionDates[2] <- testDates[i]
            }
          }
        }
        
        # Is the distance from this badger's isolate to the MRCA the minimum?
        distanceToMRCA <- as.numeric(clusterTable[row, "DistancesToMRCA"])
        if(distanceToMRCA < minMRCADistance[2]){
          minMRCADistance[2] <- distanceToMRCA
        }
        
        # Note the herd that this sampled cow was broke down on
        X <- clusterTable[row, "BreakdownX"]
        Y <- clusterTable[row, "BreakdownY"]
        if(is.na(X) == FALSE && X != "NA"){
          key <- paste(X, ":", Y, sep="")
          if(is.null(sampledHerds[[key]]) == TRUE){
            sampledHerds[[key]] <- c(0, 0, 0)
          }
        }
        
        # Examine the herds this sampled cow lived on
        Xs <- strsplit(clusterTable[row, "Xs"], split=",")[[1]]
        Ys <- strsplit(clusterTable[row, "Ys"], split=",")[[1]]
         for(i in 1:length(Xs)){
          # Don't look at NA locations
          if(is.na(Xs[i]) == FALSE && Xs[i] != "NA"){
            
            # Create key
            key <- paste(Xs[i], ":", Ys[i], sep="")
            
            # Note group, if haven't already
            if(is.null(sampledHerds[[key]]) == FALSE){

              # Add the movement of this sampled cow to the herds degree
              sampledHerds[[key]] <- sampledHerds[[key]] + c(1,0,0)
            }
          }
        }
        
      }else{
        
        # Did this cow ever react to the SICCT?
        if(is.na(clusterTable[row, "CattleTestResults"]) == FALSE){
          
          # Was this cow a reactor?
          testResults <- strsplit(clusterTable[row, "CattleTestResults"], split=",")[[1]]
          testDates <- as.Date(strsplit(
            as.character(clusterTable[row, "CattleTestDates"]), split=",")[[1]], "%d-%m-%Y")
          for(i in 1:length(testResults)){
            if(testResults[i] == "R" || testResults[i] == "SL"){
              
              # Add Reactor unsampled cow to count
              counts[4] <- counts[4] + 1
              
              # Note whether this unsampled cow's disease was detected earliest
              if(testDates[i] - earliestDetectionDates[4] < 0){
                earliestDetectionDates[4] <- testDates[i]
              }
              
              # Note each of the herds that this sampled cow lived in
              Xs <- strsplit(clusterTable[row, "Xs"], split=",")[[1]]
              Ys <- strsplit(clusterTable[row, "Ys"], split=",")[[1]]
              premisesTypes <- strsplit(clusterTable[row, "PremisesTypes"], split=",")[[1]]
              for(i in 1:length(Xs)){
                # Don't look at NA locations
                if(is.na(Xs[i]) == FALSE && Xs[i] != "NA" && is.null(premisesToIgnore[[premisesTypes[i]]]) == TRUE){
                  
                  # Create key
                  key <- paste(Xs[i], ":", Ys[i], sep="")
                  
                  # Was this a herd a sampled cow lived in? If so, add to degree
                  if(is.null(sampledHerds[[key]]) == FALSE){
                    # Add the movement of this sampled cow to the herds degree
                    sampledHerds[[key]] <- sampledHerds[[key]] + c(0,1,0)
                  }
                }
              }
            }
          }
          
        }else{
          
          # Add Negative Count to counts
          counts[6] <- counts[6] + 1
          
          # Note each of the herds that this sampled cow lived in
          Xs <- strsplit(clusterTable[row, "Xs"], split=",")[[1]]
          Ys <- strsplit(clusterTable[row, "Ys"], split=",")[[1]]
          premisesTypes <- strsplit(clusterTable[row, "PremisesTypes"], split=",")[[1]]
          for(i in 1:length(Xs)){
            # Don't look at NA locations
            if(is.na(Xs[i]) == FALSE && Xs[i] != "NA" && is.null(premisesToIgnore[[premisesTypes[i]]]) == TRUE){
              
              # Create key
              key <- paste(Xs[i], ":", Ys[i], sep="")
              
              # Was this a herd a sampled cow lived in? If so, add to degree
              if(is.null(sampledHerds[[key]]) == FALSE){
                # Add the movement of this sampled cow to the herds degree
                sampledHerds[[key]] <- sampledHerds[[key]] + c(0,0,1)
              }
            }
          }
        }
      }
    }
  }
  
  ## Examine the sampled herds
  
  # Initialise a variable to record the mean spatial distance to Woodchester Park
  meanSpatialDistance <- 0
  
  # Initialise an array to store the mean degree values for each herd sampled animals lived in
  meanDegree <- c(0,0,0)
  
  keys <- names(sampledHerds)
  for(key in keys){
    coordinates <- as.numeric(strsplit(key, split=":")[[1]])
    
    # Calculate distance to Woodchester Mansion
    meanSpatialDistance <- meanSpatialDistance + 
      euclideanDistance(mansionX, mansionY, coordinates[1], coordinates[2])
    
    # Add to meanDegree
    meanDegree <- meanDegree + sampledHerds[[key]]
  }
  meanSpatialDistance <- (meanSpatialDistance / length(keys)) / 1000
  meanDegree <- meanDegree / length(meanDegree)
  
  # Return the summary information
  summary <- data.frame(# Counts
                        "NumberSampledBadgers" = counts[1],
                        "NumberSampledCattle" = counts[2],
                        "NumberUnSampledDetectedBadgers" = counts[3],
                        "NumberUnSampledCattleReactors" = counts[4],
                        "NumberNegativeBadgers" = counts[5],
                        "NumberNegativeCattle" = counts[6],
                        
                        # Earliest Detection
                        "EarliestDetectionDateForSampledBadgers" = earliestDetectionDates[1],
                        "EarliestDetectionDateForSampledCattle" = earliestDetectionDates[2],
                        "EarliestDetectionDateForUnSampledBadgers" = earliestDetectionDates[3],
                        "EarliestDetectionDateForUnSampledCattle" = earliestDetectionDates[4],
                        
                        # Min MRCA Distance
                        "MinDistanceOfBadgersToMRCA" = minMRCADistance[1],
                        "MinDistanceOfCattleToMRCA" = minMRCADistance[2],
                        
                        # Mean Spatial Distance
                        "MeanSpatialDistanceOfSampledHerds" = meanSpatialDistance,
                        
                        # Degree of Herds Sampled Cattle Lived on
                        "MeanNumberSampledCattleMovementsToFromSampledHerds" = meanDegree[1],
                        "MeanNumberReactorCattleMovementsToFromSampledHerds" = meanDegree[2],
                        "MeanNumberNegativeCattleMovementsToFromSampledHerds" = meanDegree[3]
                        )
  
  return(summary[1, ])
}


###########
# General #
###########

removeTestDatesForCertainTestTypes <- function(table, testTypes){
  
  for(row in 1:nrow(table)){
    
    testDates <- strsplit(table[row, "CattleTestDates"], split=",")[[1]]
    testResults <- strsplit(table[row, "CattleTestResults"], split=",")[[1]]
    
    remove <- c()
    index <- 0
    
    # Do we want to ignore any of the tests?
    for(i in 1:length(testResults)){
      
      if(is.null(testTypes[[testResults[i]]]) == FALSE){
        index <- index + 1
        remove[index] <- i
      }
    }
    
    # Did we find tests to ignore
    if(length(remove) > 0){
      testDates <- testDates[-remove]
      testResults <- testResults[-remove]
      
      if(length(testDates) != 0){
        table[row, "CattleTestDates"] <- paste(testDates, collapse=",")
        table[row, "CattleTestResults"] <- paste(testResults, collapse=",")
      }else{
        table[row, "CattleTestDates"] <- NA
        table[row, "CattleTestResults"] <- NA
      }
    }
  }
  
  return(table)
}


keepLocationsWithinThresholdDistance <- function(locationInfo, thresholdInMetres, 
                                                 mansionX, mansionY){
  
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
    distance <- euclideanDistance(x1=mansionX, y1=mansionY, 
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


formatCPHForLocationInfo <- function(cphs){
  
  formattedCPHs <- c()
  for(i in 1:length(cphs)){
    
    formattedCPHs[i] <- paste(strsplit(
      strsplit(cphs[i], split="-")[[1]][1], split="/")[[1]], collapse="")
  }
  
  return(formattedCPHs)
}

convertTableToList <- function(column, table){
  
  tableAsList <- list()
  for(row in 1:nrow(table)){
    
    # Skip NAs
    if(is.na(table[row, column]) == TRUE){
      next
    }
    
    tableAsList[[table[row, column]]] <- table[row, ]
  }
  
  return(tableAsList)
}


euclideanDistance <- function(x1, y1, x2, y2){
  return(sqrt(sum((x1 - x2)^2 + (y1 - y2)^2)))
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