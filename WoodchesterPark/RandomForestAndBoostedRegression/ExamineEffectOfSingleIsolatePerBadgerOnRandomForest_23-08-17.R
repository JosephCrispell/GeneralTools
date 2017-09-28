##################
# Load libraries #
##################

library(randomForest)
library(gplots)

#######################################################
# Open the Genetic Vs Epidemiological distances table #
#######################################################

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_13-07-17/"

# Read in the table
file <- paste(path, "GeneticVsEpidemiologicalDistances/", 
              "GeneticVsEpidemiologicalDistances_10-08-17.txt", sep="")
geneticVsEpi <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

####################
# General settings #
####################

trainProp <- 0.5
colToUse <- "%IncMSE"

# Note the full names of metrics and assign them a colour
fullNames <- noteFullNames()
nameColours <- assignMetricColours(temporalCol="darkgoldenrod4", spatialCol="red",
                                   networkCol="blue")

###############
# Select data #
###############

# Select Badger-Badger comparisons
geneticVsEpi <- selectAppropriateComparisonsForSelection("BB")

# Only select small genetic distances
geneticVsEpi <- selectGeneticDistancesBelowThreshold(threshold=15)

####################################
# Select single isolate per badger #
####################################

# Read in badger sampling information - store the animal associated with each isolate
file <- paste(path, "/IsolateData/BadgerInfo_08-04-15_LatLongs_XY_Centroids.csv",
              sep="")
isolateTattoos <- noteIsolateTattoos(file)

# Note the isolates available for each badger
badgerIsolates <- noteIsolatesAvailableForEachBadger(geneticVsEpi, isolateTattoos)

# Note the isolate sequence quality
file <- paste(path, "vcfFiles/", "sequences_Prox-10_01-08-2017.fasta", sep="")
isolatePropNs <- calculatePropotionNsForEachIsolate(file)

# Select single isolate for each badger - based upon prop Ns
selectedIsolates <- selectSingleIsolatePerBadgerBasedOnQuality(badgerIsolates,
                                                               isolatePropNs)

# Only keep comparisons for selected badgers
geneticVsEpi <- geneticVsEpi[geneticVsEpi$IsolateI %in% selectedIsolates, ]
geneticVsEpi <- geneticVsEpi[geneticVsEpi$IsolateJ %in% selectedIsolates, ]

##########################
# Remove unecessary data #
##########################

# Remove irrelevant columns
geneticVsEpi <- removeColumnsIfNotRelevant(geneticVsEpi)

# Convert the columns dealing with boolean metrics to factors
geneticVsEpi <- makeBooleanColumnsFactors(table=geneticVsEpi)

# Note columns to ignore as predictors
colNamesToIgnore <- c("iSpeciesJSpecies", "IsolateI", "IsolateJ")
colsToIgnore <- which(names(geneticVsEpi) %in% colNamesToIgnore)

###########################
# Fit Random Forest model #
###########################

# Build a test and training data set for predictions
trainRows <- sample(x=1:nrow(geneticVsEpi),
                    size=floor(trainProp * nrow(geneticVsEpi)), replace=FALSE)

# Find optimal mtry parameter
optimalMtry <- findOptimalMtry(response=geneticVsEpi[trainRows, "GeneticDistance"],
                               predictors=geneticVsEpi[trainRows, -c(1, colsToIgnore)],
                               mTryInitial=3, nTrees=500, plot=TRUE)


# Train the Random Forest model
infoRF <- randomForest(geneticVsEpi[trainRows, "GeneticDistance"] ~ ., 
                       data=geneticVsEpi[trainRows, -c(1, colsToIgnore)],
                       mtry=optimalMtry, importance=TRUE, ntree=1000,
                       keep.forest=TRUE, norm.votes=FALSE, proximity=FALSE,
                       do.trace=FALSE)

# Get the Pseudo RSquared value
rSq <- round(infoRF$rsq[length(infoRF$rsq)], digits=2)

# Examine trained model prediction
predictedGeneticDistances <- predict(infoRF, geneticVsEpi[-trainRows, -c(1, colsToIgnore)])
corr <- cor(geneticVsEpi[-trainRows, "GeneticDistance"], predictedGeneticDistances)
plotPredictedVersusActual(actual=geneticVsEpi[-trainRows, "GeneticDistance"],
                          predicted=predictedGeneticDistances, rSq=rSq)

#############
# FUNCTIONS #
#############

plotPredictedVersusActual <- function(actual, predicted, rSq){
  plot(x=predicted, y=actual,
       las=1, pch=19, col=rgb(0,0,0, 0.1), xlab="Predicted", ylab="Actual")
  abline(lm(actual ~ predicted), col="red")
  correlation <- round(cor(actual, predicted), digits=2)
  
  legend("topleft", c(paste("corr =", correlation),
                      paste("Rsq = ", rSq)),
         bty="n", cex = 1)
}

selectSingleIsolatePerBadgerBasedOnQuality <- function(badgerIsolates, 
                                                       isolatePropNs){
  selectedIsolates <- c()
  
  for(tattoo in names(badgerIsolates)){
    
    # Get the isolate IDs for the current badger
    isolates <- badgerIsolates[[tattoo]]
    
    # If more than one isolate available - select one
    if(length(isolates) > 1){
      
      # Get the proportion of Ns for each isolate
      propNs <- getValuesFromList(isolates, isolatePropNs)
      
      # Note the order
      order <- order(propNs, decreasing=TRUE)
      
      # Select the isolate with highest
      selectedIsolates[length(selectedIsolates) + 1] <- isolates[order[1]]
      
    }else{
      selectedIsolates[length(selectedIsolates) + 1] <- isolates[1]
    }
  }
  
  return(selectedIsolates)
}

getValuesFromList <- function(keys, list){
  
  values <- c()
  for(i in 1:length(keys)){
    values[i] <- list[[keys[i]]]
  }
  
  return(values)
}

calculatePropotionNsForEachIsolate <- function(fastaFileName){
  
  con <- file(fastaFileName,open="r")
  fileLines <- readLines(con)
  close(con)
  
  isolatePropNs <- list()
  for(i in 2:length(fileLines)){
    
    if(grepl(x=fileLines[i], pattern=">") == TRUE){
      
      # Calculate proportion Ns for previous isolate's sequence
      if(i != 2 && grepl(x=isolate, pattern="WB") == TRUE){
        
        isolatePropNs[[isolate]] <- calculateProportionNs(sequence)
      }
      
      # Get isolate ID
      isolate <- strsplit(fileLines[i], split="_")[[1]][1]
      isolate <- substr(isolate, start=2, stop=nchar(isolate))
      
      # Reset sequence
      sequence <- ""
    }else{
      sequence <- paste(sequence, fileLines[i], sep="")
    }
  }
  
  return(isolatePropNs)
}

calculateProportionNs <- function(sequence){
  
  propNs <- 0
  nucleotides <- strsplit(sequence, split="")[[1]]
  for(nucleotide in nucleotides){
    if(nucleotide == 'N'){
      propNs <- propNs + 1
    }
  }
  propNs <- propNs / length(nucleotides)
  
  return(propNs)
}

noteIsolatesAvailableForEachBadger <- function(geneticVsEpi, isolateTattoos){
  
  isolates <- unique(c(geneticVsEpi$IsolateI, geneticVsEpi$IsolateJ))
  badgerIsolates <- list()
  for(isolate in isolates){
    
    # Check if badger associated with current isolates exists in list
    if(is.null(badgerIsolates[[isolateTattoos[[isolate]]]]) == FALSE){
      badgerIsolates[[isolateTattoos[[isolate]]]] <- c(
        badgerIsolates[[isolateTattoos[[isolate]]]], isolate)
    }else{
      badgerIsolates[[isolateTattoos[[isolate]]]] <- c(isolate)
    }
  }
  
  return(badgerIsolates)
}

noteIsolateTattoos <- function(samplingInfoFile){
  
  badgerInfo <- read.table(samplingInfoFile, header=TRUE, stringsAsFactors=FALSE,
                           sep=",")
  
  isolateTattoos <- list()
  for(row in 1:nrow(badgerInfo)){
    isolateTattoos[[badgerInfo[row, "WB_id"]]] <- badgerInfo[row, "tattoo"]
  }
  
  return(isolateTattoos)
}

selectGeneticDistancesBelowThreshold <- function(threshold){
  hist(geneticVsEpi$GeneticDistance, breaks=100,
       las=1,
       xlab="Genetic Distance (SNPs)",
       main="Inter-Isolate Genetic Distance Distribution")
  abline(v=threshold, col="red", lty=2)
  
  geneticVsEpi <- geneticVsEpi[geneticVsEpi$GeneticDistance < threshold, ]
  
  return(geneticVsEpi)
}

selectAppropriateComparisonsForSelection <- function(selection){
  if(selection != "CB" && selection != "BC"){
    geneticVsEpi <- geneticVsEpi[geneticVsEpi$iSpeciesJSpecies == selection, ]
  }else{
    geneticVsEpi <- geneticVsEpi[geneticVsEpi$iSpeciesJSpecies != "BB" &
                                   geneticVsEpi$iSpeciesJSpecies != "CC", ]
  }
  
  return(geneticVsEpi)
}

assignMetricColours <- function(temporalCol, spatialCol, networkCol){
  
  nameColours <- list(
    "SameMainGroup" = spatialCol,
    "SameSampledGroup" = spatialCol,                           
    "SameInfectedGroup" = spatialCol,
    "PeriodSpentAliveTogether" = temporalCol,                   
    "PeriodSpentInfectedTogether" = temporalCol,
    "PeriodSpentInSameGroup" = temporalCol,                     
    "TimeBetweenInfectionDetection" = temporalCol,
    "TimeBetweenSampling" = temporalCol,                        
    "TimeBetweenBreakdown" = temporalCol,
    "DistanceBetweenMainGroups" = spatialCol,
    "DistanceBetweenSampledGroups" = spatialCol,               
    "DistanceBetweenInfectedGroups" = spatialCol,
    "NMovementsBetweenMainGroups" = networkCol,                
    "NMovementsBetweenSampledGroups" = networkCol,
    "NMovementsBetweenInfectedGroups" = networkCol,            
    "SameAnimal" = "black",
    "ShortestPathLengthMain" = networkCol,                     
    "MeanNMovementsOnEdgesOfShortestPathMain" = networkCol,
    "ShortestPathLengthSampled" = networkCol,                  
    "MeanNMovementsOnEdgesOfShortestPathSampled" = networkCol,
    "ShortestPathLengthInfected" = networkCol,                 
    "MeanNMovementsOnEdgesOfShortestPathInfected" = networkCol,
    "NSharedAnimalsBetweenMainGroups" = networkCol,            
    "NSharedAnimalsBetweenSampledGroups" = networkCol,
    "NSharedAnimalsBetweenInfectedGroups" = networkCol,
    "ShortestPathLengthEXCLMain" = networkCol,                     
    "MeanNMovementsOnEdgesOfShortestPathEXCLMain" = networkCol,   
    "ShortestPathLengthEXCLSampled" = networkCol,               
    "MeanNMovementsOnEdgesOfShortestPathEXCLSampled" = networkCol,
    "ShortestPathLengthEXCLInfected" = networkCol,
    "MeanNMovementsOnEdgesOfShortestPathEXCLInfected" = networkCol,
    "CentroidDistBetweenMain" = spatialCol,
    "CentroidDistBetweenSamp" = spatialCol,
    "HostRelatedness" = "black"
  )
  
  return(nameColours)
}

noteFullNames <- function(){
  fullNames <- list(
    "SameMainGroup" = "Isolates taken from same main group (yes/no)",
    "SameSampledGroup" = "Isolates taken from same sampled group (yes/no)",                           
    "SameInfectedGroup" = "Isolates taken from same infected group (yes/no)",
    "PeriodSpentAliveTogether" = "Number of days overlap between recorded lifespans of the sampled animals",                   
    "PeriodSpentInfectedTogether" = "Number of days overlap between infected lifespans of the sampled animals",
    "PeriodSpentInSameGroup" = "Number of days that sampled animals spent in same group",                     
    "TimeBetweenInfectionDetection" = "Number of days between infection detection dates",
    "TimeBetweenSampling" = "Number of days between sampling dates",                        
    "TimeBetweenBreakdown" = "Number of days between breakdown dates",
    "DistanceBetweenMainGroups" = "Spatial distance (km) between main groups",
    "DistanceBetweenSampledGroups" = "Spatial distance (km) between sampled groups",               
    "DistanceBetweenInfectedGroups" = "Spatial distance (km) between infected groups",
    "NMovementsBetweenMainGroups" = "Number of recorded animal movements between main groups of sampled animals",                
    "NMovementsBetweenSampledGroups" = "Number of recorded animal movements between sampled groups of sampled animals",
    "NMovementsBetweenInfectedGroups" = "Number of recorded animal movements between infected groups of sampled animals",            
    "SameAnimal" = "Isolates taken from same badger (yes/no)",
    "ShortestPathLengthMain" = "Shortest path length between main groups of sampled animals",                     
    "MeanNMovementsOnEdgesOfShortestPathMain" = "Mean number of animals dispersing along edges of shortest path between main groups",
    "ShortestPathLengthSampled" = "Shortest path length between sampled groups of sampled animals",                  
    "MeanNMovementsOnEdgesOfShortestPathSampled" = "Mean number of animals dispersing along edges of shortest path between sampled groups",
    "ShortestPathLengthInfected" = "Shortest path length between infected groups of sampled animals",                 
    "MeanNMovementsOnEdgesOfShortestPathInfected" = "Mean number of animals dispersing along edges of shortest path between infected groups",
    "NSharedAnimalsBetweenMainGroups" = "Number of animals recorded in both main groups of sampled animals",            
    "NSharedAnimalsBetweenSampledGroups" = "Number of animals recorded in both sampled groups of sampled animals",
    "NSharedAnimalsBetweenInfectedGroups" = "Number of animals recorded in both infected groups of sampled animals",
    "ShortestPathLengthEXCLMain" = "Shortest path length between main herds of sampled animals (Some Herds Excluded)",                     
    "MeanNMovementsOnEdgesOfShortestPathEXCLMain" = "Mean number of animals dispersing along edges of shortest path between main herds (Some Herds Excluded)",   
    "ShortestPathLengthEXCLSampled" = "Shortest path length between sampled herds of sampled animals (Some Herds Excluded)",               
    "MeanNMovementsOnEdgesOfShortestPathEXCLSampled" = "Mean number of animals dispersing along edges of shortest path between sampled herds (Some Herds Excluded)",
    "ShortestPathLengthEXCLInfected" = "Shortest path length between infected herds of sampled animals (Some Herds Excluded)",
    "MeanNMovementsOnEdgesOfShortestPathEXCLInfected" = "Mean number of animals dispersing along edges of shortest path between main herds (Some Herds Excluded)",
    "CentroidDistBetweenMain" = "Distance from centroid of closest land parcel to badgers main sett",
    "CentroidDistBetweenSamp" = "Distance from centroid of closest land parcel to badgers sampled sett",
    "HostRelatedness" = "Genetic relatedness of sampled badgers"
  )
  return(fullNames)
}

getVariableColours <- function(variableImportance, nameColours){
  
  rowNames <- rownames(variableImportance)
  
  colours <- c()
  for(index in 1:nrow(variableImportance)){
    colours[index] <- nameColours[[rowNames[index]]]
  }
  
  return(colours)
}

returnVariableColour <- function(name, spatial, temporal, network, default){
  
  spatialPatterns <- c("DistanceBetween", "Same[a-zA-Z]+Group")
  temporalPatterns <- c("PeriodSpent", "TimeBetween")
  networkPatterns <- c("NMovementsBetween", "ShortestPath", "MeanNMovements", "NSharedAnimals")
  
  colour <- default
  if(checkForMatch(name, spatialPatterns) == TRUE){
    colour <- spatial
  }else if(checkForMatch(name, temporalPatterns) == TRUE){
    colour <- temporal
  }else if(checkForMatch(name, networkPatterns) == TRUE){
    colour <- network
  }
  
  return(colour)
}

checkForMatch <- function(x, patterns){
  match <- FALSE
  for(pattern in patterns){
    match <- grepl(x=x, pattern=pattern)
    if(match == TRUE){
      break
    }
  }
  return(match)
}

getFullVariableNames <- function(array, fullNames){
  
  output <- c()
  for(i in 1:length(array)){
    
    output[i] <- fullNames[[array[i]]]
    
    if(is.null(fullNames[[array[i]]]) == TRUE){
      print(array[i])
    }
  }
  
  return(output)
}

findOptimalMtry <- function(response, predictors, mTryInitial, nTrees, plot){
  
  tuneOutput <- tuneRF(predictors, response, mtryStart=mTryInitial,
                       ntreeTry=nTrees, stepFactor=1.5, improve=0.0001, trace=FALSE,
                       plot=FALSE)
  
  optimalMtry <- as.integer(rownames(tuneOutput)[tuneOutput[,2] == min(tuneOutput[,2])])
  
  if(plot == TRUE){
    plot(tuneOutput, las=1)
    abline(v=optimalMtry, col="red", lty=2)
  }
  
  return(optimalMtry)
}

makeBooleanColumnsFactors <- function(table){
  
  colNames <- colnames(table)
  
  for(col in 1:ncol(table)){
    
    if(grepl(x=colNames[col], pattern="Same") == TRUE &
       grepl(x=colNames[col], pattern="PeriodSpentIn") == FALSE ){
      table[, col] <- as.factor(table[, col])
    }
  }
  
  return(table)
}

removeColumnsIfNotRelevant <- function(table){
  
  # For particular comparisons: Badger-Badger, Badger-Cattle, Cattle-Cattle
  # some epidemiological metrics aren't relevant and column will be filled 
  # with -1
  
  colsToRemove <- c()
  index <- 0
  
  for(col in 3:(ncol(table)-2)){
    
    if(sd(table[, col]) == 0){
      index <- index + 1
      colsToRemove[index] <- col
      cat(paste("Removed: ", colnames(table)[col], "\n", sep=""))
    }
  }
  return(table[, -colsToRemove])
}
