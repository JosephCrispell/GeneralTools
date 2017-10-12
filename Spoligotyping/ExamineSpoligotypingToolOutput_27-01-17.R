#########################################
# Read in the Spoligotyping Tool Output #
#########################################

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_13-07-17/"

# Open the file
file <- paste(path, "Mislabelling/Spoligotyping/",
              "SpoligotypeMatches_28-09-2017.txt", sep="")
matchTable <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
nColToSkip <- 4

# Break the input table up into readable parts
spoligotypeIndices <- indexArray(colnames(matchTable)[(nColToSkip+1):ncol(matchTable)])
isolates <- parseFileNames(matchTable$File)
metaDataSpoligotypes <- matchTable$AssignedSpoligotype
averageReadDepthForIsolates <- matchTable$AverageDepth
proportionNsForIsolates <- matchTable$ProportionNs
isolateMatchingInfo <- getTheMatchingInformationForEachIsolate(matchTable, 
                                                               isolates, nColToSkip)

# Note other names available for some spoligotypes
ahvlaCodes <- list("SB0140"=9, "SB0272"=10, "SB0263"=17, "SB0275"=15)

####################
# Get Sequence IDs #
####################

# Cattle Isolates
file <- paste(path, "IsolateData/", "LinkTable_StrainID-SequenceNo_20-06-16.csv", sep="")
linkTable <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)
cattleIsolateSequenceIds <- noteSequenceIDsOfCattleIsolates(linkTable)

#####################################
# Get the Isolates' Genome Coverage #
#####################################

file <- paste(path, "vcfFiles/",
              "IsolateVariantPositionCoverage_RESCUED_27-09-2017.txt", sep="")
coverage <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
isolateCoverage <- getIsolateCoverage(coverage)

##########################################################
# Examine the Isolate's Spoligotype Matching Information #
##########################################################

# Create a summary table
summary <- summariseIsolateSpoligotypeMatchingInfo(isolates,
                                                   metaDataSpoligotypes,
                                                   spoligotypeIndices,
                                                   isolateMatchingInfo,
                                                   averageReadDepthForIsolates,
                                                   proportionNsForIsolates)



# Add the sequenceIDs of the cattle isolates
summary <- addCattleIsolateSequenceIDs(summary, cattleIsolateSequenceIds)

# Add AHVLA Spoligotype codes - where available
summary <- addAHVLACodes(summary, ahvlaCodes)

# Add the genome coverage
summary <- addIsolateCoverage(summary, isolateCoverage)
summary$Coverage <- round(summary$Coverage, digits=2)

# Round the Average Read Depth and ProportionNs columns
summary$AverageDepth <- round(summary$AverageDepth, digits=2)
summary$ProportionNs <- round(summary$ProportionNs, digits=2)

# Get a subset of isolates where assigned and best match spoligotypes don't match
misMatches <- summary[is.na(summary$AHVLAType) == FALSE, ]
misMatches <- misMatches[misMatches$AHVLAType != misMatches$BestMatch, ]

# Re-order the table
misMatches <- misMatches[order(misMatches$AHVLAType), ]
misMatches$Reason <- rep("Unknown", nrow(misMatches))
rownames(misMatches) <- seq(1, nrow(misMatches), 1)
misMatches <- addReasonGuesses(misMatches)
# Add some guesses
#misMatches[misMatches$ProportionNs > 0.5, "Reason"] <- "Low Sequencing Quality"
#misMatches[which(grepl(x=misMatches$Isolate, 
#                       pattern="WB174|TB1775|TB1391|TB1398|TB1783")),
#           "Reason"] <- "Close Enough?"
#misMatches[which(grepl(x=misMatches$Isolate, 
#                       pattern="TB1471|TB1491")),
#           "Reason"] <- "Not close to Anything"
#misMatches <- misMatches[order(misMatches$Reason, decreasing=TRUE), ]

#########################################################################
# Add Mean genetic distance to isolates that matched that are same type #
#########################################################################

# Read in the FASTA file
file <- paste(path, "vcfFiles/sequences_Prox-10_27-09-2017.fasta", sep="")
sequences <- readFasta(file)

# Build genetic distance matrix
geneticDistances <- buildGeneticDistanceMatrix(sequences)

# Note index of each isolate from match table in sequences
isolateSequenceIndices <- getIndicesOfIsolatesInArray(isolates, names(sequences))

# Get info for isolates that matched spoligotypes
matched <- summary[summary$AHVLAType == summary$BestMatch, ]
matchedTypes <- getIsolateTypes(matched$BestMatch)

# Calculate mean distance of each mismatched isolate between it and 
# - Isolates carrying its sample info type
# - Isolates carrying its best match type
misMatches$AssignedRatio <- rep(NA, nrow(misMatches))
misMatches$BestRatio <- rep(NA, nrow(misMatches))
for(row in 1:nrow(misMatches)){
  
  misMatches[row, "AssignedRatio"] <- calculateMeanGeneticDistanceToIsolatesOfType(
    isolate=misMatches[row, "Isolate"],
    type=strsplit(misMatches[row, "AHVLAType"], split=" ")[[1]][1],
    geneticDistances=geneticDistances,
    matchedIsolates=matched$Isolate,
    matchedIsolateTypes=matchedTypes,
    isolatesIndicesInMatrix=isolateSequenceIndices
  )
  
  misMatches[row, "BestRatio"] <- calculateMeanGeneticDistanceToIsolatesOfType(
    isolate=misMatches[row, "Isolate"],
    type=strsplit(misMatches[row, "BestMatch"], split=" ")[[1]][1],
    geneticDistances=geneticDistances,
    matchedIsolates=matched$Isolate,
    matchedIsolateTypes=matchedTypes,
    isolatesIndicesInMatrix=isolateSequenceIndices
  )
}


#######################
# Write out the table #
#######################

# Re-order the columns
misMatches <- misMatches[, c(1,2,3,4,5,6,7,8,9,11,12,10)]

file <- paste(path, "Mislabelling/Spoligotyping/",
              "SpoligotypeMatches_Summary_28-09-2017.txt", sep="")
write.table(x=misMatches, file, quote=FALSE, sep="\t", row.names=FALSE)

#############
# Functions #
#############

addReasonGuesses <- function(misMatches){
  
  for(row in 1:nrow(misMatches)){
    
    # Get the Type matching information
    parts <- as.numeric(strsplit(misMatches[row, "AHVLATypeInfo"], split=":")[[1]])
    nMissingInfo <- parts[1]
    nMisMatchInfo <- parts[2]
    parts <- as.numeric(strsplit(misMatches[row, "BestMatchInfo"], split=":")[[1]])
    nMissingBest <- parts[1]
    nMisMatchBest <- parts[2]
    
    # Close enough - no more than 3 missing
    if(nMisMatchInfo == 0 && nMissingInfo <= 3){
      misMatches[row, "Reason"] <- "Close Enough?"
      next
    }
    
    # Note close to anything - more than 3 missing against both
    if(nMisMatchInfo == 0 && nMisMatchBest == 0 &&
       nMissingInfo > 3 && nMissingBest > 3){
      misMatches[row, "Reason"] <- "Not Close to anything"
      next
    }
    
    # Mislabelled - Higher ratio against metadata type than best match
    # Ratio of the mean within / mean between type distance
    parts <- as.numeric(strsplit(misMatches[row, "AssignedRatio"], split=":")[[1]])
    AssignedRatio <- parts[1] / parts[2]
    parts <- as.numeric(strsplit(misMatches[row, "BestRatio"], split=":")[[1]])
    BestRatio <- parts[1] / parts[2]
    if(is.na(BestRatio) == FALSE && is.na(AssignedRatio) == FALSE &&
       AssignedRatio > BestRatio){
      misMatches[row, "Reason"] <- "Mislabelled - better distances"
      next
    }
    
    # If mismatches present against metadata and not for Best match
    if(nMisMatchInfo > 0 && nMisMatchBest == 0){
      misMatches[row, "Reason"] <- "Mislabelled - mismatches"
      next
    }
    
    # If poor region coverage
    if(misMatches[row, "ProportionNs"] > 0.2){
      misMatches[row, "Reason"] <- "Poor sequencing coverage"
      next
    }
  }
  
  return(misMatches)
}

getIsolateTypes <- function(typeInfo){
  
  output <- c()
  for(i in 1:length(typeInfo)){
    
    output[i] <- strsplit(typeInfo[i], split=" ")[[1]][1]
  }
  
  return(output)
}

calculateMeanGeneticDistanceToIsolatesOfType <- function(isolate, type, geneticDistances, matchedIsolates, 
                                                         matchedIsolateTypes, isolatesIndicesInMatrix){
  
  # Initialise a variable to calculate the mean genetic distance when comparing
  # to isolates of the same type and of different
  meanWithin <- NA
  meanBetween <- NA
  withinCounts <- 0
  betweenCounts <- 0
  
  # Calculate mean genetic distance of distances only to isolates of same type
  for(i in 1:length(matchedIsolates)){
    
    # Check if same or different type and add to appropriate running total
    if(type == matchedIsolateTypes[i]){
      if(is.na(meanWithin) == TRUE){
        meanWithin <- 0
      }
      
      meanWithin <- meanWithin + geneticDistances[isolatesIndicesInMatrix[[isolate]], 
                                                  isolatesIndicesInMatrix[[matchedIsolates[i]]]]
      withinCounts <- withinCounts + 1
    }else{
      if(is.na(meanBetween) == TRUE){
        meanBetween <- 0
      }
      meanBetween <- meanBetween + geneticDistances[isolatesIndicesInMatrix[[isolate]], 
                                                    isolatesIndicesInMatrix[[matchedIsolates[i]]]]
      betweenCounts <- betweenCounts + 1
    }
  }
  
  #print("----------------------------------------------------")
  #print(paste("Isolate", isolate, "with type", type))
  #print(paste("Found", withinCounts, "within type distances"))
  #print(paste("Found", betweenCounts, "between type distances"))
  
  # Complete the mean calculation
  meanWithin <- meanWithin / withinCounts
  meanBetween <- meanBetween / betweenCounts
  
  # Return the ratio of the two
  return(paste(round(meanWithin, digits=1), ":", round(meanBetween, digits=1), sep=""))
}

getIndicesOfIsolatesInArray <- function(isolates, array){
  indices <- list()
  
  for(i in 1:length(isolates)){
    
    indices[isolates[i]] <- which(array == isolates[i])
  }
  
  return(indices)
}

buildGeneticDistanceMatrix <- function(sequences){
  
  matrix <- matrix(nrow=length(sequences), ncol=length(sequences))
  
  
  keys <- names(sequences)
  rownames(matrix) <- keys
  colnames(matrix) <- keys
  
  for(i in 1:length(sequences)){
    
    for(j in 1:length(sequences)){
      
      if(i >= j){
        next
      }
      
      distance <- geneticDistance(sequences[[keys[i]]], sequences[[keys[j]]])
      matrix[i, j] <- distance
      matrix[j, i] <- distance
    }
  }
  
  return(matrix)
}

geneticDistance <- function(a, b){
  
  distance <- 0
  
  for(i in 1:length(a)){
    
    if(a[i] != "N" && b[i] != "N" && a[i] != b[i]){
      
      distance <- distance + 1
    }
  }
  
  return(distance)
}

readFasta <- function(fileName){
  
  # Store all file lines
  connection <- file(fileName, open="r")
  fileLines <- readLines(connection)
  close(connection)
  
  # Initialise a list to store the fasta sequences
  sequences <- list()
  
  # Examine each line - skip first line
  for(i in 2:length(fileLines)){
    
    # Check if sequence header
    if(startsWith(fileLines[i], prefix=">") == TRUE){
      
      # Store previous sequence
      if(i != 2){
        sequences[[name]] <- strsplit(sequence, split="")[[1]]
      }
      
      # Get sequence name
      name <- substr(fileLines[i], start=2, stop=nchar(fileLines[i]))
      name <- strsplit(name, split="_")[[1]][1]
      
      # Reset sequence
      sequence <- ""
    }else{
      sequence <- paste(sequence, fileLines[i], sep="")
    }
  }
  
  # Store last sequence
  sequences[[name]] <- strsplit(sequence, split="")[[1]]
  
  return(sequences)
}

addIsolateCoverage <- function(summary, isolateCoverage){
  
  summary$Coverage <- rep(NA, nrow(summary))
  
  for(row in 1:nrow(summary)){
    
    if(is.null(isolateCoverage[[summary[row, "Isolate"]]]) == FALSE){
      
      summary[row, "Coverage"] <- isolateCoverage[[summary[row, "Isolate"]]][1]
    }
  }
  
  return(summary)
}

getIsolateCoverage <- function(coverage){
  
  # Initialise a list to store the isolate coverage information
  isolateCoverage <- list()
  
  for(row in 1:nrow(coverage)){
    
    id <- strsplit(coverage[row, "Isolate"], split="_")[[1]][1]
    isolateCoverage[[id]] <- coverage[row, "Coverage"]
  }
  
  return(isolateCoverage)
}

addAHVLACodes <- function(summary, ahvlaCodes){
  
  for(row in 1:nrow(summary)){
    
    # Metadata spoligotype
    if(is.null(ahvlaCodes[[summary[row, "AHVLAType"]]]) == FALSE){
      
      summary[row, "AHVLAType"] <- paste(summary[row, "AHVLAType"], " (", 
                                         ahvlaCodes[[summary[row, "AHVLAType"]]], ")",
                                         sep="")
    }
    
    # Best Match
    if(is.null(ahvlaCodes[[summary[row, "BestMatch"]]]) == FALSE){
      
      summary[row, "BestMatch"] <- paste(summary[row, "BestMatch"], " (", 
                                         ahvlaCodes[[summary[row, "BestMatch"]]], ")",
                                         sep="")
    }
  }
  
  return(summary)
}

addCattleIsolateSequenceIDs <- function(summary, cattleIsolateSequenceIds){
  
  summary$SequenceID <- rep(NA, nrow(summary))
  
  for(row in 1:nrow(summary)){
    
    # Skip badgers
    if(grepl(x=summary[row, "Isolate"], pattern="WB") == TRUE){
      next
    }
    
    summary[row, "SequenceID"] <- cattleIsolateSequenceIds[[summary[row, "Isolate"]]]
  }
  
  return(summary)
}

noteSequenceIDsOfCattleIsolates <- function(linkTable){
  
  isolateSequenceIds <- list()
  for(row in 1:nrow(linkTable)){
    isolateSequenceIds[[linkTable[row, "Seq.number"]]] <- linkTable[row, "Strain.ID"]
  }
  
  return(isolateSequenceIds)
}

summariseIsolateSpoligotypeMatchingInfo <- function(isolates, metaDataSpoligotypes,
                                                    spoligotypeIndices,
                                                    isolateMatchingInfo,
                                                    averageReadDepthForIsolates,
                                                    proportionNsForIsolates){
  # Initialise a summary table
  summary <- data.frame(Isolate=as.character(isolates), 
                        AHVLAType=as.character(metaDataSpoligotypes), 
                        AHVLATypeInfo=rep(NA, length(isolates)),
                        BestMatch=rep(NA, length(isolates)),
                        BestMatchInfo=rep(NA, length(isolates)),
                        AverageDepth=averageReadDepthForIsolates,
                        ProportionNs=proportionNsForIsolates,
                        stringsAsFactors=FALSE)
  
  for(row in 1:nrow(summary)){
    
    # Get match information for AHVLA assigned spoligotype
    matchInfo <- getMatchingInformationForSpoligotype(summary[row, "AHVLAType"],
                                                      spoligotypeIndices,
                                                      summary[row, "Isolate"],
                                                      isolateMatchingInfo)
    summary[row, "AHVLATypeInfo"] <- paste(matchInfo, collapse=":")
    
    # Find the spoligotype that best matches the WGS data for the isolate
    summary[row, "BestMatch"] <- findBestMatch(summary[row, "Isolate"],
                                               isolateMatchingInfo,
                                               names(spoligotypeIndices))
    
    # Get match information for best match
    matchInfo <- getMatchingInformationForSpoligotype(summary[row, "BestMatch"],
                                                      spoligotypeIndices,
                                                      summary[row, "Isolate"],
                                                      isolateMatchingInfo)
    summary[row, "BestMatchInfo"] <- paste(matchInfo, collapse=":")
  }
  
  return(summary)
}

findBestMatch <- function(isolate, isolateMatchingInfo, spoligotypes){
  
  # Get the isolates matching information
  matchingInfo <- isolateMatchingInfo[[isolate]]
  
  # Examine the matching information for each spoligotype
  bestMatchindex <- -1
  minNMissing <- 999999
  for(i in 1:ncol(matchingInfo)){
    
    # Skip comparisons where mismatches were present
    if(matchingInfo[2, i] != 0){
      next
    }
    
    # Record whether we have found a better match
    if(matchingInfo[1, i] < minNMissing){
      bestMatchIndex <- i
      minNMissing <- matchingInfo[1, i]
    }
  }
  
  return(spoligotypes[bestMatchIndex])
}

getMatchingInformationForSpoligotype <- function(spoligotype, spoligotypeIndices,
                                                 isolate, isolateMatchingInfo){
  
  return(isolateMatchingInfo[[isolate]][, spoligotypeIndices[[spoligotype]]])
}

getTheMatchingInformationForEachIsolate <- function(matchTable, isolates, 
                                                    nColsToIgnore){
  
  # Create a list to store a table which records the nMissing and nMismatches
  isolateMatchingInfo <- list()
  
  # Initialise a matrix to store the matching information
  matchingInfo <- matrix(ncol=ncol(matchTable) - nColsToIgnore, nrow=2)
  
  # Examine each isolate
  for(row in 1:nrow(matchTable)){
    
    # Examine the matching information for each spoligotype
    matchingInfo <- matrix(ncol=ncol(matchTable) - nColsToIgnore, nrow=2)
    for(i in (nColsToIgnore + 1):ncol(matchTable)){
      parts <- as.numeric(strsplit(matchTable[row, i], split=":")[[1]])
      matchingInfo[, i - nColsToIgnore] <- parts
    }
    
    # Store the matching information for the current isolate
    isolateMatchingInfo[[isolates[row]]] <- matchingInfo
  }
  
  return(isolateMatchingInfo)
}

parseFileNames <- function(array){
  
  # TB1385_S1_1.vcf
  output <- c()
  for(index in 1:length(array)){
    
    output[index] <- strsplit(x=array[index], split="_")[[1]][1]
  }
  
  return(output)
}  
  
indexArray <- function(array){
  
  output <- list()
  
  for(index in 1:length(array)){
    output[array[index]] <- index
  }
  
  return(output)
}