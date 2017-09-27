#########################################
# Read in the Spoligotyping Tool Output #
#########################################

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_13-07-17/"

# Open the file
file <- paste(path, "Mislabelling/Spoligotyping/",
              "SpoligotypeMatches_27-07-2017.txt", sep="")
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
              "isolateCoverageSummary_DP-20_27-07-2017.txt", sep="")
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

# Add some guesses
misMatches[misMatches$ProportionNs > 0.5, "Reason"] <- "Low Sequencing Quality"
misMatches[which(grepl(x=misMatches$Isolate, 
                       pattern="WB174|TB1775|TB1391|TB1398|TB1783")),
           "Reason"] <- "Close Enough?"
misMatches[which(grepl(x=misMatches$Isolate, 
                       pattern="TB1471|TB1491")),
           "Reason"] <- "Not close to Anything"
misMatches <- misMatches[order(misMatches$Reason, decreasing=TRUE), ]

# Write out the table
file <- paste(path, "Mislabelling/Spoligotyping/",
              "SpoligotypeMatches_Summary_31-01-2017.txt", sep="")
write.table(x=misMatches, file, quote=FALSE, sep="\t", row.names=FALSE)

#############
# Functions #
#############

addIsolateCoverage <- function(summary, isolateCoverage){
  
  summary$Coverage <- rep(NA, nrow(summary))
  
  for(row in 1:nrow(summary)){
    
    if(is.null(isolateCoverage[[summary[row, "Isolate"]]]) == FALSE){
      
      summary[row, "Coverage"] <- isolateCoverage[[summary[row, "Isolate"]]][2]
    }
  }
  
  return(summary)
}

getIsolateCoverage <- function(coverage){
  
  # Initialise a list to store the isolate coverage information
  isolateCoverage <- list()
  
  for(row in 1:nrow(coverage)){
    
    id <- strsplit(coverage[row, "IsolateID"], split="_")[[1]][1]
    isolateCoverage[[id]] <- coverage[row, c(2,3)]
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