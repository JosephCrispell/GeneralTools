path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/"

#######################
# Herd Size Information
file <- paste(path, "CattleMovementData/cts_locations_herdSizes_2012-07-01.csv", sep="")
herdSizeInfo <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep=",", fill=TRUE)

###########################
# Herd Location Information
file <- paste(path, "CattleMovementData/20160314_joe_cts_locations.csv", sep="")
locationInfo <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep=",", fill=TRUE)

# Get a list of herds that are within threshold distance of Woodchester Park
mansionX <- 380909
mansionY <- 201377
threshold <- 15000

herdsAroundWP <- noteLocationsWithinThresholdDistance(locationInfo=locationInfo, 
                                                      thresholdInMetres=threshold,
                                                      mansionX=mansionX,
                                                      mansionY=mansionY)

#####################################################
# Get the herd size for those around Woodchester Park

herdsAroundWP <- getHerdSizesForHerds(herdSizeInfo, herdsAroundWP)

#################################
# Plot the herd size distribution

herdSizes <- as.numeric(getHerdSizes(herdsAroundWP))



file <- paste(path, "HerdSizes_28-03-17.pdf")
pdf(file)

par(mar=c(5.1,5.1,4.1,2.1)) # Bottom, Left, Top, Right

hist(herdSizes, breaks=100, main = "Herd Size around Woodchester Park", las=1, 
     ylab="Frequency", xlab="Number Cattle Present", cex.main=2, cex.axis=1.5, cex.lab=1.5)

dev.off()

#############
# FUNCTIONS #
#############

getHerdSizes <- function(herdsSizeInfo){
  
  # Get the location Ids
  locationIds <- names(herdsSizeInfo)
  
  # Initialise an array to store the herd sizes
  herdSizes <- c()
  index <- 0
  
  count <- 0
  
  # Note the herd sizes
  for(key in locationIds){
    index <- index + 1
    
    # Check that herd size information is available
    if(grepl(x=herdsSizeInfo[[key]], pattern="/") == TRUE){
      cat(paste("\"", key, "\",\"", herdsSizeInfo[[key]], "\"\n", sep=""))
    }else{
      count <- count + 1
      herdSizes[index] <- herdsSizeInfo[[key]]
    }
  }

  print(count)
  
  return(herdSizes)
}

getHerdSizesForHerds <- function(herdSizeInfo, herdsAroundWP){
  
  # Examine the herd size information
  for(row in 1:nrow(herdSizeInfo)){
    
    # Is the herd present around WP?
    if(is.null(herdsAroundWP[[as.character(herdSizeInfo[row, "location_id"])]]) == FALSE){
      herdsAroundWP[[as.character(herdSizeInfo[row, "location_id"])]] <- herdSizeInfo[row, "num_cattle"]
    }
  }
  
  return(herdsAroundWP)
}

noteLocationsWithinThresholdDistance <- function(locationInfo, thresholdInMetres, 
                                                 mansionX, mansionY){
  
  locationInfo$DistanceToWoodchester <- rep(0, nrow(locationInfo))
  
  # Initialise an array to record the rows to keep
  herdsToKeep <- list()

  # Examine each location
  for(row in 1:nrow(locationInfo)){
    
    # Skip if no location information available
    if(is.na(locationInfo[row, "x"]) == TRUE || is.na(locationInfo[row, "y"]) == TRUE){
      next
    }
    
    # Calculate distance to Woodchester Mansion
    distance <- euclideanDistance(x1=mansionX, y1=mansionY, 
                                  x2=locationInfo[row, "x"], y2=locationInfo[row, "y"])

    # Keep herd if distance is <=threshold
    if(distance <= thresholdInMetres){
      herdsToKeep[[as.character(locationInfo[row, "location_id"])]] <- locationInfo[row, "cph"]
    }
    
    if(row %% 10000 == 0){
      print(paste("Finished reading row ", row, " of ", nrow(locationInfo), sep=""))
    }
  }
  
  return(herdsToKeep)
}

euclideanDistance <- function(x1, y1, x2, y2){
  return(sqrt(sum((x1 - x2)^2 + (y1 - y2)^2)))
}
