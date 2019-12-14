#### Load the cattle spoligotype data ####

# Set the path variable
path <- "/home/josephcrispell/storage/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/"

# Read in the spoligotype table
spoligotypes <- read.table(paste0(path, "Spoligo_APHA_2013-09-09.csv"), header=TRUE, sep=",",
                           stringsAsFactors=FALSE, fill=TRUE)

# Quickly calculate the prevalence of each strain
table(spoligotypes$Genotype) / nrow(spoligotypes)

#### Select data for cattle within 15km of Woodchester Park ####

# Note the coordinates of the badger centre
badgerCentre <- c(381761.7, 200964.3)

# Calculate the distance to the badger centre
spoligotypes <- calculateDistanceToBadgerCentre(badgerCentre, spoligotypes)

# Select only those cattle within X km of Woodchester Park
threshold <- 10000
selected <- spoligotypes[spoligotypes$Distance <= threshold, ]

# Examine how many other spoligotypes were present other than 17 (SB0263)
table(selected$Genotype)
table(selected$Genotype) / nrow(selected)
length(unique(selected$Genotype))

#### Read in the badger spoligotype information ####

# Read in the isolate data table
badgerInfo <- read.table(paste0(path, "IsolateData/BadgerInfo_08-04-15_LatLongs_XY_Centroids.csv"),
                         header=TRUE, sep=",", stringsAsFactors=FALSE)

# Count the spoligotypes observed
table(badgerInfo$AHVLASpoligo)

#### FUNCTIONS ####

calculateDistanceToBadgerCentre <- function(badgerCentre, isolateInfo){
  
  isolateInfo$Distance <- rep(NA, nrow(isolateInfo))
  for(row in 1:nrow(isolateInfo)){
    
    # Skip isolates with an unknown location
    if(is.na(isolateInfo[row, "Mapx"]) == FALSE){
      isolateInfo[row, "Distance"] <- euclideanDistance(x1=badgerCentre[1], y1=badgerCentre[2],
                                                        x2=isolateInfo[row, "Mapx"],
                                                        y2=isolateInfo[row, "Mapy"])
    }  
  }
  
  return(isolateInfo)
}

euclideanDistance <- function(x1, y1, x2, y2){
  return(sqrt((x1 - x2)^2 + (y1 - y2)^2))
}