# Load the ape package
library(ape)
library(geiger) # For the tips function

#########################
# Get a list of Badgers #
#########################

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/"

# Read in the newick tree
file <- paste(path, "InvestigatingCattleMislabelling/",
              "mlTree_Prox-10_plusRef_rmResequenced_SNPCov-0.1_28-10-16.tree", sep="")
tree <- read.tree(file=file)
tips <- tree$tip.label

# Create an array of the ids of the badgers used
badgerIsolates <- tips[grepl(pattern="WB", x=tips)]

###############################
# Get Badger Isolate MetaData #
###############################

# Badger Isolates
file <- paste(path, "IsolateData/",
              "BadgerInfo_08-04-15_LatLongs_XY.csv", sep="")
isolateInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

#################################################
# Get the MetaData for the Badger Isolates Used #
#################################################

# Convert the vector of badger isolates used to a list
listOfBadgerIsolates <- convertVectorToList(badgerIsolates)

# Note which rows of the isolate info table to keep
keep <- c()
index <- 0

for(row in 1:nrow(isolateInfo)){
  
  # Was the current isolate used?
  if(is.null(listOfBadgerIsolates[[isolateInfo[row, "WB_id"]]]) == FALSE){
    index <- index + 1
    keep[index] <- row
  }
}

# Print out the metadata for the badger isolates that were used
file <- paste(path, "BadgerIsolateMetaData_163IsolatesUsed_10-01-17.csv", sep="")
write.table(isolateInfo[keep, ], file=file, quote=FALSE, sep=",", row.names=FALSE)

#############
# FUNCTIONS #
#############

convertVectorToList <- function(vector){
  
  result <- list()
  for(i in 1:length(vector)){
    result[[vector[i]]] <- i
  }
  
  return(result)
}