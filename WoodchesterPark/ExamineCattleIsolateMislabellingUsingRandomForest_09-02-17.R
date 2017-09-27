#############
# Libraries #
#############

library(ape)
library(geiger)
library(randomForest)

########
# Path #
########

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/"

#################################
# Open the ML Phylogenetic Tree #
#################################

# Read in the newick tree
file <- paste(path, "allVCFs-IncludingPoor/vcfFiles/",
              "mlTree_Prox-10_plusRef_rmResequenced_SNPCov-0.1_28-10-16.tree", sep="")
tree <- read.tree(file=file)

# Remove all the badger isolates
tree <- removeBadgersFromTree(tree)

#################
# Define Clades #
#################

# Plot the tree to define the clades
# file <- paste(path, "testTree_nodeNumbers.pdf")
# pdf(file, height=20, width=20)
# plot.phylo(tree, type="fan", show.tip.label=FALSE)
# nodelabels()
# dev.off()

# Convert branch lengths into SNPs
tree$edge.length <- tree$edge.length * 9464

# Define the clades on the phylogenetic tree
nodesDefiningClades <- c(95, 86, 162, 156, 107)

# Plot the phylogenetic tree with defined clades
cladeColours <- c("red", "blue", "green", "cyan", "darkorchid4")
branchColours <- defineBranchColoursOfClades(tree, nodesDefiningClades, cladeColours, "lightgrey")
plot.phylo(tree, show.tip.label=FALSE, "fan",
           edge.color=branchColours, edge.width=3)
text(x=c(66.34555,137.01424,136.34755,91.67961,-126.32662),
     y=c(126.6943,56.02561,-41.97719,-71.97804,-18.64319),
     labels=c("1", "2", "3", "4", "5"),
     col=cladeColours, cex=2)

# Note the clades that each isolate is found in
isolateClades <- noteIsolateClades(tree, nodesDefiningClades)

#########################
# Build Predictor Table #
#########################

## Get the sampling information
file <- paste(path, "IsolateData/",
              "CattleIsolateInfo_LatLongs_plusID_outbreakSize_Coverage_AddedTB1453-TB1456.csv", sep="")
samplingInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)
samplingInfo$IsolateID <- samplingInfo$StrainId

# Make sure all County names are upper case
samplingInfo$County <- toupper(samplingInfo$County)

# Split the breakdown ID into the CPH and date
samplingInfo <- splitBreakdownIDIntoCPHAndDate(samplingInfo)

## Get the coverage information
file <- paste(path, "allVCFs-IncludingPoor/vcfFiles/",
              "isolateGenomeCoverageSummary_28-10-16.txt", sep="")
coverageTable <- read.table(file, header=TRUE, stringsAsFactors=FALSE)
coverageTable$IsolateID <- getIsolateIDFromFileNames(coverageTable$IsolateID)
coverageTable <- coverageTable[grepl(pattern="TB", x=coverageTable$IsolateID) == TRUE, ]

## Combine all of the above into a single predictor table
indexedIsolates <- indexArray(names(isolateClades))
predictorTable <- data.frame(IsolateID=as.character(names(isolateClades)))

# Sampling information
columnsToKeep <- c("DateCultured","ReasonForSlaughter","SkinTestType","Mapx","Mapy","LesionsFound",
                   "BreakdownCPH", "BreakdownDate","County","Genotype","ComplVNTR","OutbreakSize")
predictorTable <- addColumnsToPredictorTable(predictorTable=predictorTable, table=samplingInfo,
                                             columnsToAdd=columnsToKeep, linkColumn="IsolateID",
                                             indexedIsolates=indexedIsolates)

# Coverage information
columnsToKeep <- c("PercentageCoverage")
predictorTable <- addColumnsToPredictorTable(predictorTable=predictorTable, table=coverageTable,
                                             columnsToAdd=columnsToKeep, linkColumn="IsolateID",
                                             indexedIsolates=indexedIsolates)

#############################
# Fit a Random Forest Model #
#############################

## Use the isolate clades as the response
cladesIsolatesFoundIn <- as.factor(getOrderedIsolatesClades(isolateClades, indexedIsolates))
predictorTable$Clade <- as.factor(cladesIsolatesFoundIn)

## Remove the Isolate ID column
isolateIds <- predictorTable$IsolateID
predictorTable <- predictorTable[, grepl(pattern="IsolateID", colnames(predictorTable)) == FALSE]

## Replace empty cells with NA
predictorTable[predictorTable == ""] <- "Unknown"

## Replace NA Outbreak size values
predictorTable[is.na(predictorTable$OutbreakSize), "OutbreakSize"] <- -1

## Ensure each predictor is in the correct format
predictorTable$DateCultured <- as.Date(predictorTable$DateCultured, format="%d/%m/%Y")
predictorTable$ReasonForSlaughter <- as.factor(predictorTable$ReasonForSlaughter)
predictorTable$SkinTestType <- as.factor(predictorTable$SkinTestType)
predictorTable$Mapx <- as.numeric(predictorTable$Mapx)
predictorTable$Mapy <- as.numeric(predictorTable$Mapy)
predictorTable$LesionsFound <- as.factor(predictorTable$LesionsFound)
predictorTable$BreakdownCPH <- as.factor(predictorTable$BreakdownCPH)
predictorTable$BreakdownDate <- as.Date(predictorTable$BreakdownDate, format="%Y-%m-%d")
predictorTable$County <- as.factor(predictorTable$County)
predictorTable$Genotype <- as.factor(predictorTable$Genotype)
predictorTable$ComplVNTR <- as.factor(predictorTable$ComplVNTR)
predictorTable$OutbreakSize <- as.numeric(predictorTable$OutbreakSize)
predictorTable$PercentageCoverage <- as.numeric(predictorTable$PercentageCoverage)

## Fit a Random Forest model

# Tune the model
cols <- c(1:ncol(predictorTable))[grepl(pattern="Clade", colnames(predictorTable)) == FALSE]
initialMtry <- 3
nTrees <- 1000
tuneOutput <- tuneRF(predictorTable[, cols],
                     predictorTable$Clade, mtryStart=initialMtry,
                     ntreeTry=nTrees, stepFactor=1.5, improve=0.0001, 
                     trace=TRUE, plot=TRUE)
mTry <- getOptimalMtry(tuneOutput)

# Set the input parameters
nTrees <- 10000

# Run the model
infoRF <- randomForest(Clade~., data=predictorTable,
                       proximity=FALSE, mtry=mTry, importance=TRUE,
                       ntree=nTrees, do.trace=FALSE, keep.forest=TRUE,
                       norm.votes=FALSE)
plot(infoRF)

# Plot the variable importance
varImpPlot(infoRF, cex=0.5)

# Get the isolate assignments
isolateAssignmentVotes <- as.data.frame(infoRF$votes)
isolateAssignmentVotes$Assignment <- getIsolatePredictions(infoRF$predicted)
isolateAssignmentVotes$IsolateID <- isolateIds
isolateAssignmentVotes$Clade <- predictorTable$Clade

# Get the isolates whose assignment didn't match their actual clade
poorAssignment <- isolateAssignmentVotes[isolateAssignmentVotes$Assignment !=
                                           isolateAssignmentVotes$Clade, ]

####################################################
# Fit Random Forest Model Without Spoligotype Data #
####################################################

# Remove the spoligotype column
predictorTable_NoSpoligotype <- 
  predictorTable[, grepl(pattern="Genotype", x=colnames(predictorTable)) == FALSE]

# Tune the model
cols <- c(1:ncol(predictorTable))[grepl(pattern="Clade", colnames(predictorTable)) == FALSE]
initialMtry <- 3
nTrees <- 1000
tuneOutput <- tuneRF(predictorTable[, cols],
                     predictorTable$Clade, mtryStart=initialMtry,
                     ntreeTry=nTrees, stepFactor=1.5, improve=0.0001, 
                     trace=TRUE, plot=TRUE)
mTry <- getOptimalMtry(tuneOutput)

# Set the input parameters
nTrees <- 10000

# Run the model
infoRF_NoSpoligotype <- randomForest(Clade~., data=predictorTable_NoSpoligotype,
                       proximity=FALSE, mtry=mTry, importance=TRUE,
                       ntree=nTrees, do.trace=FALSE, keep.forest=TRUE,
                       norm.votes=FALSE)
plot(infoRF_NoSpoligotype)

# Plot the variable importance
varImpPlot(infoRF_NoSpoligotype, cex=0.5)

# Get the isolate assignments
isolateAssignmentVotes_NoSpoligotype <- as.data.frame(infoRF_NoSpoligotype$votes)
isolateAssignmentVotes_NoSpoligotype$Assignment <- getIsolatePredictions(infoRF_NoSpoligotype$predicted)
isolateAssignmentVotes_NoSpoligotype$IsolateID <- isolateIds
isolateAssignmentVotes_NoSpoligotype$Clade <- predictorTable_NoSpoligotype$Clade

# Get the isolates whose assignment didn't match their actual clade
poorAssignment_NoSpoligotype <- isolateAssignmentVotes_NoSpoligotype[
                                    isolateAssignmentVotes_NoSpoligotype$Assignment !=
                                    isolateAssignmentVotes_NoSpoligotype$Clade, ]


#############
# FUNCTIONS #
#############

getOptimalMtry <- function(tuneOutput){
  
  min <- 9999999
  minIndex <- -1
  for(row in 1:nrow(tuneOutput)){
    
    if(tuneOutput[row, 2] < min){
      min <- tuneOutput[row, 2]
      minIndex <- row
    }
  }
  
  return(tuneOutput[minIndex, 1])
}

getIsolatePredictions <- function(predicted){
  isolatePredictions <- c()
  for(index in 1:length(predicted)){
    isolatePredictions[index] <- as.numeric(predicted[[index]])
  }
  
  return(isolatePredictions)
}

getOrderedIsolatesClades <- function(isolateClades, indexedIsolates){
  
  cladesIsolatesFoundIn <- c()
  
  for(isolate in names(isolateClades)){
    
    cladesIsolatesFoundIn[indexedIsolates[[isolate]]] <- isolateClades[[isolate]]
  }
  
  return(cladesIsolatesFoundIn)
}

splitBreakdownIDIntoCPHAndDate <- function(samplingInfo){
  
  samplingInfo$BreakdownDate <- rep(NA, nrow(samplingInfo))
  samplingInfo$BreakdownCPH <- rep(NA, nrow(samplingInfo))
  
  for(row in 1:nrow(samplingInfo)){
    
    parts <- strsplit(samplingInfo[row, "BreakdownID"], split="-")[[1]]
    
    samplingInfo[row, "BreakdownDate"] <- as.character(as.Date(parts[2], format="%d/%m/%Y"))
    samplingInfo[row, "BreakdownCPH"] <- parts[1]
  }
  
  return(samplingInfo)
}

addColumnsToPredictorTable <- function(predictorTable, table, columnsToAdd, linkColumn, indexedIsolates){
  
  # Add empty columns to fill
  predictorTable[, columnsToAdd] <- rep(NA, nrow(predictorTable))
  
  # Convert the IsolateID column to a character column
  table[, linkColumn] <- as.character(table[, linkColumn])
  
  # Examine each row of the table which the columns are coming from
  for(row in 1:nrow(table)){
    
    # Check if the information is for an isolate we're interested in
    if(is.null(indexedIsolates[[table[row, linkColumn]]]) == FALSE){
      
      predictorTable[indexedIsolates[[table[row, linkColumn]]], columnsToAdd] <- table[row, columnsToAdd]
    }
  }
  
  return(predictorTable)
}

indexArray <- function(array){
  
  output <- list()
  for(i in 1:length(array)){
    output[[array[i]]] <- i
  }
  
  return(output)
}

countNumberOfMovements <- function(lifeHistoryTable, premisedToIgnore){
  
  # Initialise an array to store the number of movements
  numberMovements <- c()
  
  # Examine each of the sampled cattle
  for(row in 1:nrow(lifeHistoryTable)){
    
    # Initialise a variable to record the number of movements made by each cow
    numberMovements[row] <- 0
    
    # Get the premises types that animal lived on
    premisesTypes <- strsplit(lifeHistoryTable[row, "PremisesTypes"], split=",")[[1]]
    
    # Check that the premises aren't ones we want to ignore
    for(premises in premisesTypes){
      if(is.null(premisedToIgnore[[premises]]) == TRUE){
        numberMovements[row] <- numberMovements[row] + 1
      }
    }
  }
  
  return(numberMovements)
}

getIsolateIDFromFileNames <- function(fileNames){
  isolates <- c()
  for(i in 1:length(fileNames)){
    isolates[i] <- strsplit(fileNames[i], split="_")[[1]][1]
  }
  
  return(isolates)
}


noteIsolateClades <- function(tree, nodesDefiningClades){
  
  # Initialise a list to note the clade of each isolate
  isolateClades <- list()
  
  # Examine each clade
  for(cladeIndex in 1:length(nodesDefiningClades)){
    
    # Get the tips associated with the current clade
    cladeTips <- tips(tree, node=nodesDefiningClades[cladeIndex])
    
    # Note the isolates in the current clade
    for(tip in cladeTips){
      isolateClades[[tip]] <- cladeIndex
    }
  }
  
  return(isolateClades)
}

defineBranchColoursOfClades <- function(tree, nodesDefiningClades,
                                        CladeColours, defaultColour){
  branchColours <- rep(defaultColour, dim(tree$edge)[1])
  for(i in 1:length(cladeColours)){
    clade <- tips(tree, node=nodesDefiningClades[i])
    branchesInClades <- which.edge(tree, clade)
    branchColours[branchesInClades] <- cladeColours[i]
  }
  return(branchColours)
}

removeBadgersFromTree <- function(tree){
  
  # Examine each of the tips on the phylogenetic tree
  for(tip in tree$tip.label){
    
    # Is the current tip for a badger isolate?
    if(grepl(pattern="WB", x=tip) == TRUE){
      tree <- drop.tip(tree, tip)
    }
  }
  
  return(tree)
}