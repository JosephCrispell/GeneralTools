#################
# Load Packages #
#################
suppressMessages(library(randomForest))

#####################################################
# Read in the Genetic Vs. Epidemiological Distances #
#####################################################

# Get the path to the necessary files
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_13-07-17/"

# Read Genetic V.s Epi Distances table
file <- paste(path, "Mislabelling/Badger-RF-BR/",
              "geneticVsEpiTable_28-09-2017.txt", sep="")
table <- read.table(file, header=TRUE)

#################################
# Examine the Genetic Distances #
#################################

# Set threshold 
threshold <- 15

par(mfrow=c(3,1))

hist(table$GeneticDistance, 
     las=1,
     xlab="Genetic Distance (SNPs)",
     main="Inter-Isolate Genetic Distance Distribution")

hist(table[table$GeneticDistance < 50, ]$Genetic, 
     las=1,
     xlab="Genetic Distance (SNPs)",
     main="Inter-Isolate Genetic Distance Distribution < 50 SNPs")
lines(x=c(threshold, threshold), y=c(0, 8500), col="red", lty=2)

hist(table[table$GeneticDistance < threshold, ]$Genetic, 
     las=1,
     xlab="Genetic Distance (SNPs)",
     main=paste("Inter-Isolate Genetic Distance Distribution < ",
                threshold, "15 SNPs", sep=""))

par(mfrow=c(1,1))

#######################################
# Select only Small Genetic Distances #
#######################################

table <- table[table$GeneticDistance < threshold, ]

#########################
# Load the Isolate Data #
#########################

# Badger Isolate Information
badgerInfoFile <- paste(path, "IsolateData/",
                        "BadgerInfo_08-04-15_LatLongs_XY_Centroids.csv",
                        sep="")
badgerInfo <- read.table(badgerInfoFile, header=TRUE, sep=",")
rownames(badgerInfo) <- badgerInfo[,1]

###########################################################
# Fit the Random Forest Model to test Variation Explained #
###########################################################

# Split the data into train and test sets
trainProp <- 0.5
trainRows <- sample(x=1:nrow(table),
                    size=floor(trainProp * nrow(table)), replace=FALSE)

# Define the columns of the table we're interested in
cols <- 2:26
nTrees <- 1000

# Run the Random Forest tuning algorithm to find the optimal mtry value
tuneOutput <- runRandomForestTuning(table[trainRows, ], 3, nTrees, cols)
plot(tuneOutput, las=1, type="o")

 # Calculate the mtry from the tuning output
mTry <- as.integer(rownames(tuneOutput)[tuneOutput[,2] == min(tuneOutput[,2])])
points(x=mTry, y=tuneOutput[which(tuneOutput[, 1] == mTry), 2], col="red",
       pch=20)

# Run Random Forest Algorithm
nTrees <- 1000
infoRF <- randomForest(table[trainRows, "GeneticDistance"]~., data=table[trainRows, cols],
                       proximity=FALSE, mtry=mTry, importance=TRUE,
                       ntree=nTrees, do.trace=FALSE, keep.forest=TRUE,
                       norm.votes=FALSE)
plot(infoRF, las=1)

# Add the predicted values onto the table
predictions <- predict(infoRF, table[-trainRows, cols])

# Get the Pseudo RSquared value
rSq <- round(infoRF$rsq[length(infoRF$rsq)], digits=2)

###################################################
# Fit the Random Forest Model to find Mislabelled #
###################################################

# Define the columns of the table we're interested in
cols <- 2:26
nTrees <- 1000

# Run the Random Forest tuning algorithm to find the optimal mtry value
tuneOutput <- runRandomForestTuning(table, 3, nTrees, cols)
plot(tuneOutput, las=1, type="o")

# Calculate the mtry from the tuning output
mTry <- as.integer(rownames(tuneOutput)[tuneOutput[,2] == min(tuneOutput[,2])])
points(x=mTry, y=tuneOutput[which(tuneOutput[, 1] == mTry), 2], col="red",
       pch=20)

# Run Random Forest Algorithm
nTrees <- 1000
infoRF <- randomForest(table$GeneticDistance~., data=table[, cols],
                         proximity=FALSE, mtry=mTry, importance=TRUE,
                         ntree=nTrees, do.trace=FALSE, keep.forest=FALSE,
                         norm.votes=FALSE)
plot(infoRF, las=1)

# Add the predicted values onto the table
table$predictions <- infoRF$predicted

# Calculate the mean Difference between the actual and predicted values for each isolate
meanValues <- examinePredictedVersusActual(table, badgerInfo)


#################
# Plot the Data #
#################

# Open an output PDF file
file <- paste(path, "Mislabelling/Badger-RF-BR/",
              "FittedRandomForest_28-09-17.pdf", sep="")
pdf(file)

par(mfrow=c(1,1))

# Plot the forest building output
plot(infoRF, las=1)

# Predicted Versus Actual
smoothScatter(predictions, table[-trainRows, "GeneticDistance"],
              main="Actual vs. Predicted",
              cex.main=1,
              xlab="Predicted",
              ylab="Actual",
              nrpoints=0, las=1)
abline(lm(table$GeneticDistance ~ table$predictions), col="red")
correlation <- round(cor(table[-trainRows, "GeneticDistance"], predictions), digits=2)
legend("topleft", legend=c(paste("corr =", correlation), paste("Rsq =", rSq)), bty="n", cex = 1)

# Isolate Prediction
mislabelled <- c("WB100", "WB71", "WB105", "WB107", "WB72")
plotAgainstTime(meanValues, 2, "Isolate Prediction", 
                "Median Difference",
                cex.main=1, cex.leg=1, mislabelled)



dev.off()

#############################
# Store Isolate Predictions #
#############################

### Print the mean values table for the Badger-Badger comparisons to file
outputTable <- addAnimalIdentifierColumnToMeanValuesTable(meanValues, badgerInfo)

# Add WB ID to table
outputTable$WBID <- rownames(outputTable)

# Re-order table columns
outputTable <- outputTable[, c(ncol(outputTable), ncol(outputTable)-1, 1:(ncol(outputTable) - 2))]
file <- paste(path, "Mislabelling/Badger-RF-BR/", 
              "MeanValuesTable_RandomForest_28-09-2017.csv", sep="")
write.table(outputTable, file, quote=FALSE, row.names=FALSE, sep=",")

#############
# FUNCTIONS #
#############

runRandomForestTuning <- function(inputTable, initialMtry, nTrees, cols){
  
  tuneOutput <- tuneRF(inputTable[, cols], inputTable$Genetic, mtryStart=initialMtry,
                       ntreeTry=nTrees, stepFactor=1.5, improve=0.0001, trace=TRUE,
                       plot=TRUE)
  
  return(tuneOutput)
}

addAnimalIdentifierColumnToMeanValuesTable <- function(meanValues, badgerInfo){
  
  # Initialise an additional column on the meanValues table
  meanValues$Identifier <- rep(NA, nrow(meanValues))
  
  # Index the isolates in the mean Values table
  isolateIndices <- list()
  
  # Examine each row name of the meanValues table
  rowNames <- rownames(meanValues)
  for(i in 1:nrow(meanValues)){
    
    isolateIndices[[rowNames[i]]] <- i
  }
  
  # Examine the badger info
  rowNames <- rownames(badgerInfo)
  for(row in 1:nrow(badgerInfo)){
    
    # Does the current isolate exist in the mean values table
    if(is.null(isolateIndices[[rowNames[row]]]) == FALSE){
      
      # Add the current badgers tattoo into the meanValues table
      meanValues[isolateIndices[[rowNames[row]]],ncol(meanValues)] <- as.character(badgerInfo[row, 4])
    }
  }
  
  return(meanValues)
}

plotAgainstTime <- function(inputTable, column, title, yLabel, cex.main, cex.leg,
                            highlight){
  
  ## Plotting the Mean Difference between the predicted and actual genetic distances for 
  #  each isolate
  
  # Input table structure:
  # Mean  Median  Lower Upper Range Year  Species
  # 1     2       3     4     5     6     7
  
  # Find the Boundary where 95% of data lie below
  bound95<- quantile(inputTable[,column], c(0.95))
  
  # Highlight and Label Outliers
  plot(inputTable$Date, inputTable[,column],
       #col=ifelse(inputTable[,column] > bound95,'red','black'),
       col=ifelse(inputTable$Spoligotype == "17", "blue", "red"),
       ylab=yLabel,
       main=title, cex.main=cex.main,
       pch=20,
       xlab="Year", las=1)
  
  with(subset(inputTable, inputTable[,column] >= bound95),
       text(inputTable$Date[inputTable[,column] >= bound95]+150,
            inputTable[,column][inputTable[,column] >= bound95],
            rownames(inputTable)[inputTable[,column] >= bound95], cex=0.5))
  abline(h=bound95, col="red")
  legend("topright", paste("Upper = ", round(bound95, digits=3)), bty="n",
         cex=cex.leg)
  legend("right", legend=c("17", "Other"),
         text.col=c("blue", "red"),
         bty="n", cex=cex.leg)
  
  # Additionally highlight the list of isolate IDs
  indices <- which(rownames(inputTable) %in% highlight)
  for(index in indices){
    points(x=inputTable[index, 6], y=inputTable[index, column], 
           pch=0, cex=2)
  }
  
  print(rownames(inputTable)[inputTable[,column] >= bound95])
  
}

examinePredictedVersusActual <- function(inputTable, badgerInfo){
  
  ## Store a distribution of the difference between the predicted
  ## and actual genetic distances for each isolate
  
  # Input table structure:
  # GeneticDistance iSpeciesJSpecies  EpiMetricA  EpiMetricB  ... IsolateI  IsolateJ  Predicted
  # 0               1                 2           3           ... -2        -1        ncol(table)
  
  # Ensure the Isolate ID columns are vectors of strings
  inputTable$IsolateI <- as.character(inputTable$IsolateI)
  inputTable$IsolateJ <- as.character(inputTable$IsolateJ)
  
  # Initialise a list to store the difference distributions for each isolate
  isolates <- list()
  
  # Examine each row of the input table
  for(row in 1:nrow(inputTable)){
    
    ## Check the I Isolate ID
    # Have we envountered this isolate before?
    if(is.null(isolates[[inputTable[row, ncol(inputTable) - 2]]]) == TRUE){
      
      # Calculate the difference between the predicted and actual values
      difference <- abs(inputTable[row, ncol(inputTable)] - inputTable[row, 1])
      
      # Store the array of calculated differences
      isolates[[inputTable[row, ncol(inputTable) - 2]]] <- c(difference)
      
      # We have encountered this isolate - append the difference
    }else{
      
      # Calculate the difference between the predicted and actual values
      difference <- abs(inputTable[row, ncol(inputTable)] - inputTable[row, 1])
      
      # Append the calculated difference to the array
      isolates[[inputTable[row, ncol(inputTable) - 2]]] <- c(isolates[[inputTable[row, ncol(inputTable) - 2]]], difference)
    }
    
    ## Check the J Isolate ID
    # Have we envountered this isolate before?
    if(is.null(isolates[[inputTable[row, ncol(inputTable) - 1]]]) == TRUE){
      
      # Calculate the difference between the predicted and actual values
      difference <- abs(inputTable[row, ncol(inputTable)] - inputTable[row, 1])
      
      # Store the array of calculated differences
      isolates[[inputTable[row, ncol(inputTable) - 1]]] <- c(difference)
      
      # We have encountered this isolate - append the difference  
    }else{
      
      # Calculate the difference between the predicted and actual values
      difference <- abs(inputTable[row, ncol(inputTable)] - inputTable[row, 1])
      
      # Append the calculated difference to the array
      isolates[[inputTable[row, ncol(inputTable) - 1]]] <- c(isolates[[inputTable[row, ncol(inputTable) - 1]]], difference)
    }
  }
  
  ## Prepare the badger sampling information to be added into the output table
  
  # Badger Information table structure:
  # WB_id	CB_id	Batch	tattoo	date	pm	sample	AHVLA_afno	lesions	abscesses	AHVLASpoligo
  # 1     2     3     4       5     6   7       8           9       10        11
  #
  # Social.Group.Trapped.At	AFBI_VNTRNo	AFBI_String	AFBI_Genotype	AFBI_Spoligotype	AFBI_GenSpol
  # 12                      13          14          15            16                17
  #
  # notes	SampledGrpLat	SampledGrpLong
  # 18    19            20
  
  # Convert the tables into lists - key = ID
  badgers <- setNames(split(badgerInfo, seq(nrow(badgerInfo))), rownames(badgerInfo))
  
  ## Summarise the distributions of differences for each isolate
  
  # Initialise a table to store a summary of the difference distributions for each isolate
  summaryTable <- data.frame(Mean=rep(0, length(isolates)),
                             Median=rep(0, length(isolates)),
                             Lower=rep(0, length(isolates)),
                             Upper=rep(0, length(isolates)),
                             Range=rep(0, length(isolates)),
                             Date=rep(as.Date("19/09/1900","%d/%m/%Y"), length(isolates)),
                             Spoligotype=rep(0, length(isolates)),
                             NumberDistances=rep(0, length(isolates)))
  
  # Get a list of the isolate IDs
  ids <- names(isolates)
  rownames(summaryTable) <- ids
  
  # Examine each isolates difference distribution
  for(i in 1:nrow(summaryTable)){
    
    # Get the distribution of differences for the current isolate
    differences <- isolates[[ids[i]]]
    
    # Calculate the quantiles of the distribution
    quantiles <- quantile(differences, c(0.025, 0.975))
    
    # Summarise the difference distribution
    summaryTable[i,1] <- mean(differences)
    summaryTable[i,2] <- median(differences)
    summaryTable[i,3] <- quantiles[[1]]
    summaryTable[i,4] <- quantiles[[2]]
    summaryTable[i,5] <- quantiles[[2]] - quantiles[[1]]
    summaryTable[i,8] <- length(differences)
    
    # Add in some additional sampling information
    summaryTable[i,6] <- as.Date(as.character(badgers[[ids[i]]][[5]]), "%d/%m/%Y")
    
    spoligotype <- as.character(badgers[[ids[i]]][[11]])
    if(spoligotype == "RETEST" | spoligotype == ""){
      spoligotype = "NA"
    }
    summaryTable[i,7] <- spoligotype
    
  }
  
  # Change the class of the spoligotype column
  summaryTable$Spoligotype <- as.factor(summaryTable$Spoligotype)
  
  return(summaryTable)
}