library(gbm)
library(dismo)

#####################################################
# Read in the Genetic Vs. Epidemiological Distances #
#####################################################

# Get the path to the necessary files
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/"

# Read Genetic V.s Epi Distances table
file <- paste(path, "FindMislabelled/GeneticVsEpidemiologicalDistances/",
              "geneticVsEpiTable_15-08-2016.txt", sep="")
table <- read.table(file, header=TRUE)

#################################
# Examine the Genetic Distances #
#################################

# Open an output PDF file
file <- paste(path, "FindMislabelled/GeneticVsEpidemiologicalDistances/",
              "GeneticDistanceDistribution_26-05-16.pdf", sep="")
pdf(file)

par(mfrow=c(2,1))

hist(table$GeneticDistance, breaks=150,
     las=1,
     xlab="Genetic Distance (SNPs)",
     main="Inter-Isolate Genetic Distance Distribution")
lines(x=c(15, 15), y=c(0, 8500), col="red", lty=2)

hist(table[table$GeneticDistance < 15, ]$GeneticDistance, breaks=10,
     las=1,
     xlab="Genetic Distance (SNPs)",
     main="Inter-Isolate Genetic Distance Distribution")

dev.off()

#######################################
# Select only Small Genetic Distances #
#######################################

table <- table[table$GeneticDistance < 15, ]

########################################################
# Fit a Boosted Regression to test Variation Explained #
########################################################

# Split the data into train and test sets
trainProp <- 0.5
trainRows <- sample(x=1:nrow(table),
                    size=floor(trainProp * nrow(table)), replace=FALSE)

# Define the columns to be used
colNames <- colnames(table)
cols <- 2:25

# Run Model
stepOutput <- gbm.step(data = table[trainRows, ],
                       gbm.x = colNames[cols], gbm.y = "GeneticDistance",
                       family = "poisson", tree.complexity = 5,
                       learning.rate = 0.01, bag.fraction = 0.5, max.trees = 100000)

# Use the Model to predict
predictions <- predict(stepOutput,
                       newdata = table[-trainRows, cols],
                       n.trees = stepOutput$n.trees,
                       type = "response")

################################################
# Fit a Boosted Regression to find Mislabelled #
################################################

# Define the columns to be used
colNames <- colnames(table)
cols <- 2:25

# Run Model
stepOutput <- gbm.step(data = table,
                       gbm.x = colNames[cols], gbm.y = "GeneticDistance",
                       family = "poisson", tree.complexity = 5,
                       learning.rate = 0.01, bag.fraction = 0.5, max.trees = 100000)

###############################
# Examine Isolate Predictions #
###############################

# Use the Model to predict
table$predictions <- predict(stepOutput,
                             newdata = table[, cols],
                             n.trees = stepOutput$n.trees,
                             type = "response")

# Badger Isolate Information
badgerInfoFile <- paste(path, "IsolateData/BadgerInfo_08-04-15_LatLongs.csv", sep="")
badgerInfo <- read.table(badgerInfoFile, header=TRUE, sep=",")
rownames(badgerInfo) <- badgerInfo[,1]

meanValues <- examinePredictedVersusActual(table, badgerInfo)


#################
# Plot the Data #
#################

# Open an output PDF file
file <- paste(path, "FindMislabelled/GeneticVsEpidemiologicalDistances/",
              "FittedBoostedRegression_26-05-17.pdf", sep="")
pdf(file)

par(mfrow=c(1,1))

# Predicted Versus Actual
smoothScatter(predictions, table[-trainRows, "GeneticDistance"],
              main="Actual vs. Predicted",
              cex.main=1,
              xlab="Predicted",
              ylab="Actual",
              nrpoints=0)
abline(lm(table$GeneticDistance ~ table$predictions), col="red")
correlation <- round(cor(table[-trainRows, "GeneticDistance"], predictions), digits=2)
legend("topleft", paste("corr =", correlation), bty="n", cex = 1)

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
file <- paste(path, "FindMislabelled/GeneticVsEpidemiologicalDistances/", 
              "MeanValuesTable_BoostedRegression_26-05-2016.csv", sep="")
write.table(outputTable, file, quote=FALSE, row.names=FALSE, sep=",")


#####################
# FUNCTIONS SECTION #
#####################

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