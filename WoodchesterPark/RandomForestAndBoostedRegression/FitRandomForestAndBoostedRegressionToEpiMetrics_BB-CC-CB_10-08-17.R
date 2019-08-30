################################################################
# Issues with these genetic vs. epidemiological metrics:       #
# - Presence of -1 when data not available                     #
# - Flattening a matrix - multiple data from individual animal #
################################################################

#### Load libraries ####

library(randomForest)
library(gplots)
library(gbm)
library(dismo)
library(codetools) # To fix weird error when loading RData file
library(plotteR)

#### Open the Genetic Vs Epidemiological distances table ####

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

# Set the path
path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/GeneticVsEpidemiologicalDistances/"

# Read in the table
file <- paste(path, "GeneticVsEpidemiologicalDistances_05-04-18.txt", sep="")
geneticVsEpi <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Add in random noise variables
geneticVsEpi$RandomUniform <- runif(nrow(geneticVsEpi))
geneticVsEpi$RandomBoolean <- sample(c(0,1), size=nrow(geneticVsEpi), replace=TRUE)

#### General settings ####

selection <- "CC" # CC, BB, CB
trainProp <- 0.5
colToUse <- "%IncMSE"
geneticDistanceThreshold <- 15

# Drop out genetic relatedness variable and same animal
if(selection == "BB"){
  cols <- which(colnames(geneticVsEpi) == "HostRelatedness" | colnames(geneticVsEpi) == "SameAnimal")
  geneticVsEpi <- geneticVsEpi[, -cols]
}

# Note the full names of metrics and assign them a colour
fullNames <- noteFullNames(selection)
temporalCol="darkgoldenrod4"
spatialCol="red"
networkCol="blue"
nameColours <- assignMetricColours(temporalCol=temporalCol, spatialCol=spatialCol,
                                   networkCol=networkCol)

#### Select data ####

###### Open a PDF
file <- paste0(path, "ExamineEpiVariableCorrelation_", selection, "_", date, ".pdf")
pdf(file, height=10, width=10)

par(mfrow=c(1,1))

# Subset out the selected comparisons
geneticVsEpi <- selectAppropriateComparisonsForSelection(selection, geneticVsEpi)

# Only select small genetic distances
geneticVsEpi <- selectGeneticDistancesBelowThreshold(threshold=geneticDistanceThreshold)

# Remove irrelevant columns
geneticVsEpi <- removeColumnsIfNotRelevant(geneticVsEpi)

# Convert the columns dealing with boolean metrics to factors
geneticVsEpi <- makeBooleanColumnsFactors(table=geneticVsEpi)

# Note columns to ignore as predictors
colNamesToIgnore <- c("iSpeciesJSpecies", "IsolateI", "IsolateJ")
colsToIgnore <- which(names(geneticVsEpi) %in% colNamesToIgnore)

# Examine the levels of missing data
plotEpidemiologicalMetricDistributionsWithMissingData(
  geneticVsEpi[, -c(1, colsToIgnore)], fullNames, selection)

par(mfrow=c(1,1))

#### Fit Random Forest model ####

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
predictedGeneticDistancesRF <- predict(infoRF, geneticVsEpi[-trainRows, -c(1, colsToIgnore)])
corr <- cor(geneticVsEpi[-trainRows, "GeneticDistance"], predictedGeneticDistancesRF)
plotPredictedVersusActual(actual=geneticVsEpi[-trainRows, "GeneticDistance"],
                          predicted=predictedGeneticDistancesRF, 
                          main="Predicted vs. actual genetic distances (RF)")

#### Fit the Boosted Regression model ####

# Run Model
stepOutput <- gbm.step(data = geneticVsEpi[trainRows, -colsToIgnore],
                       gbm.x = colnames(geneticVsEpi[trainRows, -colsToIgnore])[-1], gbm.y = "GeneticDistance",
                       family = "poisson", tree.complexity = 5,
                       learning.rate = 0.01, bag.fraction = 0.5,
                       max.trees = 100000)

# Use the Model to predict
predictedGeneticDistancesBR <- predict(stepOutput,
                                       newdata = geneticVsEpi[-trainRows, -colsToIgnore],
                                       n.trees = stepOutput$n.trees,
                                       type = "response")

# Examine trained model prediction
corr <- cor(geneticVsEpi[-trainRows, "GeneticDistance"], predictedGeneticDistancesBR)
plotPredictedVersusActual(actual=geneticVsEpi[-trainRows, "GeneticDistance"],
                          predicted=predictedGeneticDistancesBR, 
                          main="Predicted vs. actual genetic distances (BR)")

#### Examine the extent of the correlated variables ####

# Examine correlation between epidemiological metrics
correlationTable <- calculateCorrelationBetweenEpiMetrics(
  geneticVsEpi[, -c(1, colsToIgnore)], fullNames)

# Find clusters of highly correlated metrics
threshold <- 0.55
clusters <- noteClustersOfMetrics(correlationTable, threshold, "black")
clusterSizes <- getClusterSizes(clusters)

#### Investigate the effect of metric removal on R squared and metric rankings ####

# Remove based upon correlation clusters
variableImportance <- list()
epiMetricImportance <- noteVariableImportance(importance=infoRF$importance,
                                              variableList=variableImportance,
                                              nMetricsPerClusterToRemove=0,
                                              colToUse=colToUse)
runRandomForestAnalysesIncrementallyRemovingMetricsFromCorrelationClusters(
  rSq=rSq, corr=corr, clusterSizes=clusterSizes, clusters=clusters, 
  colsToIgnore=colsToIgnore, epiMetricImportance=epiMetricImportance, 
  trainRows=trainRows, geneticVsEpi=geneticVsEpi, colToUse=colToUse)

# Remove based upon informativeness
variableImportance <- list()
epiMetricImportance <- noteVariableImportance(importance=infoRF$importance,
                                              variableList=variableImportance,
                                              nMetricsPerClusterToRemove=0,
                                              colToUse=colToUse)
runRandomForestAnalysesIncrementallyRemovingMetricsBasedOnInformativeness(
  rSq=rSq, corr=corr, colsToIgnore=colsToIgnore, 
  importance=as.vector(infoRF$importance[, colToUse]),
  trainRows=trainRows, colToUse=colToUse, geneticVsEpi=geneticVsEpi, 
  epiMetricImportance=epiMetricImportance)


# Remove based upon proportion missing data
variableImportance <- list()
epiMetricImportance <- noteVariableImportance(importance=infoRF$importance,
                                              variableList=variableImportance,
                                              nMetricsPerClusterToRemove=0,
                                              colToUse=colToUse)
runRandomForestAnalysesIncrementallyRemovingMetricsWithMissingData(
  rSq=rSq, corr=corr, colsToIgnore=colsToIgnore,
  importance=as.vector(infoRF$importance[, colToUse]),  
  trainRows=trainRows, colToUse=colToUse, geneticVsEpi=geneticVsEpi,
  epiMetricImportance=epiMetricImportance)

#### Examine the metric importance in random forest model ####

plotVariableImportance(infoRF=infoRF, colToUseForRF=colToUse, stepOutput=stepOutput,
                       fullNames=fullNames, nameColours=nameColours,
                       temporalCol=temporalCol, spatialCol=spatialCol,
                       networkCol=networkCol, showAxes=TRUE, selection=selection,
                       geneticVsEpi=geneticVsEpi, trainRows=trainRows, colsToIgnore=colsToIgnore)


##### Close PDF
dev.off()

# Save a copy of the R data generated - so I can come back and remake some plots if necessary
save.image(file=paste0(path, "FittingRF_", selection, "_", date, ".RData"))

#### Create importance figures for each analysis (BB, CC, CB) ####

# Load libraries
library(randomForest)
library(gplots)
library(gbm)
library(dismo)
library(codetools) # To fix weird error when loading RData file
library(basicPlotteR) # spread points on boxplot
library(stringr) # Counting string occurences in string

# Remove the infoRF and stepOutput objects
remove(list= c("infoRF", "stepOutput", "trainRows", "geneticVsEpi", "colsToIgnore", "predictedGeneticDistancesRF",
               "predictedGeneticDistancesBR"))

# Set the path
path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/GeneticVsEpidemiologicalDistances/"

# Note the column to use for the importance measure
colToUseForRF <- "%IncMSE"

# Assign each metric a colour
temporalCol="darkgoldenrod4"
spatialCol="red"
networkCol="blue"
nameColours <- assignMetricColours(temporalCol=temporalCol, spatialCol=spatialCol,
                                   networkCol=networkCol)

# Set the variable importance threshold for viewing the trends
importanceThresholdRF <- 0.1

# Note the analyses types
selections <- c("BB", "CC", "CB")

# Set the date analyses were completed on
dates <- c("29-08-19", "29-08-19", "29-08-19")

# Examine the analysis for each selection
for(i in seq_along(selections)){

  # Get the analysis info
  selection <- selections[i]
  date <- dates[i]

  # Attach the RData file and retrieve the infoRF object
  attach(paste0(path, "FittingRF_", selection, "_", date, ".RData"))
  infoRF <- infoRF
  stepOutput <- stepOutput
  trainRows <- trainRows
  geneticVsEpi <- geneticVsEpi
  colsToIgnore <- colsToIgnore
  predictedGeneticDistancesRF <- predictedGeneticDistancesRF
  predictedGeneticDistancesBR <- predictedGeneticDistancesBR
  detach(paste0("file:", path, "FittingRF_", selection, "_", date, ".RData"), character.only=TRUE)

  # Note the full names
  fullNames <- noteFullNames(selection)

  # Open a pdf
  file <- paste0(path, "VariableImportance_", selection, "_", date, ".pdf")
  pdf(file, height=10, width=10)

  # Plot the variable importance
  variableImportance <- 
    plotVariableImportance(infoRF=infoRF, colToUseForRF=colToUseForRF, stepOutput=stepOutput,
                           fullNames=fullNames, nameColours=nameColours,
                           temporalCol=temporalCol, spatialCol=spatialCol,
                           networkCol=networkCol, showAxes=TRUE, selection=selection,
                           geneticVsEpi=geneticVsEpi, trainRows=trainRows, colsToIgnore=colsToIgnore)

  # Plot the predictor variable trends
  plotBRPredictorVariableTrends(stepOutput, geneticVsEpi, fullNames, variableImportance, labelCex=0.5)

  # Get and set the margins
  currentMar <- par()$mar
  par(mar=c(1.5,0.5,0.5,0.5))
  
  # Store the predictor data
  predictorData <- geneticVsEpi[trainRows, -c(1, colsToIgnore)]
  
  # Get the predictor variable names
  predictors <- colnames(predictorData)
  
  # Remove predictors with an importance below threshold
  predictors <- predictors[which(variableImportance[predictors, "RandomForestImportance"] > importanceThresholdRF)]
  
  # Create a plotting grid
  nPlotsPerRow <- ceiling(sqrt(length(predictors)))
  par(mfrow=c(nPlotsPerRow, nPlotsPerRow))
  
  # Examine each predictor
  for(i in seq_along(predictors)){
    
    # Get the trend data
    trendData <- partialPlot(infoRF, predictorData, predictors[i], plot=FALSE)
    trendData <- as.data.frame(trendData)
    
    # Remove any "-1" rows
    trendData <- trendData[which(trendData[, 1] != -1), ]
    
    # Get the raw predictor values and remove -1s
    predictorValues <- geneticVsEpi[, predictors[i]]
    responseValues <- geneticVsEpi$GeneticDistance
    keep <- predictorValues != -1
    predictorValues <- predictorValues[keep]
    responseValues <- responseValues[keep]
    
    # Plot the raw data and overlay the trend
    plotRawDataAndTrend(predictorValues, responseValues, trendData, rgb(0,0,0, 0.1))
    
    # Add predictor name to X axis label
    addPredictorNameToPlot(fullNames[[predictors[i]]], labelCex=0.5)
  }
  
  # Reset plotting parameters
  par(mar=currentMar, mfrow=c(1,1))
  
  # Plot the correlation between predicted and actual values for random forest and boosted regression
  plotPredictedVersusActual(actual=geneticVsEpi[-trainRows, "GeneticDistance"],
                            predicted=predictedGeneticDistancesRF,
                            main="Predicted vs. actual genetic distances (RF)")
  plotPredictedVersusActual(actual=geneticVsEpi[-trainRows, "GeneticDistance"],
                            predicted=predictedGeneticDistancesBR,
                            main="Predicted vs. actual genetic distances (BR)")

  # Remove the infoRF and stepOutput objects
  remove(list= c("infoRF", "stepOutput", "trainRows", "geneticVsEpi", "colsToIgnore", "predictedGeneticDistancesRF",
                 "predictedGeneticDistancesBR"))

  # Close the PDF
  dev.off()
}

#### FUNCTIONS ####

plotRawDataAndTrend <- function(predictorValues, responseValues, trendData, col){
  
  # Check if the predictor variable is factor
  if(is.factor(predictorValues)){
    
    # Check that all the levels are in response factor
    predictorValues <- droplevels(predictorValues)
    trendData[, 1] <- droplevels(trendData[, 1])
                                          
    # Plot the raw data
    plot(x=predictorValues, y=responseValues, pch=19, bty="n", ylab="", xlab="", yaxt="n", xaxt="n",
         frame=FALSE, border=rgb(0,0,0,0))
    
    # Overlay the raw points
    data <- data.frame("X"=predictorValues, "Y"=responseValues)
    spreadPointsMultiple(data=data, responseColumn="Y", categoriesColumn="X", col=col)
    
    # Overlay the partial dependence relationship
    points(x=trendData[, 1] , y=trendData[, 2], type="l", lwd=3, col="red")
    
  }else{
    
    # Plot the raw data
    plot(x=predictorValues, y=responseValues, pch=19, bty="n", ylab="", xlab="", yaxt="n", xaxt="n", col=col)
    
    # Overlay the partial dependence relationship
    points(x=trendData[, 1], y=trendData[, 2], type="l", lwd=3, col="red")
  }
}

addPredictorNameToPlot <- function(fullName, pad=0.05, labelCex){
 
  # Get the axis limits of the current plot
  axisLimits <- par()$usr
  plotWidth <- axisLimits[2] - axisLimits[1]
  plotWidth <- plotWidth - (pad*plotWidth)

  # Calculate number of characters that can fit in plotting window
  nCharInPlot <- floor(plotWidth / strwidth("X", cex=labelCex))
    
  # Add axes labels
  if(nchar(fullName) > nCharInPlot){
    
    # Split the predictor variable name over multiple lines
    fullName <- splitPredictorVariableNameOverMultipleLines(fullName, nCharInPlot)
    
    mtext(fullName, side=1, line=str_count("X\nY\n", "\n")/2, cex=labelCex)
  }else{
    mtext(fullName, side=1, line=0, cex=labelCex)
  }
}

plotBRPredictorVariableTrends <- function(stepOutput, geneticVsEpi, fullNames, variableImportance, importanceThresholdRF=0.5,
                                          labelCex=1, col=rgb(0,0,0, 0.1)){
  
  # Get the predictor variable names
  predictors <- stepOutput$var.names
  
  # Remove predictors with an importance below threshold
  predictors <- predictors[which(variableImportance[predictors, "RandomForestImportance"] > importanceThresholdRF)]
  
  # Get and set the margins
  currentMar <- par()$mar
  par(mar=c(2,0.5,0.5,0.5))
  
  # Create a plotting grid
  nPlotsPerRow <- ceiling(sqrt(length(predictors)))
  par(mfrow=c(nPlotsPerRow, nPlotsPerRow))
  
  for(predictor in predictors){
    
    # Get the trend data
    trendData <- plot(stepOutput, i.var=predictor, return.grid=TRUE)
    
    # Remove any "-1" rows
    trendData <- trendData[which(trendData[, 1] != -1), ]
    
    # Get the raw predictor values and remove -1s
    predictorValues <- geneticVsEpi[, predictor]
    responseValues <- geneticVsEpi$GeneticDistance
    keep <- predictorValues != -1
    predictorValues <- predictorValues[keep]
    responseValues <- responseValues[keep]
    
    # Plot the raw data and overlay the trend
    plotRawDataAndTrend(predictorValues, responseValues, trendData, col)
    
    # Add predictor name to X axis label
    addPredictorNameToPlot(fullNames[[predictor]], labelCex=labelCex)
  }
  
  # Reset plotting parameters
  par(mar=currentMar, mfrow=c(1,1))
}

splitPredictorVariableNameOverMultipleLines <- function(fullName, nCharInPlot){
  
  # Find the first split point
  splitPoint <- nCharInPlot
  for(characterIndex in seq(from=nCharInPlot, to=1, by=-1)){
    
    # Check if found space
    if(substr(fullName, characterIndex, characterIndex) == " "){
      splitPoint <- characterIndex
      break
    }
  }
  
  # Extract the strings before and after split
  first <- substr(fullName, 1, splitPoint)
  second <- substr(fullName, splitPoint+1, nchar(fullName))
  
  # Check if second part still more characters than width
  if(nchar(second) > nCharInPlot){
    second <- splitPredictorVariableNameOverMultipleLines(second, nCharInPlot)
  }
  
  return(paste0(first, "\n", second))
}

removeMetricWithMostMissingData <- function(geneticVsEpi, propMissingData){
  
  # Get the order the proportions
  order <- order(propMissingData)
  
  # Remove the last index in that order
  colToRemove <- colnames(geneticVsEpi)[-1][order[length(order)]]
  geneticVsEpi <- geneticVsEpi[ ,-which(colnames(geneticVsEpi) == colToRemove)]
  
  infoForRFModel <- list()
  infoForRFModel[["geneticVsEpiTable"]] <- geneticVsEpi
  infoForRFModel[["Removed"]] <- colToRemove
  
  return(infoForRFModel)
}

calculateProportionMissingData <- function(geneticVsEpi){
  
  propMissingData <- rep(0, ncol(geneticVsEpi))
  for(row in 1:nrow(geneticVsEpi)){
    
    for(col in 2:ncol(geneticVsEpi)){
      
      if(geneticVsEpi[row, col] == -1){
        propMissingData[col-1] <- propMissingData[col-1] + 1
      }
    }
  }
  
  propMissingData <- propMissingData / nrow(geneticVsEpi)
  
  return(propMissingData)
}

runRandomForestAnalysesIncrementallyRemovingMetricsWithMissingData <- function(
  rSq, corr, colsToIgnore, importance, trainRows, colToUse, geneticVsEpi, epiMetricImportance){
  
  # Initialise an output table and insert initial values from full random forest model
  par(mfrow=c(1,2))
  outputTable <- data.frame("NumberMetricsRemoved" = c(0),
                            "MetricRemoved" = c(""),
                            "PseudoRSquared" = c(rSq),
                            "Correlation" = c(corr), 
                            stringsAsFactors=FALSE)
  
  # Store the original rf model output
  infoForRFModel <- list("geneticVsEpiTable" = geneticVsEpi[, -colsToIgnore])
  
  # Initialise a value for tha max importance across the metrics
  maxImportance <- 0
  
  # Count number of columns with missing data
  propMissingData <- calculateProportionMissingData(geneticVsEpi[, -colsToIgnore])
  nColsWithMissingData <- length(which(propMissingData > 0))
  
  # Begin fitting random forest models with increasingly fewer metrics
  for(numberMetricsToRemove in 1:nColsWithMissingData){
    
    # Remove next metric
    infoForRFModel <- removeMetricWithMostMissingData(
      infoForRFModel[["geneticVsEpiTable"]], propMissingData)
    
    # Find the optimal mtry value
    optimalMtry <- findOptimalMtry(
      response=infoForRFModel[["geneticVsEpiTable"]][trainRows, "GeneticDistance"],
      predictors=infoForRFModel[["geneticVsEpiTable"]][trainRows, -1],
      mTryInitial=ifelse(ncol(infoForRFModel[["geneticVsEpiTable"]]) > 3, 3, 
                         ncol(infoForRFModel[["geneticVsEpiTable"]]) - 1),
      nTrees=500, plot=TRUE)
    
    # Fit the Random Forest model
    rfModel <- randomForest(
      infoForRFModel[["geneticVsEpiTable"]][trainRows, "GeneticDistance"] ~ .,
      data=infoForRFModel[["geneticVsEpiTable"]][trainRows, -1],
      mtry=optimalMtry, importance=TRUE, ntree=1000,
      keep.forest=TRUE, norm.votes=FALSE, proximity=FALSE,
      do.trace=FALSE)
    rSq <- rfModel$rsq[length(rfModel$rsq)]
    plot(rfModel, las=1)
    
    # Check the max importance value
    if(max(rfModel$importance[, colToUse]) > maxImportance){
      maxImportance <- max(rfModel$importance[, colToUse])
    }
    
    # Note the rankings of the variables
    epiMetricImportance <- noteVariableImportance(importance=rfModel$importance,
                                                  variableList=epiMetricImportance,
                                                  numberMetricsToRemove,
                                                  colToUse=colToUse)
    
    # Test the Random Forest model
    predictions <- predict(rfModel,
                           infoForRFModel[["geneticVsEpiTable"]][-trainRows, -1])
    corr <- cor(infoForRFModel[["geneticVsEpiTable"]][-trainRows, "GeneticDistance"],
                predictions)
    
    # Store the results
    outputTable[numberMetricsToRemove + 1, "NumberMetricsRemoved"] <- numberMetricsToRemove
    outputTable[numberMetricsToRemove + 1, "MetricRemoved"] <- infoForRFModel[["Removed"]]
    outputTable[numberMetricsToRemove + 1, "PseudoRSquared"] <- rSq
    outputTable[numberMetricsToRemove + 1, "Correlation"] <- corr
    
    # Calculate the proportion of missing data
    propMissingData <- calculateProportionMissingData(
      infoForRFModel[["geneticVsEpiTable"]])
    
    # Print progress information
    cat("###########################################################################\n")
    cat(paste("Finished fitting Random Forest model", "\n",
              "Removed ", numberMetricsToRemove," metrics\n", sep=""))
  }
  
  plot(x=outputTable$NumberMetricsRemoved, 
       y=outputTable$PseudoRSquared, type="o", las=1, cex.axis=0.9,
       xlab="Number of Metrics Removed",
       ylab="Proportion Variation Explained")
  
  par(mfrow=c(1,1))
  
  plotVariableImportanceAgainstNMetricsRemoved(
    epiMetricImportance, 
    maxNMetricsRemoved=numberMetricsToRemove,
    yLim=c(0, maxImportance))
}

runRandomForestAnalysesIncrementallyRemovingMetricsBasedOnInformativeness <- function(
  rSq, corr, colsToIgnore, importance, trainRows, colToUse, geneticVsEpi, epiMetricImportance){

  # Initialise an output table and insert initial values from full random forest model
  par(mfrow=c(1,2))
  outputTable <- data.frame("NumberMetricsRemoved" = c(0),
                            "MetricRemoved" = c(""),
                            "PseudoRSquared" = c(rSq),
                            "Correlation" = c(corr), 
                            stringsAsFactors=FALSE)
  
  # Store the original rf model output
  infoForRFModel <- list("geneticVsEpiTable" = geneticVsEpi[, -colsToIgnore])
  
  # Initialise a value for tha max importance across the metrics
  maxImportance <- 0
  
  # Begin fitting random forest models with increasingly fewer metrics
  # Note that will keep removing until two columns left
  # Random forest breaks with only one column
  nEpidemiologicalMetrics <- (ncol(geneticVsEpi) - length(colsToIgnore)) - 1
  for(numberMetricsToRemove in 1:(nEpidemiologicalMetrics - 2)){
    
    # Remove next metric
    infoForRFModel <- removeLeastInformativeMetric(infoForRFModel[["geneticVsEpiTable"]],
                                                   importance)
    
    # Find the optimal mtry value
    optimalMtry <- findOptimalMtry(
      response=infoForRFModel[["geneticVsEpiTable"]][trainRows, "GeneticDistance"],
      predictors=infoForRFModel[["geneticVsEpiTable"]][trainRows, -1],
      mTryInitial=ifelse(ncol(infoForRFModel[["geneticVsEpiTable"]]) > 3, 3, 
                         ncol(infoForRFModel[["geneticVsEpiTable"]]) - 1),
      nTrees=500, plot=TRUE)

    # Fit the Random Forest model
    rfModel <- randomForest(
      infoForRFModel[["geneticVsEpiTable"]][trainRows, "GeneticDistance"] ~ .,
      data=infoForRFModel[["geneticVsEpiTable"]][trainRows, -1],
      mtry=optimalMtry, importance=TRUE, ntree=1000,
      keep.forest=TRUE, norm.votes=FALSE, proximity=FALSE,
      do.trace=FALSE)
    rSq <- rfModel$rsq[length(rfModel$rsq)]
    plot(rfModel, las=1)
    
    # Note the variable importance
    importance <- as.vector(rfModel$importance[, colToUse])
    
    # Check the max value
    if(max(importance) > maxImportance){
      maxImportance <- max(importance)
    }
    
    # Note the rankings of the variables
    epiMetricImportance <- noteVariableImportance(importance=rfModel$importance,
                                                  variableList=epiMetricImportance,
                                                  numberMetricsToRemove,
                                                  colToUse=colToUse)
    
    # Test the Random Forest model
    predictions <- predict(rfModel,
                           infoForRFModel[["geneticVsEpiTable"]][-trainRows, -1])
    corr <- cor(infoForRFModel[["geneticVsEpiTable"]][-trainRows, "GeneticDistance"],
                predictions)
    
    # Store the results
    outputTable[numberMetricsToRemove + 1, "NumberMetricsRemoved"] <- numberMetricsToRemove
    outputTable[numberMetricsToRemove + 1, "MetricRemoved"] <- infoForRFModel[["Removed"]]
    outputTable[numberMetricsToRemove + 1, "PseudoRSquared"] <- rSq
    outputTable[numberMetricsToRemove + 1, "Correlation"] <- corr

    # Print progress information
    cat("###########################################################################\n")
    cat(paste("Finished fitting Random Forest model", "\n",
              "Removed ", numberMetricsToRemove," metrics\n", sep=""))
  }
  
  plot(x=outputTable$NumberMetricsRemoved, 
       y=outputTable$PseudoRSquared, type="o", las=1, cex.axis=0.9,
       xlab="Number of Metrics Removed",
       ylab="Proportion Variation Explained")
 
  par(mfrow=c(1,1))
  
  plotVariableImportanceAgainstNMetricsRemoved(
    epiMetricImportance, 
    maxNMetricsRemoved=numberMetricsToRemove,
    yLim=c(0, maxImportance))
}

removeLeastInformativeMetric <- function(geneticVsEpi, importance){
  
  # Get the order the importance values
  order <- order(importance)
  
  # Remove the last index in that order
  colToRemove <- colnames(geneticVsEpi)[-1][order[length(order)]]
  geneticVsEpi <- geneticVsEpi[ ,-which(colnames(geneticVsEpi) == colToRemove)]
  
  infoForRFModel <- list()
  infoForRFModel[["geneticVsEpiTable"]] <- geneticVsEpi
  infoForRFModel[["Removed"]] <- colToRemove
  
  return(infoForRFModel)
}

plotVariableImportanceAsScatterPlot <- function(variableImportance, temporalCol, spatialCol, networkCol,
                                                fullNames){
  
  # Set the plotting margins
  currentMar <- par()$mar
  par(mar=c(4.1, 4.1, 4.1, 0.1))
  
  # Plot the variable importance of random forest model against boosted regression
  plot(x=variableImportance$RandomForestImportance, y=variableImportance$BoostedRegressionImportance, 
       ylab="Relative Influence (Boosted Regresssion)", xlab="% Increase MSE (Random Forest)",
       main="Variable importance in regression models",
       pch=19, bty="n", col=rgb(0,0,0, 0.5), las=1, cex=1.5)
  
  # Get the full names of the predictor variables
  variableNames <- getFullVariableNames(rownames(variableImportance), fullNames)
  
  # Get the colours associated with each predictor
  variableColours <- getVariableColours(rownames(variableImportance),
                                        nameColours=nameColours)
  
  # Add the labels for each point
  addTextLabels(xCoords=variableImportance$RandomForestImportance, 
                yCoords=variableImportance$BoostedRegressionImportance,
                labels=variableNames, avoidPoints=TRUE, cex.label=0.6, keepLabelsInside=TRUE,
                col.label=variableColours, lty=3, col.line=rgb(0,0,0, 0.75),
                col.background=rgb(0,0,0, 0.2))
  
  # Add a legend
  legend("topleft", legend=c("Temporal", "Spatial", "Network"), 
         text.col=c(temporalCol, spatialCol, networkCol),
         bty="n")
  
  # Reset the plotting margins
  par(mar=currentMar)
}

plotVariableImportance <- function(infoRF, colToUseForRF, stepOutput, fullNames, nameColours,
                                   temporalCol, spatialCol, networkCol, showAxes, selection,
                                   geneticVsEpi, trainRows, colsToIgnore){
  
  # Get the current margin settings
  currentMar <- par()$mar
  
  # Get the variable importance from the RF model
  variableImportanceRF <- as.data.frame(infoRF$importance)
  
  # Get the variable importance from the BR model
  variableImportanceBR <- summary.gbm(stepOutput, plotit=FALSE, method=permutation.test.gbm)
  rownames(variableImportanceBR) <- variableImportanceBR$var
  
  # Combine the importance tables into one
  variableImportance <- data.frame(RandomForestImportance=variableImportanceRF[, colToUseForRF], 
                                   BoostedRegressionImportance=variableImportanceBR[rownames(variableImportanceRF), "rel.inf"],
                                   stringsAsFactors=FALSE)
  rownames(variableImportance) <- rownames(variableImportanceRF)

  # Order the table by importance - RF first, then BR
  variableImportance <- variableImportance[order(variableImportance$RandomForestImportance, 
                                                 variableImportance$BoostedRegressionImportance,
                                                 decreasing=FALSE), ]
    
  # Normalise the Random Forest and Boosted Regression importance values to vary between zero and 1
  variableImportance$RandomForestImportanceNorm <- variableImportance$RandomForestImportance / 
    max(variableImportance$RandomForestImportance)
  variableImportance$BoostedRegressionImportanceNorm <- variableImportance$BoostedRegressionImportance / 
    max(variableImportance$BoostedRegressionImportance)
  
  # Transpose the table
  transpose <- as.matrix(t(variableImportance))
  
  # Get full Variable Names
  variableNames <- getFullVariableNames(rownames(variableImportance), fullNames)
  
  # Get Variable Colours - double up colours (1 for RF and 1 for BR) for side by side bar plot
  variableColoursDoubled <- getVariableColours(doubleUpRowNames(rownames(variableImportance)),
                                               nameColours=nameColours)
  variableColours <- getVariableColours(rownames(variableImportance),
                                        nameColours=nameColours)
  
  # Create bar plot illustrating the relative variable importance from RF and BR
  par(mfrow=c(1,1))

  # Set the margin sizes according to the comparison being made
  marginSizes <- list(
    "BB" = 37,
    "CC" = 42.5,
    "CB" = 30
  )
  par(mar=c(0,marginSizes[[selection]],2,1)) # bottom, left, top, right
  
  # Check if we need extra space for the axes
  if(showAxes == TRUE){
    par(mar=c(2,marginSizes[[selection]],3,1)) # bottom, left, top, right
  }
  
  # Create the bar chart for the RF and BR values - normalised values plotted side by side
  plot <- barplot(transpose[c("RandomForestImportanceNorm", "BoostedRegressionImportanceNorm"),],
                  horiz=TRUE, beside=TRUE, xaxt="n", col=variableColoursDoubled, main="", col.axis=rgb(0,0,0,0),
                  density=c(25, 75), angle=c(45, 90))
  
  # Add a legend to illustrate which bars are the RF and BR ones
  axisLimits <- par()$usr
  legendInfo <- list(
    "BB" = c(0.25, 0.8),
    "CC" = c(0.1, 0.7),
    "CB" = c(0.3, 0.8)
  )
  legend(x=legendInfo[[selection]][1], y=axisLimits[3] + (axisLimits[4] - axisLimits[3])/5,
         legend=c(paste0("corr = ", round(cor(variableImportance$RandomForestImportance, 
                                              variableImportance$BoostedRegressionImportance), digits=2)),
                  "Boosted Regression", "Random Forest", "Temporal", "Spatial", "Network"),
         density=c(0, 75, 25, 0, 0, 0), angle=c(90, 90, 45), border=c("white", "black", "black", "white", "white", "white"), 
         bty="n", text.col=c("black", "black", "black", temporalCol, spatialCol, networkCol), cex=legendInfo[[selection]][2],
         xpd=TRUE)
  
  # Add labels for each bar
  barPositions <- plot[1,] + (plot[2, ] - plot[1, ])/2
  xLabPosition <- 0
  text(labels=variableNames, 
       col=variableColours,
       x=rep(xLabPosition,nrow(variableImportance)),
       y=barPositions,
       srt=0, pos=2, xpd=TRUE, cex=1.1)

  # If asked add an axis for the BR and the RF values
  if(showAxes == TRUE){
    
    # Add an axis for the Random Forest analysis - to the bottom
    # mgp(axisLabelPosition, tickLabelPosition, tickPosition)
    range <- range(variableImportance$RandomForestImportance)
    at <- seq(from=range[1], to=range[2], by=(range[2] - range[1])/5)
    axis(side=1, line=-0.5, cex.axis=0.75, mgp=c(0, .4, 0), at=at/range[2], labels=round(at, digits=0))
    mtext("% Increase MSE (RF)", side=1, line=0.75, cex=0.75)
    
    # Add an axis for the Boosted Regression analysis - to the top
    range <- range(variableImportance$BoostedRegressionImportance)
    at <- seq(from=range[1], to=range[2], by=(range[2] - range[1])/5)
    axis(side=3, line=-0.5, cex.axis=0.75, mgp=c(0, .4, 0), at=at/range[2], labels=round(at, digits=0))
    mtext("Relative Influence (BR)", side=3, line=0.75, cex=0.75)
  }
  
  ## Add an inset graph comparing the RF and BR values
  # addLinearComparisonInset(axisLimits, variableImportance, originX=0.3, originY=20,
  #                          xAxisProp=0.75, yAxisProp=0.2, padProp=0.05, pch=19,
  #                          col=rgb(0,0,0, 0.5), labCex=0.65)
  
  # Reset the margin sizes
  par(mar=currentMar)
  
  # Plot the variable importance as a scatter plot
  plotVariableImportanceAsScatterPlot(variableImportance, temporalCol, spatialCol, networkCol,
                                      fullNames)
  
  return(variableImportance)
}

addLinearComparisonInset <- function(axisLimits, variableImportance, 
                                     originX, originY, 
                                     xAxisProp=0.5, yAxisProp=0.5,
                                     padProp=0.05, labCex=0.75, ...){
  
  # Note the axis lengths
  xLength <- (axisLimits[2] - axisLimits[1]) * xAxisProp
  yLength <- (axisLimits[4] - axisLimits[3]) * yAxisProp
  
  # Create a frame around the inset plot
  rect(xleft=originX, ybottom=originY, xright=originX+xLength, ytop=originY+yLength,
       xpd=TRUE)
  
  # Adjust the X and Y lengths to pad the axes
  padX <- padProp * xLength
  padY <- padProp * yLength
  xLengthPadded <- xLength - (2*padX)
  yLengthPadded <- yLength - (2*padY)
  
  # Note the range of X and Y axis
  yRange <- range(variableImportance$BoostedRegressionImportance)
  xRange <- range(variableImportance$RandomForestImportance)
  
  # Calculate a unit on the X and Y axes
  xUnit <- xLengthPadded / (xRange[2] - xRange[1])
  yUnit <- yLengthPadded / (yRange[2] - yRange[1])
  
  # Add the points
  points(x=originX + padX + (variableImportance$RandomForestImportance * xUnit),
         y=originY + padY + (variableImportance$BoostedRegressionImportance * yUnit),
         xpd=TRUE, ...)
  
  # Add X axis label
  text(x=originX + (0.5*xLength), y=originY - padY,
       labels="RF variable importance", cex=labCex)
  
  # Add Y axis label
  text(x=originX - padX, y=originY + (0.5*yLength),
       labels="BR variable importance", cex=labCex, srt=90)
}

doubleUpRowNames <- function(rowNames){
  
  # Create a vector to store each element twice
  output <- c()
  
  # Examine each row name
  for(name in rowNames){
    
    # Store the current row name twice
    output <- c(output, name, name)
  }
  
  return(output)
}

runRandomForestAnalysesIncrementallyRemovingMetricsFromCorrelationClusters <- function(
  rSq, corr, clusterSizes, clusters, colsToIgnore, epiMetricImportance, trainRows,
  geneticVsEpi, colToUse){
  
  stop <- FALSE
  nMetricsPerClusterToRemove <- 1
  par(mfrow=c(1,2))
  
  outputTable <- data.frame("NumberMetricsRemovedPerCluster" = c(0),
                            "MetricsRemoved" = c(""),
                            "NumberMetricsRemoved" = c(0),
                            "PseudoRSquared" = c(rSq),
                            "Correlation" = c(corr), 
                            "ClusterSizes" = paste(clusterSizes, collapse=","), 
                            stringsAsFactors=FALSE)
  
  # Create a list to store the variable importance after variable removal
  maxImportance <- 0
  
  while(stop == FALSE){
    
    # Remove x metrics from clusters and geneticVsEpi table
    infoForRFModel <- removeVariablesByImportanceFromClustersAndTable(
      clusters=clusters,
      geneticVsEpi=geneticVsEpi[, -colsToIgnore],
      variableImportance=epiMetricImportance, 
      nToRemove=nMetricsPerClusterToRemove,
      useOriginalRfModelImportance=FALSE)
    
    # Find the optimal mtry value
    optimalMtry <- findOptimalMtry(
      response=infoForRFModel[["geneticVsEpiTable"]][trainRows, "GeneticDistance"],
      predictors=infoForRFModel[["geneticVsEpiTable"]][trainRows, -1],
      mTryInitial=3, nTrees=500, plot=TRUE)
    abline(v=optimalMtry, col="red", lty=2)
    
    # Fit the Random Forest model
    rfModel <- randomForest(
      infoForRFModel[["geneticVsEpiTable"]][trainRows, "GeneticDistance"] ~ .,
      data=infoForRFModel[["geneticVsEpiTable"]][trainRows, -1],
      mtry=optimalMtry, importance=TRUE, ntree=1000,
      keep.forest=TRUE, norm.votes=FALSE, proximity=FALSE,
      do.trace=FALSE)
    rSq <- rfModel$rsq[length(rfModel$rsq)]
    plot(rfModel, las=1)
    
    # Note the variable importance
    epiMetricImportance <- noteVariableImportance(importance=rfModel$importance,
                                                  variableList=epiMetricImportance,
                                                  nMetricsPerClusterToRemove=nMetricsPerClusterToRemove,
                                                  colToUse=colToUse)
    
    if(max(rfModel$importance[, "%IncMSE"]) > maxImportance){
      maxImportance <- max(rfModel$importance[, "%IncMSE"])
    }
    
    
    # Test the Random Forest model
    predictions <- predict(rfModel,
                           infoForRFModel[["geneticVsEpiTable"]][-trainRows, -1])
    corr <- cor(infoForRFModel[["geneticVsEpiTable"]][-trainRows, "GeneticDistance"],
                predictions)
    
    # Examine the size of the clusters following metric removal
    clusterSizes <- getClusterSizes(infoForRFModel[["clusters"]])
    if(max(clusterSizes) == 1){
      stop = TRUE
    }
    
    # Store the results
    outputTable[nMetricsPerClusterToRemove + 1, "NumberMetricsRemovedPerCluster"] <- nMetricsPerClusterToRemove
    outputTable[nMetricsPerClusterToRemove + 1, "MetricsRemoved"] <- paste(infoForRFModel[["removed"]], collapse=",")
    outputTable[nMetricsPerClusterToRemove + 1, "NumberMetricsRemoved"] <- length(infoForRFModel[["removed"]])
    outputTable[nMetricsPerClusterToRemove + 1, "PseudoRSquared"] <- rSq
    outputTable[nMetricsPerClusterToRemove + 1, "Correlation"] <- corr
    outputTable[nMetricsPerClusterToRemove + 1, "ClusterSizes"] <- paste(clusterSizes, collapse=",")
    
    # Print progress information
    cat("###########################################################################\n")
    cat(paste("Finished fitting Random Forest model", "\n",
              "Removed ", length(infoForRFModel[["removed"]])," metrics (",
              nMetricsPerClusterToRemove, " per cluster):\n",
              paste(infoForRFModel[["removed"]], collapse=","), "\nCorrelation = ",
              round(corr, digits=2), "\nRsq = ", round(rSq, digits=2), "\n", sep=""))

    # Increment the number of metrics to remove
    nMetricsPerClusterToRemove <- nMetricsPerClusterToRemove + 1
  }
  
  plot(x=outputTable$NumberMetricsRemovedPerCluster, 
       y=outputTable$PseudoRSquared, type="o", las=1, cex.axis=0.9,
       xlab="Number of Metrics Removed per Cluster",
       ylab="Proportion Variation Explained")
  
  plot(x=outputTable$NumberMetricsRemoved, 
       y=outputTable$PseudoRSquared, type="o", las=1, cex.axis=0.9,
       xlab="Number of Metrics Removed",
       ylab="Proportion Variation Explained")
  
  par(mfrow=c(1,1))
  plotVariableImportanceAgainstNMetricsRemovedPerCluster(epiMetricImportance, 
                                                         nMetricsPerClusterToRemove - 1,
                                                         c(0, maxImportance))
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

selectAppropriateComparisonsForSelection <- function(selection, geneticVsEpi){
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
    "HostRelatedness" = "black",
    "RandomUniform" = "black",
    "RandomBoolean" = "black"
  )
  
  return(nameColours)
}

noteFullNames <- function(selection){
  
  # Define the fullnames for the epidemiological metrics
  fullNames <- list(
    "SameMainGroup" = "Same main group?",
    "SameSampledGroup" = "Same sampled group?",                           
    "SameInfectedGroup" = "Same infected group?",
    "PeriodSpentAliveTogether" = "Number of days overlap between the recorded lifespans",                   
    "PeriodSpentInfectedTogether" = "Number of days overlap between the infected lifespans",
    "PeriodSpentInSameGroup" = "Number of days spent in same group",                     
    "TimeBetweenInfectionDetection" = "Number of days between infection detection dates",
    "TimeBetweenSampling" = "Number of days between sampling dates",                        
    "TimeBetweenBreakdown" = "Number of days between breakdown dates",
    "DistanceBetweenMainGroups" = "Spatial distance between main groups",
    "DistanceBetweenSampledGroups" = "Spatial distance between sampled groups",               
    "DistanceBetweenInfectedGroups" = "Spatial distance between infected groups",
    "NMovementsBetweenMainGroups" = "Number of recorded animal movements between main groups",                
    "NMovementsBetweenSampledGroups" = "Number of recorded animal movements between sampled groups",
    "NMovementsBetweenInfectedGroups" = "Number of recorded animal movements between infected groups",            
    "SameAnimal" = "Isolates from same animal?",
    "ShortestPathLengthMain" = "Shortest path length between main groups",                     
    "MeanNMovementsOnEdgesOfShortestPathMain" = "Mean number of animals traversing edges of shortest path between main groups",
    "ShortestPathLengthSampled" = "Shortest path length between sampled groups",                  
    "MeanNMovementsOnEdgesOfShortestPathSampled" = "Mean number of animals traversing edges of shortest path between sampled groups",
    "ShortestPathLengthInfected" = "Shortest path length between infected groups",                 
    "MeanNMovementsOnEdgesOfShortestPathInfected" = "Mean number of animals traversing edges of shortest path between infected groups",
    "NSharedAnimalsBetweenMainGroups" = "Number of animals recorded in both main groups",            
    "NSharedAnimalsBetweenSampledGroups" = "Number of animals recorded in both sampled groups",
    "NSharedAnimalsBetweenInfectedGroups" = "Number of animals recorded in both infected groups",
    "ShortestPathLengthEXCLMain" = "Shortest path length between main groups (some groups excluded)",                     
    "MeanNMovementsOnEdgesOfShortestPathEXCLMain" = "Mean number of animals traversing edges of shortest path between main groups (some groups excluded)",   
    "ShortestPathLengthEXCLSampled" = "Shortest path length between sampled groups (some groups excluded)",               
    "MeanNMovementsOnEdgesOfShortestPathEXCLSampled" = "Mean number of animals traversing edges of shortest path between sampled groups (some groups excluded)",
    "ShortestPathLengthEXCLInfected" = "Shortest path length between infected groups (some groups excluded)",
    "MeanNMovementsOnEdgesOfShortestPathEXCLInfected" = "Mean number of animals traversing edges of shortest path between main groups (some groups excluded)",
    "CentroidDistBetweenMain" = "Distance from closest land parcel to main group using centroids",
    "CentroidDistBetweenSamp" = "Distance from closest land parcel to sampled group using centroids",
    "HostRelatedness" = "Genetic relatedness of animals",
    "RandomUniform" = "Random uniform sample",
    "RandomBoolean" = "Random boolean sample"
  )
  
  # Replace some words in full names depending on the selection
  #   * group with herd/social group
  #   * animals/animal with cattle/badgers
  #   * animal movements with dispersal events
  if(selection == "BB"){
    
    for(key in names(fullNames)){
      
      # Replace "group" with "social group"
      fullNames[[key]] <- gsub("group", "social group", fullNames[[key]])
      
      # Replace "animal movements" with "dispersal events"
      fullNames[[key]] <- gsub("animal movements", "dispersal events", fullNames[[key]])
      
      # Replace "animal" with "badger"
      fullNames[[key]] <- gsub("animal", "badger", fullNames[[key]])
    }
    
  }else if(selection == "CC"){
    
    for(key in names(fullNames)){
      
      # Skip cattle-badger metrics
      if(grepl(fullNames[[key]], pattern="land parcel")){
        next
      }
      
      # Replace "group" with "herd" - check if comparing cattle and badgers - not really necessary since this metric never considered
      fullNames[[key]] <- gsub("group", "herd", fullNames[[key]])
      
      # Replace "animal" with "cattle"
      fullNames[[key]] <- gsub("animal", "cattle", fullNames[[key]])
      
      # Replace "cattles" with "cattle"
      fullNames[[key]] <- gsub("cattles", "cattle", fullNames[[key]])
    }
    
  }else{
    
    for(key in names(fullNames)){
      
      # Treat land parcel metrics slightly differently
      if(grepl(fullNames[[key]], pattern="land parcel")){
        fullNames[[key]] <- gsub("group", "social group", fullNames[[key]])
        next
      }
      
      # Replace "group" with "herd/social group" - check if comparing cattle and badgers - not really necessary since this metric never considered
      fullNames[[key]] <- gsub("group", "herd/social group", fullNames[[key]])
    }
  }
  
  return(fullNames)
}

getVariableColours <- function(rowNames, nameColours){

  colours <- c()
  for(index in seq_along(rowNames)){
    colours[index] <- nameColours[[rowNames[index]]]
  }
  
  return(colours)
}

plotVariableImportanceAgainstNMetricsRemoved <- function(variableImportance, 
                                                                   maxNMetricsRemoved,
                                                                   yLim){
  
  plot(x=0:maxNMetricsRemoved, y=rep(-1, maxNMetricsRemoved + 1),
       ylim=yLim, xlab="Number of Metrics Removed",
       yaxt="n", ylab="Relative Variable Importance", xaxt="n")
  axis(side=1, at=0:maxNMetricsRemoved)
  
  metrics <- names(variableImportance)
  
  for(metric in metrics){
    
    colour <- returnVariableColour(name=metric,
                                   spatial="red",
                                   temporal="darkgoldenrod4",
                                   network="blue", default="black")
    
    lines(x=variableImportance[[metric]][, 1],
          y=variableImportance[[metric]][, 2], type="o", col=colour)
  }
  
}

plotVariableImportanceAgainstNMetricsRemovedPerCluster <- function(variableImportance, 
                                                                   maxNMetricsRemoved,
                                                                   yLim){
  
  plot(x=0:maxNMetricsRemoved, y=rep(-1, maxNMetricsRemoved + 1),
       ylim=yLim, xlab="Number of Metrics Removed per Cluster",
       yaxt="n", ylab="Relative Variable Importance", xaxt="n")
  axis(side=1, at=0:maxNMetricsRemoved)
  
  metrics <- names(variableImportance)
  
  for(metric in metrics){
    
    colour <- returnVariableColour(name=metric,
                                   spatial="red",
                                   temporal="darkgoldenrod4",
                                   network="blue", default="black")
    
    lines(x=variableImportance[[metric]][, 1],
          y=variableImportance[[metric]][, 2], type="o", col=colour)
  }
  
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

noteVariableImportance <- function(importance, variableList, nMetricsPerClusterToRemove, colToUse){
  
  metrics <- rownames(importance)
  
  for(row in 1:nrow(importance)){
    
    # Check that metric exists in variable list
    if(is.null(variableList[[metrics[row]]]) == FALSE){
      
      variableList[[metrics[row]]] <- rbind(variableList[[metrics[row]]],
                                            c(nMetricsPerClusterToRemove,
                                              importance[row, colToUse]))
    }else{
      variableList[[metrics[row]]] <- matrix(data=c(nMetricsPerClusterToRemove,
                                                    importance[row, colToUse]),
                                             nrow=1, ncol=2)
    }
  }
  
  return(variableList)
}

indexArray <- function(array){
  
  output <- list()
  for(i in 1:length(array)){
    output[[array[i]]] <- i
  }
  return(output)
}

getClusterSizes <- function(clusters){
  
  sizes <- c()
  clusterNames <- names(clusters)
  
  for(i in 1:length(clusterNames)){
    
    sizes[i] <- length(clusters[[clusterNames[i]]])
  }
  
  return(sizes)
}

removeVariablesByImportanceFromClustersAndTable <- function(clusters, geneticVsEpi,
                                                    variableImportance, nToRemove,
                                                    useOriginalRfModelImportance){

  # Create an array to record which metrics have been removed
  allMetricsRemoved <- c()
  
  for(cluster in names(clusters)){
    
    # Check the number to be removed
    if(nToRemove == 0){
      next
    }
    
    variables <- clusters[[cluster]]
    
    # Get the variable importance
    varImportance <- c()
    for(i in 1:length(variables)){
      
      # Note if wanting to use original importance scores or those from previous model
      if(useOriginalRfModelImportance == TRUE){
        varImportance[i] <- variableImportance[[variables[i]]][1, 2]
      }else{
        varImportance[i] <- variableImportance[[variables[i]]][
          nrow(variableImportance[[variables[i]]]), 2]
      }
    }
    
    # Get the variable order by importance
    order <- order(varImportance)
    
    # Remove X lowest importance variables
    metricsRemoved <- c()
    if(nToRemove >= length(variables)){
      clusters[[cluster]] <- variables[order[length(variables)]]
      metricsRemoved <- variables[order[1:(length(variables) - 1)]]
    }else{
      clusters[[cluster]] <- variables[order[(nToRemove + 1):length(variables)]]
      metricsRemoved <- variables[order[1:nToRemove]]
    }

    # Remove the metrics from the geneticVsEpi table
    indicesOfMetricsToRemove <- which(names(geneticVsEpi) %in% metricsRemoved)
    geneticVsEpi <- geneticVsEpi[, -indicesOfMetricsToRemove]
    
    # Update the metrics removed vector
    allMetricsRemoved <- c(allMetricsRemoved, metricsRemoved)

  }
  
  output <- list(
    "clusters" = clusters,
    "geneticVsEpiTable" = geneticVsEpi,
    "removed" = allMetricsRemoved
  )
  
  return(output)
}

noteVariableImportanceOLD <- function(infoRF){
  
  variableImportance <- list()
  
  rowNames <- rownames(infoRF$importance)
  
  for(i in 1:nrow(infoRF$importance)){
    variableImportance[[rowNames[i]]] <- infoRF$importance[i, ]
  }
  
  return(variableImportance)
}

noteClustersOfMetrics <- function(correlationTable, threshold, defaultColour){

  # Create a list to record which clusters metrics are in
  clustersForMetrics <- list()
  cluster <- 0
  
  # Note the names of the metrics
  colNames <- colnames(correlationTable)
  
  # Examine each correlation between the metrics
  for(i in 1:nrow(correlationTable)){
    
    for(j in 1:ncol(correlationTable)){
      
      # Skip the diagonal
      if(i == j){
        next
      }
      
      # Check if correlation between the two metrics is above the threshold
      if(abs(correlationTable[i,j]) >= threshold){
        
        # Has metric i already been assigned to a cluster?
        if(is.null(clustersForMetrics[[colNames[i]]]) == FALSE && 
           is.null(clustersForMetrics[[colNames[j]]]) == TRUE){
          clustersForMetrics[[colNames[j]]] <- clustersForMetrics[[colNames[i]]]
          
          # Has metric j already been assigned to a cluster?
        }else if(is.null(clustersForMetrics[[colNames[j]]]) == FALSE &&
                 is.null(clustersForMetrics[[colNames[i]]]) == TRUE){
          clustersForMetrics[[colNames[i]]] <- clustersForMetrics[[colNames[j]]]
          
          # If neither have, create new cluster
        }else if(is.null(clustersForMetrics[[colNames[j]]]) == TRUE &&
                 is.null(clustersForMetrics[[colNames[i]]]) == TRUE){
          cluster <- cluster + 1
          clustersForMetrics[[colNames[i]]] <- cluster
          clustersForMetrics[[colNames[j]]] <- cluster
          
          # Both already in cluster - check that it is the same
        }else if(clustersForMetrics[[colNames[i]]] != clustersForMetrics[[colNames[j]]]){
          
          cat("-------------------------------------------------------------\n")
          cat(paste(colNames[i], ":", clustersForMetrics[[colNames[i]]], "\t",
                    colNames[j], ":", clustersForMetrics[[colNames[j]]], "\n", sep=""))
        }
      }
    }
  }
  
  # Summarise the clusters in figure
  plotClusterColours(cluster, clustersForMetrics, colNames, threshold, defaultColour)
  
  # Create a list of isolates in each cluster
  clusters <- list()
  for(metric in colNames){
    
    if(is.null(clustersForMetrics[[metric]]) == FALSE){
      
      # Have we encountered this cluster before?
      if(is.null(clusters[[as.character(clustersForMetrics[[metric]])]]) == FALSE){
        clusters[[as.character(clustersForMetrics[[metric]])]] <- 
          append(clusters[[as.character(clustersForMetrics[[metric]])]], metric)
      }else{
        clusters[[as.character(clustersForMetrics[[metric]])]] <- c(metric)
      }
    }
  }
  
  return(clusters)
}

plotClusterColours <- function(nClusters, clustersForMetrics, metricNames, threshold,
                               defaultColour){

  colourPalette <- colorRampPalette(c("red", "orange", "pink", "brown", "green", 
                                      "grey", "cyan", "purple", "violet", "blue"))
  colours = colourPalette(nClusters)
  
  #colours <- palette(rainbow(nClusters))
  
  metricColours <- c()

  for(i in 1:length(metricNames)){
    
    # Check if metric was assigned to cluster
    if(is.null(clustersForMetrics[[metricNames[i]]]) == TRUE){
      metricColours[i] <- defaultColour
    }else{
      metricColours[i] <- colours[clustersForMetrics[[metricNames[i]]]]
    }
  }

  plot.new()
  text(x=rep(0.5, length(metricNames)), y=seq(1, 0, -1/(length(metricNames) - 1)),
       labels=metricNames, 
       col=metricColours)
  legend("topleft", legend=(1:nClusters), text.col=colours, bty="n")
  legend("topright", legend=paste("Threshold =", threshold), bty="n")
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

calculateCorrelationBetweenEpiMetrics <- function(table, fullNames){
  
  # Note columns which are numeric
  numericColumns <- which(lapply(table, class) == "numeric")
  
  correlationTable <- round(cor(table[, numericColumns]), digits=2)
  
  plotHeatmap(correlationTable, fullNames)
 
  return(correlationTable)
}

plotHeatmap <- function(correlationTable, fullNames){
  colBreaks <- seq(-1, 1, by=0.1)
  
  # Plot the heatmap
  heatmap.2(correlationTable, # matrix is the input data
            
            # Show values
            cellnote=correlationTable,
            notecex=0.75,
            notecol="black",
            
            # Create the colour scale 
            col=colorpanel(n=length(colBreaks)-1, low="blue", mid="white", high="red"),
            breaks=colBreaks,
            
            # Turn off a density plot
            density.info="none", 
            
            # Turn off the trace
            trace="none",
            
            # Column Labels
            labCol=getFullVariableNames(colnames(correlationTable),
                                        fullNames),
            cexCol=0.4, # Change the size of the column labels
            srtCol=90, # Set the angle of the column labels (degrees from horizontal)
            offsetCol=-0.85, # Set size of space between column labels and heatmap
            
            # Change the size of the margins around the plot: c(column space, row space)
            margins = c(20, 20), 
            
            # Row labels
            labRow=getFullVariableNames(rownames(correlationTable),
                                        fullNames),
            cexRow=0.4, # Change the size of the Row labels
            offsetRow=0,
            
            # Make sure the order of the rows and columns is changed
            Rowv=TRUE, Colv=TRUE,
            
            # Don't plot any dendogram
            dendrogram="none",
            
            # Set up the Key
            key=FALSE, # Turn the key OFF
            
            # Say where to plot each part of the heatmap
            #     4     3
            #     2     1
            # 1. Heatmap
            # 2. Row Dendrogram
            # 3. Column Dendrogram
            # 4. Key
            lmat=rbind(4:3, 2:1),
            
            # Set the size of the spaces for output plots:
            #     Colour Key      |   Column Dendrogram
            #     -------------------------------------
            #     Row Dendrogram  |   Heatmap
            #
            # Note that these will be affected by the width and height you set for the PDF
            lhei=c(0.1, 10), # c(row1Width, row2Width)
            lwid=c(0.1, 10), # c(column1Width, column2Width)
            
            # Note that the input matrix is not symmetric
            symm = FALSE
  )
}

plotPredictedVersusActual <- function(actual, predicted, main="Predicted vs. Actual genetic distances"){
  
  # Plot the predicted versus actual
  plot(x=predicted, y=actual,
       las=1, pch=19, col=rgb(0,0,0, 0.1), xlab="Predicted", ylab="Actual", main=main)
  abline(lm(actual ~ predicted), col="red")
  
  # Calculate the correlation and r squared values
  correlation <- cor(actual, predicted)
  rSq <- correlation^2

  # Add a legend detailing the correlation and r squared values
  legend("topleft", c(paste("corr =", round(correlation, digits=2)),
                      paste("Rsq = ", round(rSq, digits=2))),
         bty="n", cex = 1)
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

plotEpidemiologicalMetricDistributionsWithMissingData <- function(table, fullNames,
                                                                  selection){
  
  # Get the column names of the genetic vs. epi metrics table
  colNames <- colnames(table)
  
  # Initialise an array to store the proportion of missing data for each metric
  propMissing <- calculateProportionMissingData(table)
  
  # Examine each epidemiological metric
  for(col in 1:ncol(table)){

    # Check if negative values are present in dataset
    if(length(which(as.numeric(table[, col]) < 0)) > 0){
      
      cat(paste("Found ", length(which(table[, col] < 0)), " missing entries for: ",
                colNames[col], "\n", sep=""))
      
      # Check if column is a factor
      if(is.factor(table[, col]) == TRUE){
        
        plot(table[, col], las=1, main=colNames[col])
      }else{
        hist(table[, col], las=1, main=colNames[col], xlab="value", breaks=200, 
             border="black", col="black")
      }
    }
  }
  
  # Get the order by proportion missing data
  order <- order(propMissing)
  
  # Get the full names of the metrics
  variableNames <- getFullVariableNames(colnames(table), fullNames)
  
  # Get the margin sizes - to fit in metric names
  marginSizes <- list(
    "BB" = 26,
    "CC" = 33,
    "CB" = 22,
    "BC" = 22
  )
  par(mar=c(0,marginSizes[[selection]],2,0.5)) # bottom, left, top, right
  
  plot <- barplot(propMissing[order], horiz=TRUE,
                  beside=TRUE,
                  main="Proportion Missing Data")
  
  at <- plot[,1]
  xLabPosition <- 0
  text(labels=variableNames[order],
       x=rep(xLabPosition,length(propMissing)),
       y=at, srt = 0, pos = 2, xpd = TRUE, cex=0.75)
  
  spacing <- at[2] - at[1]
  ticks <- seq(0,1,0.1)
  ticks <- ticks[ticks < max(propMissing)]
  axis(side=3, at=ticks, line=-spacing*1.2, mgp=c(3, .5, 0))
  
  # Reset margins
  par(mar=c(5.1,4.1,4.1,2.1))
}

makeBooleanColumnsFactors <- function(table){
  
  colNames <- colnames(table)
  
  for(col in 1:ncol(table)){
    
    if((grepl(x=colNames[col], pattern="Same") == TRUE &
       grepl(x=colNames[col], pattern="PeriodSpentIn") == FALSE) ||
       grepl(x=colNames[col], pattern="Boolean") == TRUE){
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

  columnsToConsider <- which(colnames(table) %in% c("GeneticDistance", "iSpeciesJSpecies",
                                                    "IsolateI", "IsolateJ") == FALSE)

  for(col in columnsToConsider){
    
    if(sd(table[, col]) == 0){
      colsToRemove[length(colsToRemove) + 1] <- col
      cat(paste("Removed: ", colnames(table)[col], "\n", sep=""))
    }
  }
  return(table[, -colsToRemove])
}
