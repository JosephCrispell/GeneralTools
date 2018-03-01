################################################################
# Issues with these genetic vs. epidemiological metrics:       #
# - Presence of -1 when data not available                     #
# - Flattening a matrix - multiple data from individual animal #
################################################################

##################
# Load libraries #
##################

library(randomForest)
library(gplots)

#######################################################
# Open the Genetic Vs Epidemiological distances table #
#######################################################

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_13-07-17/GeneticVsEpidemiologicalDistances/"

# Read in the table
file <- paste(path, "GeneticVsEpidemiologicalDistances_02-10-17.txt", sep="")
geneticVsEpi <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

####################
# General settings #
####################

selection <- "BB"
trainProp <- 0.5
colToUse <- "%IncMSE"

# Drop out genetic relatedness variable
if(selection == "BB"){
  col <- which(colnames(geneticVsEpi) == "HostRelatedness")
  geneticVsEpi <- geneticVsEpi[, -col]
}

# Note the full names of metrics and assign them a colour
fullNames <- noteFullNames()
temporalCol="darkgoldenrod4"
spatialCol="red"
networkCol="blue"
nameColours <- assignMetricColours(temporalCol=temporalCol, spatialCol=spatialCol,
                                   networkCol=networkCol)

###############
# Select data #
###############

###### Open a PDF
file <- paste(path, "ExamineEpiVariableCorrelation_", selection, "_02-10-17.pdf", sep="")
pdf(file, height=10, width=10)

par(mfrow=c(1,1))

# Subset out the selected comparisons
geneticVsEpi <- selectAppropriateComparisonsForSelection(selection)

# Only select small genetic distances
geneticVsEpi <- selectGeneticDistancesBelowThreshold(threshold=15)

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

##################################################
# Examine the extent of the correlated variables #
##################################################

# Examine correlation between epidemiological metrics
correlationTable <- calculateCorrelationBetweenEpiMetrics(
  geneticVsEpi[, -c(1, colsToIgnore)], fullNames)

# Find clusters of highly correlated metrics
threshold <- 0.55
clusters <- noteClustersOfMetrics(correlationTable, threshold, "black")
clusterSizes <- getClusterSizes(clusters)

#############################################################################
# Investigate the effect of metric removal on R squared and metric rankings #
#############################################################################

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

########################################################
# Examine the metric importance in random forest model #
########################################################

plotVariableImportance(infoRF=infoRF, colToUse=colToUse, fullNames=fullNames, 
                       nameColours=nameColours,
                       temporalCol=temporalCol, spatialCol=spatialCol,
                       networkCol=networkCol, showY=FALSE)

plotVariableImportance(infoRF=infoRF, colToUse=colToUse, fullNames=fullNames, 
                       nameColours=nameColours,
                       temporalCol=temporalCol, spatialCol=spatialCol,
                       networkCol=networkCol, showY=TRUE)


###### Close PDF
dev.off()


#############
# FUNCTIONS #
#############

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

plotVariableImportance <- function(infoRF, colToUse, fullNames, nameColours,
                                   temporalCol, spatialCol, networkCol, showY){
  
  # Get the variable importance from the RF model
  variableImportance <- as.data.frame(infoRF$importance)
  
  # Order the table by importance
  variableImportance <- variableImportance[order(variableImportance[, colToUse],
                                                 decreasing=FALSE), ]
  
  # Transpose the table
  transpose <- as.matrix(t(variableImportance))
  
  # Get full Variable Names
  variableNames <- getFullVariableNames(rownames(variableImportance), fullNames)
  
  # Get Variable Colours
  variableColours <- getVariableColours(variableImportance, nameColours=nameColours)
  
  # Create bar plot illustrating the relative variable importance from RF and BR
  par(mfrow=c(1,1))
  
  marginSizes <- list(
    "BB" = 26,
    "CC" = 33,
    "CB" = 22
  )
  
  legendPos <- list(
    "BB" = c(2, 5),
    "CC" = c(0.5, 2.5),
    "CB" = c(8, 1.5)
  )
  
  par(mar=c(0,marginSizes[[selection]],2,0.5)) # bottom, left, top, right
  
  if(showY == TRUE){
    par(mar=c(2,marginSizes[[selection]],2,0.5)) # bottom, left, top, right
    
  }
  
  plot <- barplot(transpose[-2,], horiz=TRUE, beside=TRUE,
                  xaxt="n",
                  col=variableColours,
                  main="Variable Importance",
                  col.axis="white")
  
  at <- plot[,1]
  
  xLabPosition <- 0
  text(labels=variableNames, 
       col=variableColours,
       x=rep(xLabPosition,length(variableImportance)),
       y=at,
       srt = 0, pos = 2, xpd = TRUE, cex=0.75)
  
  if(showY == TRUE){
    axis(side=1, line=-0.5, cex.axis=0.75, mgp=c(3, .25, 0))
    mtext("% Increase MSE", side=1, line=0.75)
  }
  
  
  # Add Legend
  legend(x=legendPos[[selection]][1], y=legendPos[[selection]][2],
         legend=c("Temporal", "Spatial", "Network"), 
         text.col = c(temporalCol, spatialCol, networkCol),
         bty='n', cex=1)
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

plotPredictedVersusActual <- function(actual, predicted, rSq){
  plot(x=predicted, y=actual,
       las=1, pch=19, col=rgb(0,0,0, 0.1), xlab="Predicted", ylab="Actual")
  abline(lm(actual ~ predicted), col="red")
  correlation <- round(cor(actual, predicted), digits=2)
  
  legend("topleft", c(paste("corr =", correlation),
                      paste("Rsq = ", rSq)),
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
