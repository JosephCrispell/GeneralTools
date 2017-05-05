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

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/GeneticVsEpidemiologicalDistance/ExaminingRandomForestFit_06-04-17/"

# Read in the table
file <- paste(path, "GeneticVsEpidemiologicalDistances_12-04-17.txt", sep="")
geneticVsEpi <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

#####
# General Settings

trainProp <- 0.5

fullNames <- list(
  "SameMainGroup" = "Isolates taken from same main group (yes/no)",
  "SameSampledGroup" = "Isolates taken from same sampled group (yes/no)",                           
  "SameInfectedGroup" = "Isolates taken from same infected group (yes/no)",
  "PeriodSpentAliveTogether" = "Number of days overlap between recorded lifespans of the sampled badgers",                   
  "PeriodSpentInfectedTogether" = "Number of days overlap betweeninfected lifespans of the sampled badgers",
  "PeriodSpentInSameGroup" = "Number of days that sampled badgers spent in same group",                     
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
  "NSharedAnimalsBetweenMainGroups" = "Number of badgers captured in both main groups of sampled badgers",            
  "NSharedAnimalsBetweenSampledGroups" = "Number of badgers captured in both sampled groups of sampled badgers",
  "NSharedAnimalsBetweenInfectedGroups" = "Number of badgers captured in both infected groups of sampled badgers",
  "ShortestPathLengthEXCLMain" = "Shortest path length between main herds of sampled animals (Some Herds Excluded)",                     
  "MeanNMovementsOnEdgesOfShortestPathEXCLMain" = "Mean number of animals dispersing along edges of shortest path between main herds (Some Herds Excluded)",   
  "ShortestPathLengthEXCLSampled" = "Shortest path length between main herds of sampled animals (Some Herds Excluded)",               
  "MeanNMovementsOnEdgesOfShortestPathEXCLSampled" = "Mean number of animals dispersing along edges of shortest path between sampled herds (Some Herds Excluded)",
  "ShortestPathLengthEXCLInfected" = "Shortest path length between main herds of sampled animals (Some Herds Excluded)",
  "MeanNMovementsOnEdgesOfShortestPathEXCLInfected" = "Mean number of animals dispersing along edges of shortest path between main herds (Some Herds Excluded)",
  "CentroidDistBetweenMain" = "Distance from centroid of closest land parcel to badgers main sett",
  "CentroidDistBetweenSamp" = "Distance from centroid of closest land parcel to badgers sampled sett"
)

#####

### Badger - Badger #

# Open a PDF
file <- paste(path, "ExamineEpiVariableCorrelation_02-05-17.pdf", sep="")
pdf(file, height=10, width=10)

par(mfrow=c(1,1))

# Select Data
#####


# Subset out the badger - badger comparisons
geneticVsEpi_BB <- geneticVsEpi[geneticVsEpi$iSpeciesJSpecies == "BB", ]

# Only select small genetic distances
hist(geneticVsEpi_BB$GeneticDistance, breaks=100,
     las=1,
     xlab="Genetic Distance (SNPs)",
     main="Inter-Isolate Genetic Distance Distribution")
lines(x=c(15, 15), y=c(0, 8500), col="red", lty=2)

geneticVsEpi_BB <- geneticVsEpi_BB[geneticVsEpi_BB$GeneticDistance < 15, ]

# Remove irrelevant columns
geneticVsEpi_BB <- removeColumnsIfNotRelevant(geneticVsEpi_BB)

# Convert the columns dealing with boolean metrics to factors
geneticVsEpi_BB <- makeBooleanColumnsFactors(table=geneticVsEpi_BB)

# Note columns to ignore as predictors
colNamesToIgnore <- c("iSpeciesJSpecies", "IsolateI", "IsolateJ")
colsToIgnore <- which(names(geneticVsEpi_BB) %in% colNamesToIgnore)

# Examine the levels of missing data
# -1 restricted to shortest paths and periods spent together
# For shortest paths, this means no path found
# For periods spent together - means at least one animal doesn't have start and end dates
# Decided to LEAVE IN
plotEpidemiologicalMetricDistributionsWithMissingData(
  geneticVsEpi_BB[, -colsToIgnore])


#####

# Fit Random Forest Model
#####


# Build a test and training data set for predictions
trainRows_BB <- sample(x=1:nrow(geneticVsEpi_BB),
                    size=floor(trainProp * nrow(geneticVsEpi_BB)), replace=FALSE)

# Find optimal mtry parameter
optimalMtry <- findOptimalMtry(response=geneticVsEpi_BB[trainRows_BB, "GeneticDistance"],
                               predictors=geneticVsEpi_BB[trainRows_BB, -c(1, colsToIgnore)],
                               mTryInitial=3, nTrees=500, plot=TRUE)
abline(v=optimalMtry, col="red", lty=2)

# Train the Random Forest model
infoRF_BB <- randomForest(geneticVsEpi_BB[trainRows, "GeneticDistance"] ~ ., 
                          data=geneticVsEpi_BB[trainRows, -c(1, colsToIgnore)],
                          mtry=optimalMtry, importance=TRUE, ntree=1000,
                          keep.forest=TRUE, norm.votes=FALSE, proximity=FALSE,
                          do.trace=FALSE)

# Get the Pseudo RSquared value
rSq <- round(infoRF_BB$rsq[length(infoRF_BB$rsq)], digits=2)

# Examine trained model prediction
predictedGeneticDistances <- predict(infoRF_BB, geneticVsEpi_BB[-trainRows, -c(1, colsToIgnore)])
plotPredictedVersusActual(actual=geneticVsEpi_BB[-trainRows, "GeneticDistance"],
                          predicted=predictedGeneticDistances, rSq=rSq)

# Note the metric importance
epiMetricImportance <- noteVariableImportance(infoRF_BB)

#####

# Examine Epidemiological Metric Correlations
#############################################


# Examine correlation between epidemiological metrics
correlationTable <- calculateCorrelationBetweenEpiMetrics(
  geneticVsEpi_BB[, -c(1, colsToIgnore)], fullNames)

# Find clusters of highly correlated metrics
threshold <- 0.55
clusters <- noteClustersOfMetrics(correlationTable, threshold, "black")

#####

# Examine the Effect of Removing Highly Correlated Variables
############################################################


# Investigate the effect of removing coorelated variables on Random Forest fit
stop <- FALSE
nMetricsPerClusterToRemove <- 0
par(mfrow=c(1,2))

outputTable <- data.frame("NumberMetricsRemovedPerCluster" = c(),
                          "MetricsRemoved" = c(),
                          "NumberMetricsRemoved" = c(),
                          "PseudoRSquared" = c(),
                          "Correlation" = c(), 
                          "ClusterSizes" = c(), 
                          stringsAsFactors=FALSE)

# Create a list to store the variable importance after variable removal
variableImportance <- list()
maxImportance <- 0

while(stop == FALSE){

  # Remove x metrics from clusters and geneticVsEpi table
  inforForRFModel <- removeVariablesByImportanceFromClustersAndTable(
    clusters=clusters,
    geneticVsEpi=geneticVsEpi_BB[, -colsToIgnore],
    variableImportance=epiMetricImportance,
    nToRemove=nMetricsPerClusterToRemove)

  # Find the optimal mtry value
  optimalMtry <- findOptimalMtry(
    response=inforForRFModel[["geneticVsEpiTable"]][trainRows, "GeneticDistance"],
    predictors=inforForRFModel[["geneticVsEpiTable"]][trainRows, -1],
    mTryInitial=3, nTrees=500, plot=TRUE)
  abline(v=optimalMtry, col="red", lty=2)

  # Fit the Random Forest model
  rfModel <- randomForest(
    inforForRFModel[["geneticVsEpiTable"]][trainRows, "GeneticDistance"] ~ .,
    data=inforForRFModel[["geneticVsEpiTable"]][trainRows, -1],
    mtry=optimalMtry, importance=TRUE, ntree=1000,
    keep.forest=TRUE, norm.votes=FALSE, proximity=FALSE,
    do.trace=FALSE)
  rSq <- rfModel$rsq[length(rfModel$rsq)]
  plot(rfModel, las=1)

  # Note the variable importance
  variableImportance <- noteVariableImportance(importance=rfModel$importance,
                                               variableList=variableImportance,
                                               nMetricsRemovedPerCluster=nMetricsPerClusterToRemove,
                                               colToUse="%IncMSE")
  if(max(rfModel$importance[, "%IncMSE"]) > maxImportance){
    maxImportance <- max(rfModel$importance[, "%IncMSE"])
  }
  
  
  # Test the Random Forest model
  predictions <- predict(rfModel,
                         inforForRFModel[["geneticVsEpiTable"]][-trainRows, -1])
  corr <- cor(inforForRFModel[["geneticVsEpiTable"]][-trainRows, "GeneticDistance"]
              , predictions)

  # Examine the size of the clusters following metric removal
  clusterSizes <- getClusterSizes(inforForRFModel[["clusters"]])
  if(max(clusterSizes) == 1){
    stop = TRUE
  }

  # Store the results
  outputTable[nMetricsPerClusterToRemove + 1, "NumberMetricsRemovedPerCluster"] <- nMetricsPerClusterToRemove
  outputTable[nMetricsPerClusterToRemove + 1, "MetricsRemoved"] <- paste(inforForRFModel[["removed"]], collapse=",")
  outputTable[nMetricsPerClusterToRemove + 1, "NumberMetricsRemoved"] <- length(inforForRFModel[["removed"]])
  outputTable[nMetricsPerClusterToRemove + 1, "PseudoRSquared"] <- rSq
  outputTable[nMetricsPerClusterToRemove + 1, "Correlation"] <- corr
  outputTable[nMetricsPerClusterToRemove + 1, "ClusterSizes"] <- paste(clusterSizes, collapse=",")

  # Print progress information
  cat("###########################################################################\n")
  cat("###########################################################################\n")
  cat(paste("Finished fitting Random Forest model", "\n",
            "Removed ", length(inforForRFModel[["removed"]])," metrics (",
            nMetricsPerClusterToRemove, " per cluster):\n",
            paste(inforForRFModel[["removed"]], collapse=","), "\nCorrelation = ",
            round(corr, digits=2), "\nRsq = ", round(rSq, digits=2), "\n", sep=""))
  cat("###########################################################################\n")
  cat("###########################################################################\n")

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
plotVariableImportanceAgainstNMetricsRemovedPerCluster(variableImportance, 
                                                       nMetricsPerClusterToRemove - 1,
                                                       c(0, maxImportance))
#####

dev.off()


#############
# FUNCTIONS #
#############

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

noteVariableImportance <- function(importance, variableList, nMetricsRemovedPerCluster, colToUse){
  
  metrics <- rownames(importance)
  
  for(row in 1:nrow(importance)){
    
    # Check that metric exists in variable lise
    if(is.null(variableList[[metrics[row]]]) == FALSE){
      
      variableList[[metrics[row]]] <- rbind(variableList[[metrics[row]]],
                                            c(nMetricsRemovedPerCluster,
                                              importance[row, colToUse]))
    }else{
      variableList[[metrics[row]]] <- matrix(data=c(nMetricsRemovedPerCluster,
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
                                                    variableImportance, nToRemove){

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
      varImportance[i] <- variableImportance[[variables[i]]][1]
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

noteVariableImportance <- function(infoRF){
  
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
            cexCol=0.5, # Change the size of the column labels
            srtCol=90, # Set the angle of the column labels (degrees from horizontal)
            offsetCol=-0.85, # Set size of space between column labels and heatmap
            
            # Change the size of the margins around the plot: c(column space, row space)
            margins = c(20, 20), 
            
            # Row labels
            labRow=getFullVariableNames(rownames(correlationTable),
                                        fullNames),
            cexRow=0.5, # Change the size of the Row labels
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
                      plot=plot)
  
  return(as.integer(rownames(tuneOutput)[tuneOutput[,2] == min(tuneOutput[,2])]))
}

plotEpidemiologicalMetricDistributionsWithMissingData <- function(table){
  
  colNames <- colnames(table)
  
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
    
    if(mean(table[, col]) == -1){
      index <- index + 1
      colsToRemove[index] <- col
    }
  }
  return(table[, -colsToRemove])
}
