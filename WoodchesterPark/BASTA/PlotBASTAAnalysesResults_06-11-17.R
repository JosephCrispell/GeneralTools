##########
# Set up #
##########

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_13-07-17/"

#################################
# Note the Genome Size examined #
#################################

# Get the constant site counts
constantSiteCountsFile <- paste(path, "vcfFiles/", "constantSiteCounts_29-09-2017.txt", sep="")
constantSiteCounts <- getConstantSiteCounts(constantSiteCountsFile)
genomeSize <- sum(constantSiteCounts) + 8893

#####################
# Read in the files #
#####################

# Move the path
path <- paste(path, "BASTA/Replicate1/", sep="")

# Note deme structure to use
demeStructureDates <- list(
  "2Deme"="09-10-17",
  "3Deme-outerIsBoth"="09-10-17",
  "3Deme-outerIsCattle"="09-10-17",
  "4Deme"="09-10-17",
  "6Deme-EastWest"="13-10-17",
  "6Deme-NorthSouth"="13-10-17",
  "8Deme-EastWest"="13-10-17",
  "8Deme-NorthSouth"="13-10-17"
)

# Note the population size estimation options
popEstimationTypes <- c("varying", "equal")

# Note the clock model options
clockEstimateTypes <- c("relaxed") # strict not used

# Store each of the log tables in a list
logTables <- readInBASTALogTables(demeStructureDates, popEstimationTypes, clockEstimateTypes,
                                  path)

####################
# Examine each run #
####################

# Create summary plots, note migration rates and calculate AICM scores for each model
nBootstraps <- 1000
migrationRateEstimates <- summarisePosteriorLogTables(path, logTables, code=2,
                                                      arrowFactor=20, nBootstraps,
                                                      genomeSize)

#################################
# Plot Model comparison results #
#################################

# Open a PDF
file <- paste(path, "SummaryFiguresOfModelEstimations_18-12-17.pdf", sep="")
pdf(file)

# Examine the model likelihoods
plotModelAICMScores(migrationRateEstimates, nBootstraps)

######################################################################
# Produce single summary plot of rate estimation: BADGERS <-> CATTLE #
######################################################################

# Calculate the weighted mean estimated transition rates between cattle and badger demes
# weighted by the AICM model scores
weightedMeanEstimates <- 
  calculateMeanEstimatedTransitionRatesBetweenCattleAndBadgerPopulationsWeightedByAICM(
    migrationRateEstimates)

# Plot the rates
plotTransitionRatesBetweenBadgerAndCow(badgerToCow=weightedMeanEstimates[1],
                                       cowToBadger=weightedMeanEstimates[2],
                                       code=2, arrowFactor=20)

##############################################################################
#                                                                            #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# !!!                                                                    !!! #
# !!! Combining the posteriors based upon the AICM weight probabilities? !!! #
# !!!                                                                    !!! #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
#                                                                            #
##############################################################################


# Close the pdf
dev.off()

#############
# FUNCTIONS #
#############

getConstantSiteCounts <- function(constantSiteCountsFile){
 
  # Open the file and store the file lines
  connection <- file(constantSiteCountsFile, "r")
  lines <- readLines(connection)
  close(connection)
  
  # Calculate the constant site counts
  counts <- c(0,0,0,0)
  for(line in lines){
    
    # Skip lines without counts
    if(grepl(line, pattern="Counts|A") == TRUE || line == ""){
      next
    }
    
    # Split the line into its parts
    parts <- strsplit(line, "\t")[[1]]
    
    # Add the counts to the tally
    counts <- counts + as.numeric(parts)
  }
  
  return(counts)
}

plotParameterTraces <- function(logTable, colNamesToPlot){
  
  # Define a grid of plots
  nPlots <- ceiling(sqrt(length(colNamesToPlot)))
  
  # Set the plotting window dimensions
  par(mfrow=c(nPlots, nPlots))
  
  # Set the margins
  par(mar=c(0, 0, 2, 0))
  
  # Plot a trace for each parameter
  for(col in colNamesToPlot){
    plot(logTable[, col], type="l", xlab="", xaxt="n", ylab="", yaxt="n", bty="n",
         main=col, cex.main=0.5)
  }
  
  # Reset the margins
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  
  # Reset the plotting window dimensions
  par(mfrow=c(1,1))
}

plotTransitionRatesBetweenBadgerAndCow <- function(badgerToCow, cowToBadger, code, 
                                                   arrowFactor){
  
  # Create empty plot
  createEmptyPlot()
  
  # Get deme names and assign colours
  demeNames <- c("badgers", "cattle")
  demeColours <- c("red", "blue")
  
  # Normalise the values between 0 and arrowFactor
  values <- c(badgerToCow, cowToBadger)
  values <- (values / max(values)) * arrowFactor

  # Add labels
  x <- c(0.1, 0.9)
  y <- c(0.5, 0.5)
  text(x=x, y=y, 
       labels=demeNames,
       col=demeColours, cex=2)
  
  # badger -> cow
  arrows(x0=x[1]+0.15, x1=x[2]-0.15, y0=y[1]-0.15, y1=y[2]-0.15,
         code=code, lwd=values[1])
  
  # cow -> badger
  arrows(x0=x[2]-0.15, x1=x[1]+0.15, y0=y[2]+0.15, y1=y[1]+0.15,
         code=code, lwd=values[2])
  
  legend("topleft", legend=c("Badgers -> Cattle",
                       paste("-------------- =", round(badgerToCow/cowToBadger, digits=2)),
                             "Cattle -> Badgers"),
         bty="n", cex=1)
  
  # Reset the margins
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

calculateMeanEstimatedTransitionRatesBetweenCattleAndBadgerPopulationsWeightedByAICM <-
  function(migrationRateEstimates){
    
    # Calculating the weighted mean rates of transitions between badger and cattle demes
    # across BASTA models
    # Weighting is by the AICM score
    #
    # 1. For each analysis calculate the mean transition rates between cattle and badger demes
    # 2. Store those rates for each model in array
    # 3. Get the AICM score for each model
    # 4. Convert AICM scores to weights:
    #       exp((AIC_min - AIC_i)/2)
    # 5. Normalise to sum to 1: 
    #       scores / sum(scores)
    # 6. Calculate weighted average rates:
    #       sum(ratesFromModels * modelAICMWeightsNormalised)
    #
    # Above method taken from Paul Johnson
    # This is known as Ensemble Bayesian Model Averaging
    
    # Get the analysis names
    analyses <- names(migrationRateEstimates)
    
    # Get the AICM values
    aicmScores <- c()
    shortenedNames <- c()
    for(i in 1:length(analyses)){
      aicmScores[i] <- migrationRateEstimates[[analyses[i]]][["AICM"]][1]
      
      parts <- strsplit(analyses[i], split="_")[[1]]
      shortenedNames[i] <- paste(parts[1], parts[2], sep="_")
    }
    
    # Convert the AICM scores into weights
    modelAICMWeights <- exp((min(aicmScores) - aicmScores)/2)
    
    # Normalise the AICM weights
    normalisedModelAICMWeights <- modelAICMWeights / sum(modelAICMWeights)                           
    
    # Initialise two arrays to store the mean transmission rates from each model
    modelMeanBadgerToCowRates <- c()
    modelMeanCowToBadgerRates <- c()

    # Examine each of the different model structures
    for(i in 1:length(analyses)){
      
      # Initialise arrays to store the rates from badgers to cattle and vice versa
      ratesBadgerToCow <- c()
      ratesCowToBadger <- c()
      
      # Get the demeStructure
      demeStructure <- strsplit(analyses[i], "_")[[1]][1]
      
      # Examine each rate
      for(key in names(migrationRateEstimates[[analyses[i]]])){
        
        # Ignore AICM
        if(key == "AICM"){
          next
        }
        
        # Split the key into its deme numbers
        demeNumbers <- as.numeric(strsplit(key, split="_")[[1]])
        
        # Check if current rate is between badger and cattle populations
        if(grepl(getDemeNamesForDemeStructure(demeStructure, demeNumbers[1]),
                 pattern="badger") == TRUE &&
           grepl(getDemeNamesForDemeStructure(demeStructure, demeNumbers[2]),
                 pattern="cow") == TRUE){
          
          ratesBadgerToCow[length(ratesBadgerToCow) + 1] <- 
            migrationRateEstimates[[analyses[i]]][[key]]
          
        }else if(grepl(getDemeNamesForDemeStructure(demeStructure, demeNumbers[1]),
                       pattern="cow") == TRUE &&
                 grepl(getDemeNamesForDemeStructure(demeStructure, demeNumbers[2]),
                       pattern="badger") == TRUE){
          
          ratesCowToBadger[length(ratesCowToBadger) + 1] <- 
            migrationRateEstimates[[analyses[i]]][[key]]
        }
      }
      
      # Calculate mean rates between cattle and badgers
      meanRateBadgerToCow <- mean(ratesBadgerToCow)
      modelMeanBadgerToCowRates[
        length(modelMeanBadgerToCowRates) + 1] <- 
        meanRateBadgerToCow
      meanRateCowToBadger <- mean(ratesCowToBadger)
      modelMeanCowToBadgerRates[length(modelMeanCowToBadgerRates) + 1] <- meanRateCowToBadger
    }
    
    # Calculate the weighted means for the rates between badgers and cattle
    weightedMeanBadgerToCow <- sum(modelMeanBadgerToCowRates * 
                                     normalisedModelAICMWeights)
    weightedMeanCowToBadger <- sum(modelMeanCowToBadgerRates * 
                                     normalisedModelAICMWeights)

    # Plot the AICM Scores against the rates
    par(mar=c(1, 7, 0, 0))
    plot(x=NULL, y=NULL,
         ylim=range(aicmScores), 
         xlim=range(c(modelMeanBadgerToCowRates, modelMeanCowToBadgerRates)),
         las=1, bty="n", yaxt="n", xaxt="n", ylab="", xlab="")
    axis(side=2, at=aicmScores, labels=shortenedNames, las=1,
         cex.axis=0.5, tick=FALSE, line=-1.5)
    for(i in 1:length(analyses)){
      lines(x=c(0, modelMeanBadgerToCowRates[i]),
            y=c(aicmScores[i], aicmScores[i]), lty=2)
    }
    points(x=modelMeanBadgerToCowRates, y=aicmScores,
           col=rgb(1,0,0, 0.75), pch=19, cex=2)
    points(x=modelMeanCowToBadgerRates, y=aicmScores, col=rgb(0,0,1, 0.75), pch=17, cex=2)
    
    mtext("AICM Score", side=2, line=5)
    mtext("Estimated Transition Rates", side=1, line=-0.5)
    
    legend("topright", legend=c("badgers-to-cattle", "cattle-to-badgers"),
           pch=c(19, 17), text.col=c("red", "blue"), bty="n", 
           col=c("red", "blue"))
    
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    
    return(c(weightedMeanBadgerToCow, weightedMeanCowToBadger))
}

plotModelAICMScores <- function(migrationRateEstimates, nBootstraps){
  
  # Get the deme structure names
  analyses <- names(migrationRateEstimates)
  
  # Initialise an array to store the shortened analyses names
  names <- c()
  
  # Get the AICM values
  aicmScores <- c()
  bootstrapMins <- c()
  bootstrapMaxs <- c()
  for(i in 1:length(analyses)){
    
    values <- migrationRateEstimates[[analyses[i]]][["AICM"]]
    
    aicmScores[i] <- values[1]
    bootstrapMins[i] <- values[2]
    bootstrapMaxs[i] <- values[3]
    
    parts <- strsplit(analyses[i], split="_")[[1]]
    names[i] <- paste(parts[1], parts[2], sep="_")
  }
  
  # Calculate the range of the AICM scores - y axis limits
  range <- range(aicmScores, bootstrapMins, bootstrapMaxs)
  yLim <- c(range[1] - (0.1 * (range[2] - range[1])),
            range[2] + (0.1 * (range[2] - range[1])))
  
  # Set the margins
  par(mar=c(10, 7, 0, 0.1))
  
  # Set bar colours to highlight min
  colours <- rep("gray", length(aicmScores))
  min <- min(aicmScores)
  minIndex <- which(aicmScores == min)
  colours[minIndex] <- "red"
  
  # Create an empty plot
  plot(x=NULL, y=NULL, bty="n", 
       xlim=c(0.5, length(analyses) + 1), ylim=yLim,
       xaxt="n", yaxt="n", ylab="", xlab="")
  axis(side=2, at=seq(yLim[1], yLim[2], by=(yLim[2] - yLim[1])/5), las=2)
  mtext(side=2, text="AICM Scores", line=5.5, cex=1.5)
  axis(side=1, at=1:length(analyses),
       labels=names, las=2, tick=FALSE, cex.axis=0.75)
  
  # Add the bars
  for(i in 1:length(analyses)){
    polygon(x=c(i-0.4, i-0.4, i+0.4, i+0.4),
            y=c(0, aicmScores[i], aicmScores[i], 0),
            col=colours[i])
    lines(x=c(i, i),
          y=c(bootstrapMins[i], bootstrapMaxs[i]))
  }
  
  # Add some horizontal lines
  points(x=c(1, length(analyses)),
         y=c(min, min), lty=2, col="red", type="l")
  points(x=c(1, length(analyses)),
         y=c(bootstrapMins[minIndex], bootstrapMins[minIndex]),
         lty=2, col=rgb(0,0,0, 0.5), type="l")
  points(x=c(1, length(analyses)),
         y=c(bootstrapMaxs[minIndex], bootstrapMaxs[minIndex]),
         lty=2, col=rgb(0,0,0, 0.5), type="l")

  # Reset the margins
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

summarisePosteriorLogTables <- function(path, logTables, code, arrowFactor, nBootstraps,
                                        genomeSize){
  
  # Get analysis names
  analyses <- names(logTables)
  
  # Initialise a list to store the migration rate estimates and AICM
  migrationRateEstimates <- list()
  
  # Examine each analysis
  for(analysis in analyses){
    
    # Open a pdf file
    file <- paste(path, analysis, "/", analysis, "_ResultsSummary.pdf", sep="")
    pdf(file)
    
    # Get the deme structure
    demeStructure <- strsplit(analysis, "_")[[1]][1]
    
    # Note progress
    cat(paste("\rExamining log table for: ", demeStructure, "                    ", sep=""))
    
    # Get the log table
    logTable <- logTables[[analysis]]
    
    # Note which parameters to examine posterior support for
    colsToCalculateESS <- colnames(logTable)[
      grepl(x=colnames(logTable), pattern="Sample|rateMatrixFlag|forward|Count") == FALSE]

    # Plot the ESS values of each parameter estimated
    plotParameterESSValues(logTable, colsToCalculateESS)
    
    # Plot the posterior support for each deme as source
    plotPosteriorSupportForEachDemeAsRoot(logTable, demeStructure)
    
    # Plot the population size estimates
    plotPopulationSizes(logTable, demeStructure, alpha=0.5)
    
    # Examine the substitution rate estimates
    examineSubstitutionRateEstimates(logTable, genomeSize)
    
    # Produce a migration rate estimation figure - weight by rate flags
    # Diagrams designed with code = 2 (FORWARDS) in mind
    migrationRateEstimates[[analysis]] <- 
      plotMigrationRates(logTable, demeStructure, code, arrowFactor)
    
    # Plot the parameter traces
    plotParameterTraces(logTable, colsToCalculateESS)
    
    # Calculate the acim
    migrationRateEstimates[[analysis]][["AICM"]] <- 
      calculateAICM(logTable$treeLikelihood1, nBootstraps)
    
    
    # Close the pdf file
    dev.off()
  }
  cat("\rFinished examining log tables...          ")
  
  # Reset the margins
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  
  return(migrationRateEstimates)
}

examineSubstitutionRateEstimates <- function(logTable, genomeSize){
  
  # Get the substitution rate estimates note that stored differently between strict and relaxed
  rateEstimates <- c()
  if(length(which(grepl(colnames(logTable), pattern="mutationRate") == TRUE)) == 1){
    rateEstimates <- logTable$mutationRate
  }else{
    rateEstimates <- logTable$ucedMean
  }
  
  # Define a colour palette based upon the likelihood
  rbPal <- colorRampPalette(c('red','blue'))
  colours <- rbPal(10)[as.numeric(cut(logTable$posterior,breaks = 10))]
  colours <- setAlphas(colours, 0.5)

  # Plot the substitution rate distribution versus the root height
  plot(x=rateEstimates * genomeSize, y=logTable$tree.height, pch=20, 
       xlab="Substitution Rate (per Genome per Year)", ylab="Root Height (years)",
       main="Substitution Rate versus Root Height",
       col=colours, las=1, bty="n")
  legend("topright", legend=c("High", "Low"), pch=20, col=c("red", "blue"),
         bty="n")
}

plotPopulationSizes <- function(logTable, demeStructure, alpha=0.5){
  
  # Get the population size distributions
  popSizeEstimates <- list()
  for(col in colnames(logTable)){
    
    if(grepl(col, pattern="popSize") == FALSE){
      next
    }
    
    popSizeEstimates[[strsplit(col, split="_")[[1]][2]]] <- logTable[, col]
  }
  
  # Get the deme names
  demeNames <- getDemeNamesForDemeStructure(demeStructure)
  
  # Define the limits of the population sizes and histogram breaks
  values <- c()
  for(key in names(popSizeEstimates)){
    values <- c(values, popSizeEstimates[[key]])
  }
  xLim <- range(values)
  breaks <- seq(from=xLim[1], to=xLim[2] + 5, by=5)
  
  # Create a histogram object for each distribution of sizes
  histograms <- list()
  for(key in names(popSizeEstimates)){
    
    histograms[[key]] <- hist(popSizeEstimates[[key]], plot=FALSE, breaks=breaks)
  }
  
  # Get the y axis limits
  counts <- c()
  for(key in names(histograms)){
    counts <- c(counts, histograms[[key]]$counts)
  }
  yLim <- c(0, max(counts))
  
  # Create a vector of colours
  colours <- c("red", "blue", "green", "cyan", "orange", "darkorchid4", "deeppink", "black", "brown", "darkolivegreen")
  
  # Plot the histograms
  keys <- names(histograms)
  for(i in 1:length(keys)){
    
    if(i == 1){
      plot(histograms[[keys[i]]], 
           col=setAlpha(colours[i], alpha), las=1, 
           main="Deme Effective Population Sizes", xlab="Effective Size",
           ylim=yLim, xlim=xLim)
    }else{
      plot(histograms[[keys[i]]], col=setAlpha(colours[i], alpha), add=TRUE)
    }
  }
  
  # Add legend
  legend("topright", legend=demeNames, bty="n",
         text.col=colours)
  
}

setAlphas <- function(colours, alpha){
  output <- c()
  for(i in 1:length(colours)){
    output[i] <- setAlpha(colours[i], alpha)
  }
  return(output)
}

setAlpha <- function(colour, alpha){
  
  # Convert the input colour into rgb values
  rgbValues <- col2rgb(colour)
  
  # Place rgb values within rgb function and insert alpha value
  # Note that col2rgb returns rgbvlues from 0 to 255
  rgbColour <- rgb(rgbValues["red", 1], rgbValues["green", 1], rgbValues["blue", 1],
                   alpha=alpha*255, maxColorValue=255)
  return(rgbColour)
}

calculateAICM <- function(logLikelihoodValues, nBootstraps){
  # Calculation taken from:
  # Baele et al. 2012 - Improving the accuracy of demographic and molecular clock
  # model comparison while accomodating phylogenetic uncertainty
  # Equation 10: AICM = k - 2m
  #   k = effective number of parameters = 2 * variance of the posterior log likehoods
  #   m = mean of the posterior log likelihoods
  aicm <- (2 * var(logLikelihoodValues)) - (2 * mean(logLikelihoodValues))
 
  # Conduct some bootstrapping
  bootstrapAICMs <- c()
  for(i in 1:nBootstraps){
    
    sample <- sample(logLikelihoodValues, size=length(logLikelihoodValues), replace=TRUE)
    bootstrapAICMs[i] <- (2 * var(sample)) - (2 * mean(sample))
  }
  
  # Summarise the bootstraps
  quantiles <- quantile(bootstrapAICMs, probs=c(0.025, 0.975))
  output <- c(aicm, quantiles[[1]], quantiles[[2]])
  
  return(output)
}

getValues <- function(list){
  values <- c()
  keys <- names(list)
  for(i in 1:length(keys)){
    values[i] <- list[[keys[i]]]
  }
  
  return(values)
}

setRangeOfValues <- function(values, minToSet, maxToSet){
  output <- c()
  range <- range(values)
  for(i in 1:length(values)){
    output[i] <- ((values[i] - range[1]) / (range[2] - range[1])) * (maxToSet - minToSet)
    output[i] <- output[i] + minToSet
  }
  
  return(output)
}

timesValuesInListByValue <- function(list, value){
  
  output <- list()
  for(key in names(list)){
    output[[key]] <- list[[key]] * value
  }
  
  return(output)
}

divideValuesInListByMax <- function(arrowRates){
    max <- max(getValues(arrowRates))
    output <- timesValuesInListByValue(arrowRates, 1/max)
    
  return(output)
}

setRangeOfValuesInList <- function(arrowRates, minToSet, maxToSet){
  
  range <- range(getValues(arrowRates))
  for(key in names(arrowRates)){
    arrowRates[[key]] <- ((arrowRates[[key]] - range[1]) / (range[2] - range[1]) * maxToSet) + minToSet
  }
  
  return(arrowRates)
}

getArrowRates <- function(logTable){
  
  # Initialise a list to store the arrow weights
  arrowRates <- list()
  
  # Get the column names
  colNames <- colnames(logTable)
  
  # Examine each column
  for(col in colNames){
  
    # Ignore all columns except the forward rate columns
    if(grepl(x=col, pattern="forward") == FALSE){
      next
    }
    
    # Get the directional information (i.e. 0_1, 2_0)
    parts <- strsplit(col, split="_")[[1]]
    direction <- paste(parts[length(parts) - 1], "_", parts[length(parts)], sep="")
    
    # Store a summary statistic for the rate distribution
    arrowRates[[direction]] <- median(logTable[, col])
  }
  
  return(arrowRates)
}

plotMigrationRates <- function(logTable, demeStructure, code, arrowFactor){
  
  # Get the migration rates
  output <- getArrowRates(logTable)
  
  # Normalise those rates to vary between 0 and MAX (arrow factor)
  migrationRates <- divideValuesInListByMax(output)
  migrationRates <- timesValuesInListByValue(migrationRates, arrowFactor)
  
  # Plot the demes and associated rates
  if(demeStructure == "2Deme"){
    plot2Deme(migrationRates, code)
  }else if(demeStructure == "3Deme-outerIsBoth"){
    plot3Deme(migrationRates, code, demeStructure)
  }else if(demeStructure == "3Deme-outerIsCattle"){
    plot3Deme(migrationRates, code, demeStructure)
  }else if(demeStructure == "4Deme"){
    plot4Deme(migrationRates, code)
  }else if(demeStructure == "6Deme-EastWest"){
    plot6DemeEW(migrationRates, code)
  }else if(demeStructure == "6Deme-NorthSouth"){
    plot6DemeNS(migrationRates, code)
  }else if(demeStructure == "8Deme-EastWest"){
    plot8DemeEW(migrationRates, code)
  }else if(demeStructure == "8Deme-NorthSouth"){
    plot8DemeNS(migrationRates, code)
  }else{
    cat(paste("Input deme structure not recognised:", demeStructure, "\n"))
  }
  
  return(output)
}

plot2Deme <- function(arrowWeights, code){

  demeStructure <- "2Deme"
  
  # Create empty plot
  createEmptyPlot()
  
  # Get deme names and assign colours
  demeNames <- getDemeNamesForDemeStructure(demeStructure)
  demeColours <- rep("black", length(demeNames))
  demeColours[grepl(demeNames, pattern="badger")] <- "red"
  demeColours[grepl(demeNames, pattern="cow")] <- "blue"
  
  # Add labels
  x <- c(0.1, 0.9)
  y <- c(0.5, 0.5)
  text(x=x, y=y, 
       labels=demeNames,
       col=demeColours)
  
  # badger -> cow
  if(arrowWeights[["0_1"]] != 0){
    arrows(x0=x[1]+0.15, x1=x[2]-0.15, y0=y[1]-0.15, y1=y[2]-0.15,
           code=code, lwd=arrowWeights[["0_1"]])
  }
  # cow -> badger
  if(arrowWeights[["1_0"]] != 0){
    arrows(x0=x[2]-0.15, x1=x[1]+0.15, y0=y[2]+0.15, y1=y[1]+0.15,
         code=code, lwd=arrowWeights[["1_0"]])
  }
}

plot3Deme <- function(arrowWeights, code, demeStructure){

  # Create empty plot
  createEmptyPlot()
  
  # Get deme names and assign colours
  demeNames <- getDemeNamesForDemeStructure(demeStructure)
  demeColours <- rep("black", length(demeNames))
  demeColours[grepl(demeNames, pattern="badger")] <- "red"
  demeColours[grepl(demeNames, pattern="cow")] <- "blue"
  
  # Position the labels
  x <- c(0.1, 0.9, 0.5)
  y <- c(0.1, 0.1, 0.9)
  text(x=x, y=y, 
       labels=demeNames, col=demeColours)
  
  # badger -> cow/cow-inner
  if(arrowWeights[["0_1"]] != 0){
    arrows(x0=x[1]+0.15, x1=x[2]-0.15, y0=y[1]-0.05, y1=y[2]-0.05,
           code=code, lwd=arrowWeights[["0_1"]])
  }
  # cow/cow-inner -> badger
  if(arrowWeights[["1_0"]] != 0){
    arrows(x0=x[2]-0.15, x1=x[1]+0.15, y0=y[2]+0.05, y1=y[1]+0.05,
           code=code, lwd=arrowWeights[["1_0"]])
  }
  
  # badger -> outer/cow-outer
  if(arrowWeights[["0_2"]] != 0){
    arrows(x0=x[1]-0.05, x1=x[3]-0.15, y0=y[1]+0.1, y1=y[3]-0.1,
           code=code, lwd=arrowWeights[["0_2"]])
  }
  # outer/cow-outer -> badger
  if(arrowWeights[["2_0"]] != 0){
    arrows(x0=x[3]-0.05, x1=x[1]+0.05, y0=y[3]-0.1, y1=y[1]+0.1,
           code=code, lwd=arrowWeights[["2_0"]])
  }
  
  # cow/cow-inner -> outer/cow-outer
  if(arrowWeights[["1_2"]] != 0){
    arrows(x0=x[2]+0.05, x1=x[3]+0.15, y0=y[2]+0.1, y1=y[3]-0.1,
           code=code, lwd=arrowWeights[["1_2"]])
  }
  # outer/cow-outer -> cow/cow-inner
  if(arrowWeights[["2_1"]] != 0){
    arrows(x0=x[3]+0.05, x1=x[2]-0.05, y0=y[3]-0.1, y1=y[2]+0.1,
           code=code, lwd=arrowWeights[["2_1"]])
  }
}

plot4Deme <- function(arrowWeights, code){
  
  demeStructure <- "4Deme"
  
  # Create empty plot
  createEmptyPlot()
  
  # Get deme names and assign colours
  demeNames <- getDemeNamesForDemeStructure(demeStructure)
  demeColours <- rep("black", length(demeNames))
  demeColours[grepl(demeNames, pattern="badger")] <- "red"
  demeColours[grepl(demeNames, pattern="cow")] <- "blue"
  
  # Note the x and y positions
  x <- c(0.1, 0.1, 0.9, 0.9)
  y <- c(0.1, 0.9, 0.9, 0.1)
  
  # Add labels
  text(x=x, y=y, 
       labels=demeNames, col=demeColours)

  # badger-inner -> cow-inner
  if(arrowWeights[["0_1"]] != 0){
    arrows(x0=x[1]-0.05, x1=x[2]-0.05, y0=y[1]+0.1, y1=y[2]-0.1,
          code=code, lwd=arrowWeights[["0_1"]])
  }
  # cow-inner -> badger-inner
  if(arrowWeights[["0_1"]] != 0){
    arrows(x0=x[2]+0.05, x1=x[1]+0.05, y0=y[2]-0.1, y1=y[1]+0.1,
           code=code, lwd=arrowWeights[["0_1"]])
  }
  
  # badger-outer -> cow-outer
  if(arrowWeights[["3_2"]] != 0){
    arrows(x0=x[4]+0.05, x1=x[3]+0.05, y0=y[4]+0.1, y1=y[3]-0.1,
           code=code, lwd=arrowWeights[["3_2"]])
  }
  # cow-outer -> badger-outer
  if(arrowWeights[["2_3"]] != 0){
    arrows(x0=x[3]-0.05, x1=x[4]-0.05, y0=y[3]-0.1, y1=y[4]+0.1,
           code=code, lwd=arrowWeights[["2_3"]])
  }
  
  # badger-inner -> badger-outer
  if(arrowWeights[["0_3"]] != 0){
    arrows(x0=x[1]+0.15, x1=x[4]-0.15, y0=y[1]-0.05, y1=y[4]-0.05,
           code=code, lwd=arrowWeights[["0_3"]])
  }
  # badger-outer -> badger-inner
  if(arrowWeights[["3_0"]] != 0){
    arrows(x0=x[4]-0.15, x1=x[1]+0.15, y0=y[4]+0.05, y1=y[1]+0.05,
           code=code, lwd=arrowWeights[["3_0"]])
  }
  
  # cow-inner -> cow-outer
  if(arrowWeights[["1_2"]] != 0){
    arrows(x0=x[2]+0.15, x1=x[3]-0.15, y0=y[2]+0.05, y1=y[3]+0.05,
           code=code, lwd=arrowWeights[["1_2"]])
  }
  # cow-outer -> cow-inner
  if(arrowWeights[["2_1"]] != 0){
    arrows(x0=x[3]-0.15, x1=x[2]+0.15, y0=y[3]-0.05, y1=y[2]-0.05,
           code=code, lwd=arrowWeights[["2_1"]])
  }
}

plot6DemeEW <- function(arrowWeights, code){
  
  demeStructure <- "6Deme-EastWest"
  
  # Create empty plot
  createEmptyPlot()
  
  # Get deme names and assign colours
  demeNames <- getDemeNamesForDemeStructure(demeStructure)
  demeColours <- rep("black", length(demeNames))
  demeColours[grepl(demeNames, pattern="badger")] <- "red"
  demeColours[grepl(demeNames, pattern="cow")] <- "blue"
  
  # Add labels
  x <- c(0.5, 0.5, 0.9, 0.1, 0.9, 0.1)
  y <- c(0.35, 0.65, 0.9, 0.9, 0.1, 0.1)
  text(x=x, y=y, 
       labels=demeNames, col=demeColours)
  
  # badger-inner -> cow-inner
  if(arrowWeights[["0_1"]] != 0){
    arrows(x0=x[1]-0.05, x1=x[2]-0.05, y0=y[1]+0.05, y1=y[2]-0.05,
           code=code, lwd=arrowWeights[["0_1"]])
  }
  # cow-inner -> badger-inner
  if(arrowWeights[["1_0"]] != 0){
    arrows(x0=x[2]+0.05, x1=x[1]+0.05, y0=y[2]-0.05, y1=y[1]+0.05,
           code=code, lwd=arrowWeights[["1_0"]])
  }
  
  # cow-outer-east -> badger-outer-east
  if(arrowWeights[["2_4"]] != 0){
    arrows(x0=x[3]-0.05, x1=x[5]-0.05, y0=y[3]-0.1, y1=y[5]+0.1,
            code=code, lwd=arrowWeights[["2_4"]])
  }
  # badger-outer-east -> cow-outer-east
  if(arrowWeights[["4_2"]] != 0){
    arrows(x0=x[5]+0.05, x1=x[3]+0.05, y0=y[5]+0.1, y1=y[3]-0.1,
           code=code, lwd=arrowWeights[["4_2"]])
  }
  
  # cow-outer-west -> badger-outer-west
  if(arrowWeights[["3_5"]] != 0){
    arrows(x0=x[4]+0.05, x1=x[6]+0.05, y0=y[4]-0.1, y1=y[6]+0.1,
           code=code, lwd=arrowWeights[["3_5"]])
  }
  # badger-outer-west -> cow-outer-west
  if(arrowWeights[["5_3"]] != 0){
    arrows(x0=x[6]-0.05, x1=x[4]-0.05, y0=y[6]+0.1, y1=y[4]-0.1,
           code=code, lwd=arrowWeights[["5_3"]])
  }
  
  # cow-outer-east -> cow-outer-west
  if(arrowWeights[["2_3"]] != 0){
    arrows(x0=x[3]-0.2, x1=x[4]+0.2, y0=y[3]-0.05, y1=y[4]-0.05,
           code=code, lwd=arrowWeights[["2_3"]])
  }
  # cow-outer-west -> cow-outer-east
  if(arrowWeights[["3_2"]] != 0){
    arrows(x0=x[4]+0.2, x1=x[3]-0.2, y0=y[4]+0.05, y1=y[3]+0.05,
           code=code, lwd=arrowWeights[["3_2"]])
  }
  
  # badger-outer-east -> badger-outer-west
  if(arrowWeights[["4_5"]] != 0){
    arrows(x0=x[5]-0.2, x1=x[6]+0.2, y0=y[5]+0.05, y1=y[6]+0.05,
           code=code, lwd=arrowWeights[["4_5"]])
  }
  # badger-outer-west -> badger-outer-east
  if(arrowWeights[["5_4"]] != 0){
    arrows(x0=x[6]+0.2, x1=x[5]-0.2, y0=y[6]-0.05, y1=y[5]-0.05,
           code=code, lwd=arrowWeights[["5_4"]])
  }
  
  # badger-inner -> badger-outer-east
  if(arrowWeights[["0_4"]] != 0){
    arrows(x0=x[1]+0.2, x1=x[5]-0.1, y0=y[1]-0.05, y1=y[5]+0.1,
           code=code, lwd=arrowWeights[["0_4"]])
  }
  # badger-outer-east -> badger-inner
  if(arrowWeights[["4_0"]] != 0){
    arrows(x0=x[5]-0.2, x1=x[1]+0.1, y0=y[5]+0.1, y1=y[1]-0.05,
           code=code, lwd=arrowWeights[["4_0"]])
  }
 
  # badger-inner -> badger-outer-west
  if(arrowWeights[["0_5"]] != 0){
    arrows(x0=x[1]-0.2, x1=x[6]+0.1, y0=y[1]-0.05, y1=y[6]+0.1,
           code=code, lwd=arrowWeights[["0_5"]])
  }
  # badger-outer-west -> badger-inner
  if(arrowWeights[["5_0"]] != 0){
    arrows(x0=x[6]+0.2, x1=x[1]-0.1, y0=y[6]+0.1, y1=y[1]-0.05,
           code=code, lwd=arrowWeights[["5_0"]])
  }
  
  # cow-inner -> cow-outer-east
  if(arrowWeights[["1_2"]] != 0){
    arrows(x0=x[2]+0.2, x1=x[3]-0.1, y0=y[2]+0.05, y1=y[3]-0.1,
           code=code, lwd=arrowWeights[["1_2"]])
  }
  # cow-outer-east -> cow-inner
  if(arrowWeights[["2_1"]] != 0){
    arrows(x0=x[3]-0.2, x1=x[2]+0.1, y0=y[3]-0.1, y1=y[2]+0.05,
           code=code, lwd=arrowWeights[["2_1"]])
  }

  # cow-inner -> cow-outer-west
  if(arrowWeights[["1_3"]] != 0){
    arrows(x0=x[2]-0.2, x1=x[4]+0.1, y0=y[2]+0.05, y1=y[4]-0.1,
           code=code, lwd=arrowWeights[["1_3"]])
  }
  # cow-outer-west -> cow-inner
  if(arrowWeights[["3_1"]] != 0){
    arrows(x0=x[4]+0.2, x1=x[2]-0.1, y0=y[4]-0.1, y1=y[2]+0.05,
           code=code, lwd=arrowWeights[["3_1"]])
  }
}

plot6DemeNS <- function(arrowWeights, code){
  
  demeStructure <- "6Deme-NorthSouth"
  
  # Create empty plot
  createEmptyPlot()
  
  # Get deme names and assign colours
  demeNames <- getDemeNamesForDemeStructure(demeStructure)
  demeColours <- rep("black", length(demeNames))
  demeColours[grepl(demeNames, pattern="badger")] <- "red"
  demeColours[grepl(demeNames, pattern="cow")] <- "blue"
  
  # Add labels
  x <- c(0.35, 0.7, 0.9, 0.9, 0.1, 0.1)
  y <- c(0.5, 0.5, 0.9, 0.1, 0.9, 0.1)
  text(x=x, y=y, 
       labels=demeNames, col=demeColours)
  
  # badger-inner -> cow-inner
  if(arrowWeights[["0_1"]] != 0){
    arrows(x0=x[1]+0.1, x1=x[2]-0.1, y0=y[1]-0.05, y1=y[2]-0.05,
           code=code, lwd=arrowWeights[["0_1"]])
  }
  # cow-inner -> badger-inner
  if(arrowWeights[["1_0"]] != 0){
      arrows(x0=x[2]-0.1, x1=x[1]+0.1, y0=y[2]+0.05, y1=y[1]+0.05,
             code=code, lwd=arrowWeights[["1_0"]])
  }

  # cow-outer-north -> badger-outer-north
  if(arrowWeights[["2_4"]] != 0){
    arrows(x0=x[3]-0.2, x1=x[5]+0.2, y0=y[3]-0.05, y1=y[5]-0.05,
           code=code, lwd=arrowWeights[["2_4"]])
  }
  # cow-outer-north -> badger-outer-north
  if(arrowWeights[["4_2"]] != 0){
    arrows(x0=x[5]+0.2, x1=x[3]-0.2, y0=y[5]+0.05, y1=y[3]+0.05,
           code=code, lwd=arrowWeights[["4_2"]])
  }
  
  # cow-outer-south -> badger-outer-south
  if(arrowWeights[["3_5"]] != 0){
    arrows(x0=x[4]-0.2, x1=x[6]+0.2, y0=y[4]+0.05, y1=y[6]+0.05,
           code=code, lwd=arrowWeights[["3_5"]]) 
  }
  # badger-outer-south -> cow-outer-south 
  if(arrowWeights[["5_3"]] != 0){
    arrows(x0=x[6]+0.2, x1=x[4]-0.2, y0=y[6]-0.05, y1=y[4]-0.05,
           code=code, lwd=arrowWeights[["5_3"]])
  }
  
  # badger-outer-north -> badger-outer-south
  if(arrowWeights[["4_5"]] != 0){
    arrows(x0=x[5]+0.05, x1=x[6]+0.05, y0=y[5]-0.1, y1=y[6]+0.1,
           code=code, lwd=arrowWeights[["4_5"]]) 
  }
  # badger-outer-south -> badger-outer-north
  if(arrowWeights[["5_4"]] != 0){
    arrows(x0=x[6]-0.05, x1=x[5]-0.05, y0=y[6]+0.1, y1=y[5]-0.1,
           code=code, lwd=arrowWeights[["5_4"]])
  } 
  
  # cow-outer-north -> cow-outer-south
  if(arrowWeights[["2_4"]] != 0){
    arrows(x0=x[3]-0.05, x1=x[4]-0.05, y0=y[3]-0.1, y1=y[4]+0.1,
           code=code, lwd=arrowWeights[["2_4"]]) 
  }
  # cow-outer-south -> cow-outer-north
  if(arrowWeights[["3_2"]] != 0){
    arrows(x0=x[4]+0.05, x1=x[3]+0.05, y0=y[4]+0.1, y1=y[3]-0.1,
           code=code, lwd=arrowWeights[["3_2"]])
  }
  
  # badger-inner -> badger-outer-north
  if(arrowWeights[["0_4"]] != 0){
    arrows(x0=x[1]-0.1, x1=x[5]+0.1, y0=y[1]+0.05, y1=y[5]-0.1,
           code=code, lwd=arrowWeights[["0_4"]])
  }
  # badger-outer-north -> badger-inner
  if(arrowWeights[["4_0"]] != 0){
    arrows(x0=x[5]+0.2, x1=x[1], y0=y[5]-0.1, y1=y[1]+0.05,
           code=code, lwd=arrowWeights[["4_0"]])
  }
  
  # badger-inner -> badger-outer-south
  if(arrowWeights[["0_5"]] != 0){
    arrows(x0=x[1]-0.1, x1=x[6]+0.1, y0=y[1]-0.05, y1=y[6]+0.1,
           code=code, lwd=arrowWeights[["0_5"]])
  }
  # badger-outer-south -> badger-inner
  if(arrowWeights[["5_0"]] != 0){
    arrows(x0=x[6]+0.2, x1=x[1], y0=y[6]+0.1, y1=y[1]-0.05,
           code=code, lwd=arrowWeights[["5_0"]])
  }
  
  # cow-inner -> cow-outer-north
  if(arrowWeights[["1_2"]] != 0){
    arrows(x0=x[2]+0.05, x1=x[3]-0.1, y0=y[2]+0.05, y1=y[3]-0.1,
           code=code, lwd=arrowWeights[["1_2"]])
  }
  # cow-outer-north -> cow-inner
  if(arrowWeights[["2_1"]] != 0){
    arrows(x0=x[3]-0.2, x1=x[2]-0.05, y0=y[3]-0.1, y1=y[2]+0.05,
           code=code, lwd=arrowWeights[["2_1"]])
  }
  
  # cow-inner -> cow-outer-south
  if(arrowWeights[["1_3"]] != 0){
    arrows(x0=x[2]+0.05, x1=x[4]-0.1, y0=y[2]-0.05, y1=y[4]+0.1,
           code=code, lwd=arrowWeights[["1_3"]])
  }
  # cow-outer-south -> cow-inner
  if(arrowWeights[["3_1"]] != 0){
    arrows(x0=x[4]-0.2, x1=x[2]-0.05, y0=y[4]+0.1, y1=y[2]-0.05,
           code=code, lwd=arrowWeights[["3_1"]])
  }
}

plot8DemeEW <- function(arrowWeights, code){
  
  demeStructure <- "8Deme-EastWest"
  
  # Create empty plot
  createEmptyPlot()
  
  # Get deme names and assign colours
  demeNames <- getDemeNamesForDemeStructure(demeStructure)
  demeColours <- rep("black", length(demeNames))
  demeColours[grepl(demeNames, pattern="badger")] <- "red"
  demeColours[grepl(demeNames, pattern="cow")] <- "blue"
  
  # Add labels
  x <- c(0.25,  0.25,  0.75,  0.75,  0.9, 0.9, 0.1, 0.1)
  y <- c(0.65,  0.35,  0.65,  0.35,  0.9, 0.1, 0.9, 0.1)
  
  x <- c(0.75,  0.25,  0.75,  0.25,  0.9, 0.1, 0.9, 0.1)
  y <- c(0.35,  0.35,  0.65,  0.65,  0.9, 0.9, 0.1, 0.1)
  
  text(x=x, y=y, labels=demeNames, col=demeColours)

  # cow-inner-east -> cow-inner-west
  if(arrowWeights[["2_3"]] != 0){
    arrows(x0=x[3]-0.175, x1=x[4]+0.175, y0=y[3]-0.05, y1=y[4]-0.05,
           code=code, lwd=arrowWeights[["2_3"]])
  }
  # cow-inner-west -> cow-inner-east
  if(arrowWeights[["3_2"]] != 0){
    arrows(x0=x[4]+0.175, x1=x[3]-0.175, y0=y[4]+0.05, y1=y[3]+0.05,
           code=code, lwd=arrowWeights[["3_2"]])
  }
  
  # badger-inner-east -> badger-inner-west
  if(arrowWeights[["0_1"]] != 0){
    arrows(x0=x[1]-0.175, x1=x[2]+0.175, y0=y[1]+0.05, y1=y[2]+0.05,
           code=code, lwd=arrowWeights[["0_1"]])
  }
  # badger-inner-west -> badger-inner-east
  if(arrowWeights[["1_0"]] != 0){
    arrows(x0=x[2]+0.175, x1=x[1]-0.175, y0=y[2]-0.05, y1=y[1]-0.05,
           code=code, lwd=arrowWeights[["1_0"]])
  }
    
  # badger-inner-west -> cow-inner-west
  if(arrowWeights[["1_3"]] != 0){
    arrows(x0=x[2]-0.05, x1=x[4]-0.05, y0=y[2]+0.05, y1=y[4]-0.05,
           code=code, lwd=arrowWeights[["1_3"]])
  }
  # cow-inner-west -> badger-inner-west
  if(arrowWeights[["3_1"]] != 0){
    arrows(x0=x[4]+0.05, x1=x[2]+0.05, y0=y[4]-0.05, y1=y[2]+0.05,
           code=code, lwd=arrowWeights[["3_1"]])
  }

  # badger-inner-east -> cow-inner-east
  if(arrowWeights[["0_2"]] != 0){
    arrows(x0=x[1]+0.05, x1=x[3]+0.05, y0=y[1]+0.05, y1=y[3]-0.05,
           code=code, lwd=arrowWeights[["0_2"]])
  }
  # cow-inner-east -> badger-inner-east
  if(arrowWeights[["2_0"]] != 0){
    arrows(x0=x[3]-0.05, x1=x[1]-0.05, y0=y[3]-0.05, y1=y[1]+0.05,
           code=code, lwd=arrowWeights[["2_0"]])
  }
  
  # cow-outer-east -> badger-outer-east
  if(arrowWeights[["4_6"]] != 0){
    arrows(x0=x[5]+0.025, x1=x[7]+0.025, y0=y[5]-0.05, y1=y[7]+0.05,
           code=code, lwd=arrowWeights[["4_6"]])
  }
  # badger-outer-east -> cow-outer-east
  if(arrowWeights[["6_4"]] != 0){
    arrows(x0=x[7]+0.1, x1=x[5]+0.1, y0=y[7]+0.05, y1=y[5]-0.05,
           code=code, lwd=arrowWeights[["6_4"]])
  }
  
  # cow-outer-west -> badger-outer-west
  if(arrowWeights[["5_7"]] != 0){
    arrows(x0=x[6]-0.025, x1=x[8]-0.025, y0=y[6]-0.05, y1=y[8]+0.05,
           code=code, lwd=arrowWeights[["5_7"]])
  }
  # badger-outer-west -> cow-outer-west
  if(arrowWeights[["7_5"]] != 0){
    arrows(x0=x[8]-0.1, x1=x[6]-0.1, y0=y[8]+0.05, y1=y[6]-0.05,
           code=code, lwd=arrowWeights[["7_5"]])
  }
 
  # cow-outer-east -> cow-outer-west
  if(arrowWeights[["4_5"]] != 0){
    arrows(x0=x[5]-0.175, x1=x[6]+0.175, y0=y[5]-0.05, y1=y[6]-0.05,
           code=code, lwd=arrowWeights[["4_5"]])
  }
  # cow-outer-west -> cow-outer-east
  if(arrowWeights[["5_4"]] != 0){
    arrows(x0=x[6]+0.175, x1=x[5]-0.175, y0=y[6]+0.05, y1=y[5]+0.05,
           code=code, lwd=arrowWeights[["5_4"]])
  }
  
  # badger-outer-east -> badger-outer-west
  if(arrowWeights[["6_7"]] != 0){
    arrows(x0=x[7]-0.175, x1=x[8]+0.175, y0=y[7]+0.05, y1=y[8]+0.05,
           code=code, lwd=arrowWeights[["6_7"]])
  }
  # badger-outer-west -> badger-outer-east
  if(arrowWeights[["7_6"]] != 0){
    arrows(x0=x[8]+0.175, x1=x[7]-0.175, y0=y[8]-0.05, y1=y[7]-0.05,
           code=code, lwd=arrowWeights[["7_6"]])
  }

  # badger-inner-west -> badger-outer-west
  if(arrowWeights[["1_7"]] != 0){
    arrows(x0=x[2]-0.05, x1=x[8]+0.025, y0=y[2]-0.05, y1=y[8]+0.05,
           code=code, lwd=arrowWeights[["1_7"]])
  }
  # badger-outer-west -> badger-inner-west
  if(arrowWeights[["7_1"]] != 0){
    arrows(x0=x[8]+0.125, x1=x[2]+0.05, y0=y[8]+0.05, y1=y[2]-0.05,
           code=code, lwd=arrowWeights[["7_1"]])
  }
  # cow-inner-west -> cow-outer-west
  if(arrowWeights[["3_5"]] != 0){
    arrows(x0=x[4]-0.05, x1=x[6]+0.025, y0=y[4]+0.05, y1=y[6]-0.05,
           code=code, lwd=arrowWeights[["3_5"]])
  }
  # cow-outer-west -> cow-inner-west
  if(arrowWeights[["5_3"]] != 0){
    arrows(x0=x[6]+0.125, x1=x[4]+0.05, y0=y[6]-0.05, y1=y[4]+0.05,
           code=code, lwd=arrowWeights[["5_3"]])
  }
  
  # cow-inner-east -> cow-outer-east
  if(arrowWeights[["2_4"]] != 0){
    arrows(x0=x[3]+0.05, x1=x[5], y0=y[3]+0.05, y1=y[5]-0.05,
           code=code, lwd=arrowWeights[["2_4"]])
  }
  # cow-outer-east -> cow-inner-east
  if(arrowWeights[["4_2"]] != 0){
    arrows(x0=x[5]-0.1, x1=x[3]-0.05, y0=y[5]-0.05, y1=y[3]+0.05,
           code=code, lwd=arrowWeights[["4_2"]])
  }
  
  # badger-inner-east -> badger-outer-east
  if(arrowWeights[["0_6"]] != 0){
    arrows(x0=x[1]+0.05, x1=x[7], y0=y[1]-0.05, y1=y[7]+0.05,
           code=code, lwd=arrowWeights[["0_6"]])
  }
  # badger-outer-east -> badger-inner-east
  if(arrowWeights[["6_0"]] != 0){
    arrows(x0=x[7]-0.1, x1=x[1]-0.05, y0=y[7]+0.05, y1=y[1]-0.05,
           code=code, lwd=arrowWeights[["6_0"]])
  }
}

plot8DemeNS <- function(arrowWeights, code){
  
  demeStructure <- "8Deme-NorthSouth"
  
  # Create empty plot
  createEmptyPlot()
  
  # Get deme names and assign colours
  demeNames <- getDemeNamesForDemeStructure(demeStructure)
  demeColours <- rep("black", length(demeNames))
  demeColours[grepl(demeNames, pattern="badger")] <- "red"
  demeColours[grepl(demeNames, pattern="cow")] <- "blue"
  
  # Add labels
#  x <- c(0.25, 0.25, 0.75, 0.75, 0.9,  0.9,  0.1,  0.1)
#  y <- c(0.9,  0.1,  0.9,  0.1,  0.65, 0.35, 0.65, 0.35)

#  x <- c(0.1, 0.1, 0.9, 0.9, 0.75,  0.75,  0.25,  0.25)
#  y <- c(0.65,  0.35,  0.65,  0.35,  0.9, 0.1, 0.9, 0.1)
  
  x <- c(0.25, 0.25, 0.75, 0.75, 0.9,  0.9,  0.1,  0.1)
  y <- c(0.65,  0.35,  0.65,  0.35,  0.9, 0.1, 0.9, 0.1)
  
  text(x=x, y=y, labels=demeNames, col=demeColours)
  
  # badger-inner-north -> badger-inner-south
  if(arrowWeights[["0_1"]] != 0){
    arrows(x0=x[1]+0.05, x1=x[2]+0.05, y0=y[1]-0.05, y1=y[2]+0.05,
           code=code, lwd=arrowWeights[["0_1"]])
  }
  # badger-inner-south -> badger-inner-north
  if(arrowWeights[["1_0"]] != 0){
    arrows(x0=x[2]-0.05, x1=x[1]-0.05, y0=y[2]+0.05, y1=y[1]-0.05,
           code=code, lwd=arrowWeights[["1_0"]])
  }
  
  # cow-inner-north -> cow-inner-south
  if(arrowWeights[["2_3"]] != 0){
    arrows(x0=x[3]-0.05, x1=x[4]-0.05, y0=y[3]-0.05, y1=y[4]+0.05,
           code=code, lwd=arrowWeights[["2_3"]])
  }
  # cow-inner-south -> cow-inner-north
  if(arrowWeights[["3_2"]] != 0){
    arrows(x0=x[4]+0.05, x1=x[3]+0.05, y0=y[4]+0.05, y1=y[3]-0.05,
           code=code, lwd=arrowWeights[["3_2"]])
  }
  
  # badger-inner-north -> cow-inner-north
  if(arrowWeights[["0_2"]] != 0){
    arrows(x0=x[1]+0.175, x1=x[3]-0.15, y0=y[1]+0.05, y1=y[3]+0.05,
           code=code, lwd=arrowWeights[["0_2"]])
  }
  # cow-inner-north -> badger-inner-north
  if(arrowWeights[["2_0"]] != 0){
    arrows(x0=x[3]-0.15, x1=x[1]+0.175, y0=y[3]-0.05, y1=y[1]-0.05,
           code=code, lwd=arrowWeights[["2_0"]])
  }
  
  # badger-inner-south -> cow-inner-south
  if(arrowWeights[["1_3"]] != 0){
    arrows(x0=x[2]+0.175, x1=x[4]-0.15, y0=y[2]-0.05, y1=y[4]-0.05,
           code=code, lwd=arrowWeights[["1_3"]])
  }
  # cow-inner-south -> badger-inner-south
  if(arrowWeights[["3_1"]] != 0){
    arrows(x0=x[4]-0.15, x1=x[2]+0.175, y0=y[4]+0.05, y1=y[2]+0.05,
           code=code, lwd=arrowWeights[["3_1"]])
  }
  
  # cow-outer-south -> badger-outer-south
  if(arrowWeights[["5_7"]] != 0){
    arrows(x0=x[6]-0.175, x1=x[8]+0.2, y0=y[6]+0.05, y1=y[8]+0.05,
           code=code, lwd=arrowWeights[["5_7"]])
  }
  # badger-outer-south -> cow-outer-south
  if(arrowWeights[["7_5"]] != 0){
    arrows(x0=x[8]+0.2, x1=x[6]-0.175, y0=y[8]-0.05, y1=y[6]-0.05,
           code=code, lwd=arrowWeights[["7_5"]])
  }
  
  # cow-outer-north -> badger-outer-north
  if(arrowWeights[["4_6"]] != 0){
    arrows(x0=x[5]-0.175, x1=x[7]+0.2, y0=y[5]-0.05, y1=y[7]-0.05,
           code=code, lwd=arrowWeights[["4_6"]])
  }
  # badger-outer-north -> cow-outer-north
  if(arrowWeights[["6_4"]] != 0){
    arrows(x0=x[7]+0.2, x1=x[5]-0.175, y0=y[7]+0.05, y1=y[5]+0.05,
           code=code, lwd=arrowWeights[["6_4"]])
  }
  
  # badger-outer-north -> badger-outer-south
  if(arrowWeights[["6_7"]] != 0){
    arrows(x0=x[7]-0.025, x1=x[8]-0.025, y0=y[7]-0.05, y1=y[8]+0.05,
           code=code, lwd=arrowWeights[["6_7"]])
  }
  # badger-outer-south -> badger-outer-north
  if(arrowWeights[["7_6"]] != 0){
    arrows(x0=x[8]-0.1, x1=x[7]-0.1, y0=y[8]+0.05, y1=y[7]-0.05,
           code=code, lwd=arrowWeights[["7_6"]])
  }
  
  # cow-outer-north -> cow-outer-south
  if(arrowWeights[["5_6"]] != 0){
    arrows(x0=x[5], x1=x[6], y0=y[5]-0.05, y1=y[6]+0.05,
           code=code, lwd=arrowWeights[["5_6"]])
  }
  # cow-outer-south -> cow-outer-north
  if(arrowWeights[["6_5"]] != 0){
    arrows(x0=x[6]+0.075, x1=x[5]+0.075, y0=y[6]+0.05, y1=y[5]-0.05,
           code=code, lwd=arrowWeights[["6_5"]])
  }
  
  # badger-inner-north -> badger-outer-north
  if(arrowWeights[["0_6"]] != 0){
    arrows(x0=x[1]-0.05, x1=x[7]+0.025, y0=y[1]+0.05, y1=y[7]-0.05,
           code=code, lwd=arrowWeights[["0_6"]])
  }
  # badger-outer-north -> badger-inner-north
  if(arrowWeights[["6_0"]] != 0){
    arrows(x0=x[7]+0.125, x1=x[1]+0.05, y0=y[7]-0.05, y1=y[1]+0.05,
           code=code, lwd=arrowWeights[["6_0"]])
  }
  
  # badger-inner-south -> badger-outer-south
  if(arrowWeights[["1_7"]] != 0){
    arrows(x0=x[2]-0.05, x1=x[8]+0.025, y0=y[2]-0.05, y1=y[8]+0.05,
           code=code, lwd=arrowWeights[["1_7"]])
  }
  # badger-outer-south -> badger-inner-south
  if(arrowWeights[["7_1"]] != 0){
    arrows(x0=x[8]+0.125, x1=x[2]+0.05, y0=y[8]+0.05, y1=y[2]-0.05,
           code=code, lwd=arrowWeights[["7_1"]])
  }
  
  # cow-inner-north -> cow-outer-north
  if(arrowWeights[["2_4"]] != 0){
    arrows(x0=x[3]+0.05, x1=x[5]-0.025, y0=y[3]+0.05, y1=y[5]-0.05,
           code=code, lwd=arrowWeights[["2_4"]])
  }
  # cow-outer-north -> cow-inner-north
  if(arrowWeights[["4_2"]] != 0){
    arrows(x0=x[5]-0.125, x1=x[3]-0.05, y0=y[5]-0.05, y1=y[3]+0.05,
           code=code, lwd=arrowWeights[["4_2"]])
  }
  
  # cow-inner-south -> cow-outer-south
  if(arrowWeights[["3_5"]] != 0){
    arrows(x0=x[4]+0.05, x1=x[6]-0.025, y0=y[4]-0.05, y1=y[6]+0.05,
           code=code, lwd=arrowWeights[["3_5"]])
  }
  # cow-outer-south -> cow-inner-south
  if(arrowWeights[["5_3"]] != 0){
    arrows(x0=x[6]-0.125, x1=x[4]-0.05, y0=y[6]+0.05, y1=y[4]-0.05,
           code=code, lwd=arrowWeights[["5_3"]])
  }
}

calculateForwardMigrationRates <- function(logTable){
  
  # Get the names of the backward in time migration rate estimates
  migrationRateCols <- colnames(logTable)[
    grepl(colnames(logTable), pattern = "migModel.rateMatrix")]
  
  # For each backward in time migration rate caculate the forward migration rate
  # FMR_ab = BMR_ba * (Nb / Na)
  #   MR: Migration rate (F - Forward, B - Backward)
  #   N: Effective population size
  #   Demes: a, b
  # Equation taken from second paragraph of Methods section in:
  # De Maio et al. 2015 - New routes to phylogeography ...
  # BMR_ba = FMR_ab * (Na / Nb)
  # FMR_ab = BMR_ba / (Na / Nb)
  # FMR_ab = BMR_ba * (Nb / Na)
  
  for(backwardMigrationRateCol in migrationRateCols){
    
    # Get the demes involved
    parts <- strsplit(backwardMigrationRateCol, split="_")[[1]]
    a <- parts[2]
    b <- parts[3]
    
    # Calculate forward rate
    forwardMigrationRateCol <- paste("migModel.forwardRateMatrix_", b, "_", a, sep="")
    logTable[, forwardMigrationRateCol] <-
      logTable[, backwardMigrationRateCol] *
      (logTable[, paste("migModel.popSize_", a, sep="")] * 
         logTable[, paste("migModel.popSize_", b, sep="")])
  }
  
  return(logTable)
}

createEmptyPlot <- function(){
  par(mar=c(0,0,0,0))
  plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", bty="n",
       ylab="", xlab="")
}

plotPosteriorSupportForEachDemeAsRoot <- function(logTable, demeStructure){

  # Reset Margins
  par(mar=c(7, 0, 3, 0))
  
  barplot(table(logTable$treePrior.rootColor), yaxt="n", 
          names=getDemeNamesForDemeStructure(demeStructure),
          main="Assignment of Demes to Root State",
          cex.names=0.5, las=2)
  
  # Reset Margins
  par(mar=c(5.1, 4.1, 4.1, 2.1))}

plotParameterESSValues <- function(logTable, colNamesToPlot){
  
  essValues <- rep(NA, length(colNamesToPlot))

  for(i in 1:length(colNamesToPlot)){

    # Remove NAs, if present
    values <- as.numeric(logTable[, colNamesToPlot[i]])
    values <- values[is.na(values) == FALSE]
    
    # Note and skip those where the value is always the same
    range <- range(values)
    if(range[1] == range[2]){
      #cat(paste("Parameter \"", colNamesToPlot[i],
      #          "\" always has single value: ", range[1], "\n", sep=""))
      next
    }
    
    essValues[i] <- calculateEffectiveSampleSize(values)
  }

  # Set the margins
  par(mar=c(0,11,2,0.5)) # bottom, left, top, right
  
  barplot(essValues, las=1, names=colNamesToPlot, horiz=TRUE, xaxt="n", cex.names=0.3,
          main="Parameter Effective Sample Sizes")
  abline(v=100, lty=2, col="red")
  abline(v=1000, lty=2, col="blue")
  text(x=c(250, 1200), y=c(0, 0), labels=c("100", "1000"), cex=0.5, col=c("red", "blue"))
  
  # Reset Margins
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

calculateEffectiveSampleSize <- function(posteriorSample){
  
  # Calculate the mean of the sample
  sampleMean <- mean(posteriorSample)
  nSamples <- length(posteriorSample)
  
  # Create an array to store the correlation values between a[index] vs. a[index + lag]
  autoCorrelationValues <- c()
  
  for(lag in 1:(nSamples - 1)){
    
    # Calculate the auto-correlation value for the current lag period
    # Taken from: http://www.itl.nist.gov/div898/handbook/eda/section3/autocopl.htm
    a <- posteriorSample[1:(nSamples - lag)] - sampleMean
    b <- posteriorSample[(1:(nSamples - lag) + lag)] - sampleMean
    
    Ch <- (1/nSamples) * sum(a * b)
    C0 <- sum((posteriorSample - sampleMean)^2) / nSamples
    autoCorrelationValues[lag] <- Ch / C0
    
    # Check if auto-correlation has dropped below zero
    if(lag != 1 && 
       (autoCorrelationValues[lag - 1] > 0 && autoCorrelationValues[lag] <= 0 || 
        autoCorrelationValues[lag - 1] < 0 && autoCorrelationValues[lag] >= 0)){
      break;
    }
    
    # Monitor progress
    #if(lag %% 1000 == 0){
    #  cat(paste("Calculate correlation based upon gap of ", lag, ". Max gap = ", (nSamples - 1), "\n", sep=""))
    #}
  }
  
  # Calculate the Effective Sample Size
  # Taken from: http://people.duke.edu/~ccc14/sta-663-2016/16C_PyMC3.html
  ess <- nSamples / (1 + 2 * (sum(autoCorrelationValues)))
  
  return(ess);
}

getDemeNamesForDemeStructure <- function(demeStructure, number=NULL){

  # Species-InOrOut-Location
  # badger  inner   east
  # cow     outer   west
  #                 north
  #                 south
  demeNames <- list(
    "2Deme"=c("badger", "cow"),
    
    "3Deme-outerIsBoth"=c("badger-inner", "cow-inner", "outer"),
    
    "3Deme-outerIsCattle"=c("badger-inner", "cow-inner", "cow"),
    
    "4Deme"=c("badger-inner", "cow-inner", "cow-outer", "badger-outer"),
    
    "6Deme-EastWest"=c("badger-inner", "cow-inner", "cow-outer-east", "cow-outer-west", "badger-outer-east", "badger-outer-west"),
    
    "6Deme-NorthSouth"=c("badger-inner", "cow-inner", "cow-outer-north", "cow-outer-south", "badger-outer-north", "badger-outer-south"),
    
    "8Deme-EastWest"=c("badger-inner-east", "badger-inner-west", "cow-inner-east", "cow-inner-west",
                       "cow-outer-east", "cow-outer-west", "badger-outer-east", "badger-outer-west"),
    
    "8Deme-NorthSouth"=c("badger-inner-north", "badger-inner-south", "cow-inner-north", "cow-inner-south",
                         "cow-outer-north", "cow-outer-south", "badger-outer-north", "badger-outer-south")
  )

  output = demeNames[[demeStructure]]
  if(is.null(number) == FALSE){
    output = demeNames[[demeStructure]][number + 1]
  }
  
  return(output)
}

readInBASTALogTables <- function(demeStructureDates, popEstimationTypes, clockEstimateTypes,
                                 path, burnInProp=0.1){
  
  # Store each of the log tables in a list
  logTables <- list()
  
  for(demeStructure in names(demeStructureDates)){
    
    for(popEstimationType in popEstimationTypes){
      
      for(clockEstimationType in clockEstimateTypes){
        
        # Build run defining prefix
        prefix <- paste(demeStructure, "_", popEstimationType, "_",  clockEstimationType, "_",
                        demeStructureDates[[demeStructure]], sep="")
        
        # Print progress information
        cat(paste("\rReading: ", prefix, ".log                    ", sep=""))
        
        # Create file name
        file <- paste(path, prefix, "/", prefix, ".log", sep="")
        
        # Read in the file as table
        logTable <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
        
        # Replace "N" in table with NAs
        logTable[logTable == "N"] <- NA
        
        # Remove the burn-in
        burnIn <- round(burnInProp * nrow(logTable), digits=0)
        logTable <- logTable[burnIn:nrow(logTable), ]
        
        # Calculate the forward rates
        logTable <- calculateForwardMigrationRates(logTable)
        
        # Store the table
        logTables[[prefix]] <- logTable
      }
    }
  }
  cat("\rFinished reading in log tables...                              ")
  
  return(logTables)
}
