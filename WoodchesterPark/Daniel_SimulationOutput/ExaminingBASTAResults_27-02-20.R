#### Preparation ####

# Load the required libraries
library(ape) # Reading in FASTA
library(lubridate) # Converting to decimal dates
library(basicPlotteR) # Set alpha of colours

# Set the path
path <- file.path("~", "Desktop", "BuildingBASTAXML")

# Set the date
date <- format(Sys.Date(), "%d-%m-%y")

# Note the states/demes considered - ALPHABETICAL ORDER
demes <- c("badger", "cow")

# Note the number of sites in genome examined
genomeSize <- 4345492

#### Read in the BASTA .log file ####

# Note the files in BASTA directory
files <- list.files(file.path(path, "BASTA_equal_relaxed_25-02-20"), full.names=TRUE)

# Read in the log table
logFile <- files[grepl(files, pattern=".log")]
logTable <- readLogFile(logFile, burnInProp=0.1)

#### Count transitions between demes estimated in .trees file ####

# Note the path to the JAVA jar tool that counts transitions in .trees file
pathToTransitionCountingJarFile <- file.path(path, "CountTransitions_07-06-19.jar")

# Note the name of the .trees file
treesFile <- files[grepl(files, pattern=".trees")]

# Count the number of transitions between demes
countsFile <- paste0(substr(treesFile, 1, nchar(treesFile)-6), "_TransitionCounts.txt")
countTransitionsOnPosteriorTrees(countsFile, logFile, treesFile, demes, pathToTransitionCountingJarFile)

#### Create summary plots for log table ####



#### FUNCTIONS - transition counts ####

countTransitionsOnPosteriorTrees <- function(countsFile, logFile, treesFile, demes, pathToJarFile){
  
  # Check that transition counts file doesn't exist
  if(file.exists(countsFile) == FALSE){
    
    # Collapse states into string
    demes <- paste(demes, collapse=",")
    
    # Build the command to run the CountTransitions_DATE.jar file
    command <- paste0("java -jar ", pathToJarFile, " ", treesFile, " ", logFile, " ", demes)
    
    # Count the number of transitions between states for the current posterior tree distribution
    system(command, wait=TRUE)
  }
}

#### FUNCTIONS - log file ####

plotMigrationRates <- function(logTable, demes, code, arrowFactor){
  
  # Get the migration rates
  output <- getArrowRates(logTable)
  
  # Calculate the median for each migration rate being estimate
  medians <- list()
  for(rate in names(output)[names(output) != "AICM"]){
    medians[[rate]] <- median(output[[rate]], na.rm=TRUE)
  }
  
  # Normalise those rates to vary between 0 and MAX (arrow factor)
  migrationRates <- divideValuesInListByMax(medians)
  migrationRates <- timesValuesInListByValue(migrationRates, arrowFactor)
  
  # Check for "NaN" or NA values
  for(key in names(migrationRates)){
    if(is.nan(migrationRates[[key]]) == TRUE || is.na(migrationRates[[key]]) == TRUE){
      migrationRates[[key]] <- 0
    }
  }
  
  # Plot the demes and associated rates
  plotDemesAndMigrationRates(migrationRates, code)
  
  return(output)
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
  colours <- sapply(head(colours), setAlpha, alpha=0.5, simplify=TRUE)
  
  # Plot the substitution rate distribution versus the root height
  plot(x=rateEstimates * genomeSize, y=logTable$tree.height, pch=20, 
       xlab="Substitution Rate (per Genome per Year)", ylab="Root Height (years)",
       main="Substitution Rate versus Root Height",
       col=colours, las=1, bty="n")
  legend("topright", legend=c("High", "Low"), pch=20, col=c("red", "blue"),
         bty="n")
}

plotPopulationSizes <- function(logTable, demeNames, alpha=0.5){
  
  # Get the population size distributions
  popSizeEstimates <- list()
  for(col in colnames(logTable)){
    
    if(grepl(col, pattern="popSize") == FALSE){
      next
    }
    
    popSizeEstimates[[strsplit(col, split="_")[[1]][2]]] <- logTable[, col]
  }
  
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

plotParameterESSValues <- function(logTable, colNamesToPlot){
  
  essValues <- rep(NA, length(colNamesToPlot))
  
  for(i in 1:length(colNamesToPlot)){
    
    # Get the posterior values from the current column
    values <- as.numeric(logTable[, colNamesToPlot[i]])
    
    # Remove transition rate values when rate flag is set to zero
    if(grepl(colNamesToPlot[i], pattern="migModel.rateMatrix") == TRUE){
      
      # Get the deme numbers
      parts = strsplit(colNamesToPlot[i], split="_")[[1]]
      a <- parts[2]
      b <- parts[3]
      
      # Build the rate Flag column name
      flagCol <- paste("migModel.rateMatrixFlag_", a, "_", b, sep="")
      
      # Remove values where rate flag == 0
      values <- values[logTable[, flagCol] != 0]
    }
    
    # Remove NAs, if present
    values <- values[is.na(values) == FALSE]
    
    # Note and skip those where the value is always the same
    range <- range(values)
    if(range[1] == range[2]){
      #cat(paste("Parameter \"", colNamesToPlot[i],
      #          "\" always has single value: ", range[1], "\n", sep=""))
      next
    }
    
    if(length(essValues) > 1){
      essValues[i] <- calculateEffectiveSampleSize(values)
    }else{
      essValues[i] <- length(essValues[i])
    }
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

plotPosteriorSupportForEachDemeAsRoot <- function(logTable, demes){
  
  # Reset Margins
  currentMar <- par()$mar
  par(mar=c(7, 0, 3, 0))
  
  # Count the number of times each deme was assigned to the root
  rootCounts <- table(logTable$treePrior.rootColor)
  
  # Add in the root counts for the demes that weren't counted
  counts <- c()
  for(index in seq_along(demes)){
    counts[index] <- ifelse(as.character(index-1) %in% names(rootCounts), rootCounts[[as.character(index-1)]], 0)
  }
  
  # Plot the counts
  barplot(counts, yaxt="n", 
          names=demes,
          main="Assignment of Demes to Root State",
          cex.names=0.5, las=2)
  
  # Reset Margins
  par(mar=currentMar)
}

summarisePosteriorLogTable <- function(logFile, logTable, demes, code=2, arrowFactor, nBootstraps=1000,
                                       genomeSize, date){
  
  # Open a pdf file
  plotsFile <- paste0(substr(logFile, 1, nchar(logFile)-4), ".pdf")
  pdf(plotsFile)
    
  # Note which parameters to examine posterior support for
  colsToCalculateESS <- colnames(logTable)[
    grepl(x=colnames(logTable), pattern="Sample|rateMatrixFlag|forward|Count") == FALSE]
    
  # Plot the ESS values of each parameter estimated
  plotParameterESSValues(logTable, colsToCalculateESS)
    
  # Plot the posterior support for each deme as source
  plotPosteriorSupportForEachDemeAsRoot(logTable, demes)
    
  # Plot the population size estimates
  plotPopulationSizes(logTable, demes, alpha=0.5)
    
  # Examine the substitution rate estimates
  examineSubstitutionRateEstimates(logTable, genomeSize)
    
  # Produce a migration rate estimation figure - weight by rate flags
  # Diagrams designed with code = 2 (FORWARDS) in mind
  migrationRateEstimates <- plotMigrationRates(logTable, demeStructure, code, arrowFactor)
    
  # Plot the forward migration rate posterior distributions and their flags
  plotMigrationRatePosteriors(logTable, demeStructure)
    
  # Plot the parameter traces
  plotParameterTraces(logTable, colsToCalculateESS)
    
  # Calculate the acim
  migrationRateEstimates[[analysis]][["AICM"]] <- 
    calculateAICM(logTable$treeLikelihood1, nBootstraps)
    
  # Close the pdf file
  dev.off()
  
  return(migrationRateEstimates)
}


calculateForwardMigrationRates <- function(logTable){
  
  # Get the names of the backward in time migration rate estimates
  migrationRateCols <- colnames(logTable)[
    grepl(colnames(logTable), pattern = "migModel.rateMatrix_")]
  
  # For each backward in time migration rate calculate the forward migration rate
  # FMR_ab = BMR_ba * (Nb / Na)
  #   MR: Migration rate (F - Forward, B - Backward)
  #   N: Effective population size
  #   Demes: a, b
  # Equation taken from second paragraph of Methods section in:
  # De Maio et al. 2015 - New routes to phylogeography ...
  #   BMR_ba = FMR_ab * (Na / Nb)
  #   ->  FMR_ab = BMR_ba / (Na / Nb)
  #     ->  FMR_ab = BMR_ba * (Nb / Na)
  #
  ### NOTE: Backward rates are multiplied by rate flag before being used ###
  # - Converts estimates to 0 when flag = 0 (i.e. rate turned off). If rate isn't likely, then they'll be a lot of zeros that will drag down estimate
  # - Set to NA when flag = 0
  
  for(backwardMigrationRateCol in migrationRateCols){
    
    # Get the demes involved
    parts <- strsplit(backwardMigrationRateCol, split="_")[[1]]
    a <- parts[2]
    b <- parts[3]
    
    #### Multiply the backward rates by the rate flag column ####
    backwardRate <- logTable[, backwardMigrationRateCol] * logTable[, paste("migModel.rateMatrixFlag_", a, "_", b, sep="")]
    
    # Get the estimate population sizes for a and b
    popASizes <- logTable[, paste("migModel.popSize_", a, sep="")]
    popBSizes <- logTable[, paste("migModel.popSize_", b, sep="")]
    
    # Calculate forward rate
    forwardMigrationRateCol <- paste("migModel.forwardRateMatrix_", b, "_", a, sep="")
    logTable[, forwardMigrationRateCol] <- backwardRate * (popASizes/popBSizes)
    
    # Convert the rates to NAs when flag set to zero
    logTable[, forwardMigrationRateCol][logTable[, paste("migModel.rateMatrixFlag_", a, "_", b, sep="")] == 0] <- NA
  }
  
  return(logTable)
}

readLogFile <- function(logFile, burnInProp=0.1){
  
  # Read in the file as table
  logTable <- read.table(logFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  
  # Replace "N" in table with NAs
  logTable[logTable == "N"] <- NA
  
  # Remove the burn-in
  burnIn <- round(burnInProp * nrow(logTable), digits=0)
  logTable <- logTable[burnIn:nrow(logTable), ]
  
  # Calculate the forward rates
  logTable <- calculateForwardMigrationRates(logTable)
  
  return(logTable)
}
