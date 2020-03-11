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

#### Read in the BASTA .log files ####

# Note the number of replicates that were done
nReplicates <- 3

# Note the file prefix - without replicate
prefix <- "BASTA_equal_relaxed"

# Note the date the BASTA xmls were created
date <- "02-03-20"

# Read in the log tables into a single log table
logTables <- readLogFiles(path, nReplicates, prefix, date)

#### Create summary plots for each log table ####

migrationRateEstimates <- summarisePosteriorLogTable(logTables, demes, code=2, arrowFactor=20, nBootstraps=1000,
                                                     genomeSize, date)

#### Count transitions between demes estimated in .trees file ####

# Note the path to the JAVA jar tool that counts transitions in .trees file
pathToTransitionCountingJarFile <- file.path(path, "CountTransitions_07-06-19.jar")

# Count the transitions between demes for each posterior distribution of trees
transitionCounts <- countTransitionsForEachReplicate(names(logTables), burnInProp=0.1, demes, 
                                                     pathToTransitionCountingJarFile)

#### Plot the estimated transition counts and rates distributions ####

pdf(paste0(path, "BASTAResults_", date, ".pdf"), width=10, height=17)

par(mfrow=c(2,1))

# Summarise the lineage transition rate estimates
plotSummaryOfLineageTransitionRates(migrationRateEstimates,"4DemeCumbria", nReplicates, date, flagCex=1, varyingOnly=TRUE,
                                    label="a")

# Summarise the transition count distributions
plotTransitionCountDistributionSummaries(transitionCountTables, varyingOnly=TRUE, label="a")

dev.off()

#### Report the AICM score ####

# Examine the AICM scores for the varying and equal models
pdf(paste0(path, "Cumbria/Figures/BASTAResults_AICM_", date, ".pdf"))
plotAICMFromVaryingAndEqualModels(migrationRateEstimates)
dev.off()

### NOT SURE ABOUT NOT INCLUDING REPLICATES!!!??!?!?!?!?!?!

#### FUNCTIONS - plot transition rates ####

summariseLineageTransitionRates <- function(migrationRateEstimates){
  
  # Get the replicate names
  replicates <- names(migrationRateEstimates)
  
  # Get the names of the different transition events
  eventLabels <- names(migrationRateEstimates[[1]])
  eventLabels <- eventLabels[-length(eventLabels)] # Ignoring AICM score stored at end
  
  # Initialise a list to store a summary of each event from each replicate
  eventSummaries <- list("names"=replicates)
  for(label in eventLabels){
    eventSummaries[[label]] <- list("medians"=c(), "uppers"=c(), "lowers"=c(), "flags"=c())
  }

  # Examine each of the replicates
  for(replicateIndex in seq_along(replicates)){
    
    # Examine each of the events for the current replicate
    for(event in eventLabels){
      
      # Get the migration rate values
      values <- migrationRateEstimates[[replicates[replicateIndex]]][[event]]
      
      # Estimate the flag proportion support
      flagSupport <- sum(ifelse(is.na(values), 0, 1)) / length(values)
      
      # Calculate median of the values
      median <- median(values, na.rm=TRUE)
      
      # Calculate the upper and lower bounds
      quantiles <- quantile(values, probs=c(0.025, 0.975), na.rm=TRUE)
      
      # Store all the values calculated 
      eventSummaries[[event]]$medians[replicateIndex] <- median
      eventSummaries[[event]]$uppers[replicateIndex] <- quantiles[2]
      eventSummaries[[event]]$lowers[replicateIndex] <- quantiles[1]
      eventSummaries[[event]]$flags[replicateIndex] <- flagSupport
    }
  }
  
  return(eventSummaries)
}

plotSummaryOfLineageTransitionRates <- function(migrationRateEstimates, date, demeNames, flagCex=1,
                                                priors=FALSE, varyingOnly=FALSE, label=NULL){
  
  # Get and set the current margins
  currentMar <- par()$mar
  par(mar=c(10,4.1,4,0.5))
  
  # Summarise the lineage transition rates
  rateSummaries <- summariseLineageTransitionRates(migrationRateEstimates)
  
  # Note the number of replicates
  nReplicates <- length(migrationRateEstimates)
  
  # Note the locations of the bars
  replicatesPad <- seq(from=-0.4, to=0.4, length.out=nReplicates)
  
  # Note the Y axis ticks
  yRange <- range(c(lowers, uppers))
  at <- c(yRange[1], 0.001, 0.01, 0.1, 1)
  
  # Create an empty plot
  plot(x=NULL, y=NULL, ylim=log(range(at)), xlim=c(0.5, 5.5), bty="n", las=1,
       xaxt="n", ylab="Per lineage transition rate per year", xlab="", yaxt="n")
  
  # Create plot label
  main <- "Estimated using genomic data"
  if(priors){
    main <- "Estimated without genomic data (sampling from priors)"
  }
  title(main, col.main=ifelse(priors, "red", "black"), line=1, cex.main=1.5)
  
  # Add plot label
  if(is.null(label) == FALSE){
    mtext(label, side=3, line=1, at=axisLimits[1], cex=2.5)
  }
  
  # Add the Y labels
  axis(side=2, at=log(at), labels=round(at, digits=3), las=1, xpd=TRUE, line=-0.2)
  
  # Get the axisLimits
  axisLimits <- par()$usr
  yLength <- axisLimits[4] - axisLimits[3]
  
  # Add the X axis labels
  labels <- c("badgers-to-cattle", "cattle-to-badgers")
  axis(side=1, at=c(1, 2), labels=labels, las=2)
  yBracketHeight <- axisLimits[3] - (0.3*yLength)
  yBracketHeights <- c(yBracketHeight+0.2,
                       yBracketHeight, yBracketHeight,
                       yBracketHeight-0.2,
                       yBracketHeight, yBracketHeight,
                       yBracketHeight+0.2)
  points(x=c(1, 1, 1.5, 1.5, 1.5, 2, 2), y=yBracketHeights, xpd=TRUE, type="l")
  points(x=c(3, 3, 3.5, 3.5, 3.5, 4, 4), y=yBracketHeights, xpd=TRUE, type="l")
  axis(side=1, at=c(1.5, 3.5), labels=c("Cumbria", "Northern Ireland (NI)"), tick=FALSE, line=8)
  
  # Loop through each analysis
  for(index in seq_along(names)){
    
    # Note the indices of the to and from IDs
    fromIndex <- 6
    toIndex <- 7
    if(priors){
      fromIndex <- 7
      toIndex <- 8
    }
    
    # Get the location on the X axis for the current analysis
    parts <- strsplit(names[index], split="_")[[1]]
    replicate <- as.numeric(parts[4])
    from <- getDemeNamesForDemeStructure("4DemeCumbria", as.numeric(parts[fromIndex]))
    to <- getDemeNamesForDemeStructure("4DemeCumbria", as.numeric(parts[toIndex]))
    pad <- replicatesPad[replicate] 
    position <- ifelse(parts[2] == "varying", rateLocations[[paste0(from, "-to-", to)]] - pad,
                       rateLocations[[paste0(from, "-to-", to)]] + pad)
    
    # Plot the a summary of the distribution
    points(x=c(position, position), y=log(c(lowers[index], uppers[index])), type="l", lwd=2,
           col=ifelse(parts[2] == "varying", "grey", "black"))
    points(x=position, y=log(medians[index]), pch=19)
    
    # Plot the rate flag value (proportion of steps rate turned on for)
    text(x=position, y=log(uppers[index]), labels=paste0(" ", round(flags[index], digits=2)), srt=90, adj=0, cex=flagCex)
  }
  
  # Add legend
  if(varyingOnly == FALSE){
    legend("topleft", legend=c("varying", "equal"), text.col=c("grey", "black"), bty="n")
  }
  
  # Reset the margins
  par(mar=currentMar)
}

#### FUNCTIONS - transition counts ####

countTransitionsForEachReplicate <- function(logFiles, burnInProp=0.1, demes, pathToJarFile){
  
  # Initialise a variable to store the combined transition counts tables
  transitionCountsTables <- list()
  
  # Examine each replicate
  for(logFile in logFiles){
    
    # Note the trees file name for the current replicate
    treesFile <- paste0(substr(logFile, 1, nchar(logFile) - 4), ".trees")
    
    # Build the output file for the counts
    countsFile <- paste0(substr(logFile, 1, nchar(logFile) - 4), "_TransitionCounts.txt")
    
    # Count the transitions between demes on the current posterior distribution of trees
    countTransitionsOnPosteriorTrees(countsFile, logFile, treesFile, demes, pathToJarFile)
    
    # Read in the transition count tables
    transitionCountTable <- readInBASTATransitionCounts(countsFile, burnInProp)
    
    # Store the transition counts table
    transitionCountsTables[[logFile]] <- transitionCountTable
  }
  
  return(transitionCountsTables)
}

readInBASTATransitionCounts <- function(countsFile, burnInProp=0.1){
          
  # Read in the file as table
  countTable <- read.table(countsFile, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
          
  # Re-order the table by the posterior sample
  countTable <- countTable[order(countTable$Sample), ]
          
  # Remove the burn-in
  burnIn <- round(burnInProp * nrow(countTable), digits=0)
  countTable <- countTable[burnIn:nrow(countTable), ]
  
  return(countTable)
}

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

readLogFiles <- function(path, nReplicates, prefix, date, samplingFromPriors=FALSE, burnInProp=0.1){
  
  # Initialise a store each log table for each replicate
  logTables <- list()
  
  # Examine each replicate
  for(replicate in seq_len(nReplicates)){
    
    # Note the log file name for the current replicate
    logFile <- file.path(path, paste0(prefix, "_", replicate, "_", date), 
                         paste0(prefix, "_", replicate, "_", date, ".log"))
    
    # Check if sampling from priors
    if(samplingFromPriors){
      logFIle <- file.path(path, paste0(prefix, "_PRIOR_", replicate, "_", date), 
                           paste0(prefix, "_PRIOR_", replicate, "_", date, ".log"))
    }
    
    # Read in the log table
    logTables[[logFile]] <- logTable
  }
  
  return(logTables)
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

plotMigrationRatePosteriors <- function(logTable, demeNames, alpha=0.5){
  
  # Get the forward migration rate estimates and their associated flags
  rateEstimates <- list()
  rateFlags <- list()
  for(col in colnames(logTable)){
    
    if(grepl(col, pattern="migModel.forwardRateMatrix_") == TRUE){
      
      # Get the demes involved
      parts <- strsplit(col, split="_")[[1]]
      a <- parts[2]
      b <- parts[3]
      
      # Skip rate if never estimate e.d. cattle-outer -> badger-inner
      if(length(unique(logTable[, col])) == 1){
        next
      }
      
      # Build migration rate name
      rateName <- paste(demeNames[as.numeric(a) + 1], "->", demeNames[as.numeric(b) + 1])
      
      # Store the posterior distribution
      rateEstimates[[rateName]] <- logTable[, col]
      
      # Find the rate flag values for this rate (note that flags will be for backward rate from b to a)
      rateFlags[[rateName]] <- logTable[, paste("migModel.rateMatrixFlag_", b, "_", a, sep="")]
    }
  }
  
  # Define a grid of plots
  nPlots <- ceiling(sqrt(length(rateEstimates)))
  
  # Set the plotting window dimensions
  par(mfrow=c(nPlots, nPlots))
  
  # Set the margins
  origMar <- par("mar")
  par(mar=c(2, 0, 2, 0))
  
  # Plot a trace for each parameter
  for(key in names(rateEstimates)){
    
    values <- rateEstimates[[key]]
    values[is.na(values)] <- 0
    
    hist(values, xlab="", ylab="", yaxt="n", bty="n",
         main=key, cex.main=0.5)
  }
  
  # Reset the plotting window dimensions
  par(mfrow=c(1,1))
  
  # Change the margin sizes
  par(mar=c(20, 4.1, 4.1, 2.1))
  
  # Plot the rate flags
  plot(x=NULL, y=NULL, 
       las=1, bty="n",
       main="Forward migration rate posterior flags", xlab="", ylab="Proportion MCMC steps flag = ON",
       ylim=c(0,1), xlim=c(1, length(rateFlags)), xaxt="n")
  
  # Add in the points
  for(i in 1:length(rateFlags)){
    
    points(x=i, y=mean(rateFlags[[i]]))
  }
  
  # Get the axis limits
  axisLimits <- par("usr")
  
  # Add X axis
  axis(side=1, at=1:length(rateFlags), labels=names(rateFlags), cex.axis=0.75, las=2)
  
  # Reset the margin sizes
  par(origMar)
  
}

createEmptyPlot <- function(){
  par(mar=c(0,0,0,0))
  plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", bty="n",
       ylab="", xlab="")
}

plotDemesAndMigrationRates <- function(arrowWeights, code=2){
  
  # Create empty plot
  createEmptyPlot()
  
  # Note the x and y positions
  x <- c(0.15, 0.85)
  y <- c(0.5, 0.5)
  
  # Add labels
  text(x=x, y=y, labels=c("badgers", "cattle"), col=c("red", "blue"))
  
  # badger -> cow
  if(arrowWeights[["0_1"]] != 0){
    arrows(x0=x[1]+0.1, x1=x[2]-0.1, y0=y[1]-0.1, y1=y[2]-0.1,
           code=code, lwd=arrowWeights[["0_1"]])
  }
  # cow -> badger
  if(arrowWeights[["1_0"]] != 0){
    arrows(x0=x[2]-0.1, x1=x[1]+0.1, y0=y[2]+0.1, y1=y[1]+0.1,
           code=code, lwd=arrowWeights[["1_0"]])
  }
}

timesValuesInListByValue <- function(list, value){
  
  output <- list()
  for(key in names(list)){
    output[[key]] <- list[[key]] * value
  }
  
  return(output)
}

divideValuesInListByMax <- function(arrowRates){
  max <- max(unlist(arrowRates), na.rm=TRUE)
  output <- timesValuesInListByValue(arrowRates, 1/max)
  
  return(output)
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
    
    # Skip rate if never estimate e.d. cattle-outer -> badger-inner
    if(length(unique(logTable[, col])) <= 2){
      next
    }
    
    # Get the directional information (i.e. 0_1, 2_0)
    parts <- strsplit(col, split="_")[[1]]
    direction <- paste(parts[length(parts) - 1], "_", parts[length(parts)], sep="")
    
    # Store a summary statistic for the rate distribution
    arrowRates[[direction]] <- logTable[, col]
  }
  
  return(arrowRates)
}

plotMigrationRates <- function(logTable, code, arrowFactor){
  
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

summarisePosteriorLogTable <- function(logTables, demes, code=2, arrowFactor, nBootstraps=1000,
                                       genomeSize, date){
  
  # Initialise a list to store the estimate migration rates from each log table
  migrationRateEstimates <- list()
  
  # Examine each of the log tables
  for(logFile in names(logTables)){
    
    # Get the log table for the current logFile
    logTable <- logTables[[logFile]]
    
    # Open an PDF for the plots
    plotsFile <- paste0(substr(logFile, 1, nchar(logFile) - 4), ".pdf")
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
    migrationRateEstimates[[logFile]] <- plotMigrationRates(logTable, code, arrowFactor)
    
    # Plot the forward migration rate posterior distributions and their flags
    plotMigrationRatePosteriors(logTable, demes)
    
    # Plot the parameter traces
    plotParameterTraces(logTable, colsToCalculateESS)
    
    # Calculate the acim
    migrationRateEstimates[[logFile]]$AICM <- calculateAICM(logTable$treeLikelihood1, nBootstraps)
    
    # Close the pdf file
    dev.off()
  }

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
