#### Load the data ####

# Set the path
path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/BASTA/Simulations/119-03_2Deme_varying_strict_28-11-18/"

# Note the trees file
treesFile <- paste0(path, "2Deme_equal_relaxed_10-04-18.trees")

# Load the log file
logFile <- paste0(path, "2Deme_varying_strict_28-11-18.log")
log <- read.table(logFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Replace "N" in table with NAs
log[log == "N"] <- NA

# Calculate the forward rates
log <- calculateForwardMigrationRates(log)

#### Count the transitions ####

# Run the Java code
pathToJarFile <- "/home/josephcrispell/Desktop/Research/Java/ExecutableJarFiles/CountTransitions_22-11-18.jar"
countTransitionsOnPosteriorTrees(pathToJarFile, treesFile, logFile)

# Load the transition counts file
countsFile <- paste0(substr(logFile, 1, nchar(logFile)-4), "_TransitionCounts.txt")
counts <- read.table(countsFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)

#### Plot the results ####

# Open a pdf
outputFile <- paste0(substr(logFile, 1, nchar(logFile)-4), "_summary.pdf")
pdf(outputFile, width=14)
par(mfrow=c(1,2))

# Remove the burn-in
burnInProp <- 0.1
log <- log[round(0.1 * nrow(log), digits=0):nrow(log), ]
counts <- counts[round(0.1 * nrow(counts), digits=0):nrow(counts), ]

# Plot the inter-species transition rates
plotInterSpeciesTransitionRateEstimates(log, logFile)

# Plot the transition counts
plotTransitionCounts(counts)

dev.off()

#### FUNCTIONS ####

plotTransitionCounts <- function(counts, withoutPoints=FALSE){
  
  # Note the column names
  columns <- c("Count_BB", "Count_BC", "Count_CB", "Count_CC")
  
  # Plot a boxplot of the inter-species transition rates
  boxplot(counts[, columns], pch=19, outcol=rgb(0,0,0, 0.75), frame=FALSE, xaxt="n", las=1,
          ylab="Number of transitions", ylim=c(0, max(counts[, columns])),
          outline=withoutPoints, main="Number of estimate transitions on posterior trees")
  
  # Add an X axis
  axis(side=1, at=1:length(columns), labels=columns, tick=TRUE)
  
  # Add jittered points on top of boxplot
  if(withoutPoints == FALSE){
    stripchart(counts[, "Count_BB"], at=1, vertical=TRUE, jitter=0.3, 
               method="jitter", add=TRUE, pch=19, col=rgb(0,0,0, 0.002))
    stripchart(counts[, "Count_BC"], at=2, vertical=TRUE, jitter=0.3, 
               method="jitter", add=TRUE, pch=19, col=rgb(0,0,0, 0.002))
    stripchart(counts[, "Count_CB"], at=3, vertical=TRUE, jitter=0.3, 
               method="jitter", add=TRUE, pch=19, col=rgb(0,0,0, 0.002))
    stripchart(counts[, "Count_CC"], at=4, vertical=TRUE, jitter=0.3, 
               method="jitter", add=TRUE, pch=19, col=rgb(0,0,0, 0.002))
  }

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
    
    "3Deme-outerIsBadger"=c("badger-inner", "cow", "unsampled"),
    
    "3Deme-outerIsCattle"=c("badger", "cow-inner", "cow-outer"),
    
    "4Deme"=c("badger-inner", "cow-inner", "cow-outer", "unsampled"),
    
    "6Deme-EastWest"=c("badger-inner", "cow-inner", "cow-outer-east", "cow-outer-west", "unsampled-1", "unsampled-2"),
    
    "6Deme-NorthSouth"=c("badger-inner", "cow-inner", "cow-outer-north", "cow-outer-south", "unsampled-1", "unsampled-2"),
    
    "8Deme-EastWest"=c("badger-inner-east", "badger-inner-west", "cow-inner-east", "cow-inner-west",
                       "cow-outer-east", "cow-outer-west", "unsampled-1", "unsampled-2"),
    
    "8Deme-NorthSouth"=c("badger-inner-north", "badger-inner-south", "cow-inner-north", "cow-inner-south",
                         "cow-outer-north", "cow-outer-south", "unsampled-1", "unsampled-2")
  )
  
  output = demeNames[[demeStructure]]
  if(is.null(number) == FALSE){
    output = demeNames[[demeStructure]][number + 1]
  }
  
  return(output)
}

getMigrationRates <- function(logTable){
  
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
    if(length(unique(logTable[, col])) == 1){
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

plotInterSpeciesTransitionRateEstimates <- function(log, logFile, withoutPoints=FALSE){
  
  # Initialise two arrays to store samples of the posterior sums from each model - sample size relative to AICM weight
  badgerToCowRatePosteriorSumSamples<- c()
  cowToBadgerRatePosteriorSumSamples <- c()
  
  # Initialise arrays to store the rates from badgers to cattle and vice versa
  sumRatesBadgerToCow <- rep(0, nrow(log))
  sumRatesCowToBadger <- rep(0, nrow(log))
    
  # Get the demeStructure
  parts <- strsplit(logFile, "/")[[1]]
  demeStructure <- strsplit(parts[length(parts)], "_")[[1]][1]
  
  # Get the rate estimates from the log table
  migrationRateEstimates <- getMigrationRates(log)
  
  # Examine each rate
  for(key in names(migrationRateEstimates)){
    
    # Get the values for the current migration rate
    values <- migrationRateEstimates[[key]]
    values[is.na(values)] <- 0
    
    # Split the key into its deme numbers
    demeNumbers <- as.numeric(strsplit(key, split="_")[[1]])
    
    # Check if current rate is between badger and cattle populations
    if(grepl(getDemeNamesForDemeStructure(demeStructure, demeNumbers[1]),
             pattern="badger") == TRUE &&
       grepl(getDemeNamesForDemeStructure(demeStructure, demeNumbers[2]),
             pattern="cow") == TRUE){
      
      # Add the current rate values to the growing sum
      sumRatesBadgerToCow <- sumRatesBadgerToCow + values
      
    }else if(grepl(getDemeNamesForDemeStructure(demeStructure, demeNumbers[1]),
                   pattern="cow") == TRUE &&
             grepl(getDemeNamesForDemeStructure(demeStructure, demeNumbers[2]),
                   pattern="badger") == TRUE){
      
      # Add the current rate values to the growing sum
      sumRatesCowToBadger <- sumRatesCowToBadger + values
    }
  }
  
  # Remove any values that are exactly zero
  # - These will result when flag=0 across badger-to-cattle/cattle-to-badger rates estimated
  sumRatesBadgerToCow <- sumRatesBadgerToCow[sumRatesBadgerToCow != 0]
  sumRatesCowToBadger <- sumRatesCowToBadger[sumRatesCowToBadger != 0]
    
  # Plot a boxplot of the inter-species transition rates
  boxplot(sumRatesBadgerToCow, sumRatesCowToBadger, pch=19, outcol=rgb(0,0,0, 0.75), frame=FALSE, xaxt="n", las=1,
          ylab="Per lineage transition rate per year", ylim=c(0, max(c(sumRatesBadgerToCow, sumRatesCowToBadger))),
          outline=withoutPoints, main="Estimated inter-species transition rates")
  
  # Add an X axis
  axis(side=1, at=c(1,2), labels=c("Badger->Cow", "Cow->Badger"), tick=TRUE)
  
  # Add jittered points on top of boxplot
  if(withoutPoints == FALSE){
    stripchart(sumRatesBadgerToCow, at=1, vertical=TRUE, jitter=0.3, 
               method="jitter", add=TRUE, pch=19, col=rgb(0,0,0, 0.02))
    stripchart(sumRatesCowToBadger, at=2, vertical=TRUE, jitter=0.3, 
               method="jitter", add=TRUE, pch=19, col=rgb(0,0,0, 0.02))
  }
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

countTransitionsOnPosteriorTrees <- function(pathToJarFile, treesFile, logFile){
  
  # Build the command to be used
  command <- paste0("java -jar ", pathToJarFile, " ", treesFile, " ", logFile)

  # Run the command
  system(command, wait=TRUE)
}