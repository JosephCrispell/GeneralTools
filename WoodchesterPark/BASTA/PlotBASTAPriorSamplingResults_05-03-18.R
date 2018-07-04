#--------------------------#
#### Read in Log Tables ####
#--------------------------#

# Set the path
path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/BASTA/PriorSampling_10-04-18/"

# Note deme structure to use
demeStructures <- list(
  "2Deme"=2,
  "3Deme-outerIsBoth"=3,
  "3Deme-outerIsBadger"=3,
  "3Deme-outerIsCattle"=3,
  "4Deme"=4,
  "6Deme-EastWest"=6,
  "6Deme-NorthSouth"=6,
  "8Deme-EastWest"=8,
  "8Deme-NorthSouth"=8
)

# Note the date when all analyses were created
date <- "10-04-18"

# Note the population size estimation options
popEstimationTypes <- c("varying", "equal")

# Note the clock model options
clockEstimateTypes <- c("relaxed") # strict not used

# Read in the log tables
logTables <- readInBASTALogTables(date, demeStructures, popEstimationTypes, clockEstimateTypes, path, ignoreIfFlagged=FALSE)

#------------------------#
#### Examine each run ####
#------------------------#

summarisePosteriorLogTables(path, logTables, arrowFactor=20, date=date)


#-----------------#
#### FUNCTIONS ####
#-----------------#

readInBASTALogTables <- function(date, demeStructures, popEstimationTypes,
                                 clockEstimateTypes, path, burnInProp=0.1, ignoreIfFlagged=FALSE){
  
  # Store each of the log tables in a list
  logTables <- list()
  
  for(demeStructure in names(demeStructures)){
    
    for(popEstimationType in popEstimationTypes){
      
      for(clockEstimationType in clockEstimateTypes){
        
        # Build run defining prefix
        prefix <- paste(demeStructure, "_", popEstimationType, "_",
                        clockEstimationType, "_PRIOR_", date, sep="")
        
        # Print progress information
        cat(paste("\rReading: ", prefix, ".log\t\t\t\t\t\t", sep=""))
          
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
        logTable <- calculateForwardMigrationRates(logTable, ignoreIfFlagged)
          
        # Store the tables as a single log table
        logTables[[prefix]] <- logTable
      }
    }
  }
  cat("\rFinished reading in log tables...\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n")
  
  return(logTables)
}

calculateForwardMigrationRates <- function(logTable, ignore=FALSE){
  
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
    
    # Convert the zeros to NAs, if requested. If ignore == TRUE, we want to ignore values when flag = 0
    if(ignore == TRUE){
      logTable[, forwardMigrationRateCol][logTable[, forwardMigrationRateCol] == 0] <- NA
    }
  }
  
  return(logTable)
}

summarisePosteriorLogTables <- function(path, logTables, code=2, arrowFactor, nBootstraps,
                                        genomeSize, date){
  
  # Get analysis names
  analyses <- names(logTables)
  
  # Initialise a list to store the migration rate estimates and AICM
  migrationRateEstimates <- list()
  
  # Create directory for summary plots
  dir.create(paste(path, "SummaryPlots_", date, "/", sep=""), showWarnings = FALSE)
  
  # Examine each analysis
  for(analysis in analyses){
    
    # Open a pdf file
    file <- paste(path, "SummaryPlots_", date, "/", analysis, "_ResultsSummary.pdf",
                  sep="")
    pdf(file)
    
    # Get the deme structure
    demeStructure <- strsplit(analysis, "_")[[1]][1]
    
    # Note progress
    cat(paste("\rExamining log table for: ", analysis, "\t\t\t\t\t\t\t", sep=""))
    
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
    examineSubstitutionRateEstimates(logTable, genomeSize=1)
    
    # Produce a migration rate estimation figure - weight by rate flags
    # Diagrams designed with code = 2 (FORWARDS) in mind
    plotMigrationRates(logTable, demeStructure, code, arrowFactor)
    
    # Plot the forward migration rate posterior distributions and their flags
    plotMigrationRatePosteriors(logTable, demeStructure)
    
    # Plot the parameter traces
    plotParameterTraces(logTable, colsToCalculateESS)

    # Close the pdf file
    dev.off()
  }
  cat("\rFinished examining log tables...\t\t\t\t\t\tn")
  
  # Reset the margins
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

createEmptyPlot <- function(){
  par(mar=c(0,0,0,0))
  plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", bty="n",
       ylab="", xlab="")
}

timesValuesInListByValue <- function(list, value){
  
  output <- list()
  for(key in names(list)){
    output[[key]] <- list[[key]] * value
  }
  
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

divideValuesInListByMax <- function(arrowRates){
  max <- max(getValues(arrowRates), na.rm=TRUE)
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
    
    # Get the directional information (i.e. 0_1, 2_0)
    parts <- strsplit(col, split="_")[[1]]
    direction <- paste(parts[length(parts) - 1], "_", parts[length(parts)], sep="")
    
    # Store a summary statistic for the rate distribution
    arrowRates[[direction]] <- median(logTable[, col], na.rm=TRUE)
  }
  
  return(arrowRates)
}

plotMigrationRatePosteriors <- function(logTable, demeStructure, alpha=0.5){
  
  # Get the deme names
  demeNames <- getDemeNamesForDemeStructure(demeStructure)
  
  # Get the forward migration rate estimates and their associated flags
  rateEstimates <- list()
  rateFlags <- list()
  for(col in colnames(logTable)){
    
    if(grepl(col, pattern="migModel.forwardRateMatrix_") == TRUE){
      
      # Get the demes involved
      parts <- strsplit(col, split="_")[[1]]
      a <- parts[2]
      b <- parts[3]
      
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
  par(mar=c(12, 4.1, 4.1, 2.1))
  
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
  axis(side=1, at=1:length(rateFlags), labels=names(rateFlags), cex=0.5)
  
  # Reset the margin sizes
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  
}

plotMigrationRates <- function(logTable, demeStructure, code, arrowFactor){
  
  # Get the migration rates
  output <- getArrowRates(logTable)
  
  # Normalise those rates to vary between 0 and MAX (arrow factor)
  migrationRates <- divideValuesInListByMax(output)
  migrationRates <- timesValuesInListByValue(migrationRates, arrowFactor)
  
  # Check for "NaN" or NA values
  for(key in names(migrationRates)){
    if(is.nan(migrationRates[[key]]) == TRUE || is.na(migrationRates[[key]]) == TRUE){
      migrationRates[[key]] <- 0
    }
  }
  
  # Plot the demes and associated rates
  if(demeStructure == "2Deme"){
    plot2Deme(migrationRates, code)
  }else if(grepl(pattern="3Deme", x=demeStructure) == TRUE){
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
  
  # badger/badger-inner -> cow/cow-inner
  if(arrowWeights[["0_1"]] != 0){
    arrows(x0=x[1]+0.15, x1=x[2]-0.15, y0=y[1]-0.05, y1=y[2]-0.05,
           code=code, lwd=arrowWeights[["0_1"]])
  }
  # cow/cow-inner -> badger/badger-inner
  if(arrowWeights[["1_0"]] != 0){
    arrows(x0=x[2]-0.15, x1=x[1]+0.15, y0=y[2]+0.05, y1=y[1]+0.05,
           code=code, lwd=arrowWeights[["1_0"]])
  }
  
  if(demeStructure != "3Deme-outerIsCattle"){
    # badger/badger-inner -> outer/badger-outer
    if(arrowWeights[["0_2"]] != 0){
      arrows(x0=x[1]-0.05, x1=x[3]-0.15, y0=y[1]+0.1, y1=y[3]-0.1,
             code=code, lwd=arrowWeights[["0_2"]])
    }
    # outer/badger-outer -> badger/badger-inner
    if(arrowWeights[["2_0"]] != 0){
      arrows(x0=x[3]-0.05, x1=x[1]+0.05, y0=y[3]-0.1, y1=y[1]+0.1,
             code=code, lwd=arrowWeights[["2_0"]])
    }
  }
  
  # cow/cow-inner -> outer/cow-outer/badger-outer
  if(arrowWeights[["1_2"]] != 0){
    arrows(x0=x[2]+0.05, x1=x[3]+0.15, y0=y[2]+0.1, y1=y[3]-0.1,
           code=code, lwd=arrowWeights[["1_2"]])
  }
  # outer/cow-outer/badger-outer -> cow/cow-inner
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

plotPosteriorSupportForEachDemeAsRoot <- function(logTable, demeStructure){
  
  # Reset Margins
  par(mar=c(7, 0, 3, 0))
  
  barplot(table(logTable$treePrior.rootColor), yaxt="n", 
          names=getDemeNamesForDemeStructure(demeStructure),
          main="Assignment of Demes to Root State",
          cex.names=0.5, las=2)
  
  # Reset Margins
  par(mar=c(5.1, 4.1, 4.1, 2.1))}

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
