#### Preparation ####

# Load the required libraries
library(ape) # Reading in FASTA
library(lubridate) # Converting to decimal dates

# Set the path
path <- file.path("~", "Desktop", "BuildingBASTAXML")

# Set the date
date <- format(Sys.Date(), "%d-%m-%y")

# Note the states/demes considered - ALPHABETICAL ORDER
states <- c("badger", "cow")

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
countTransitionsOnPosteriorTrees(countsFile, logFile, treesFile, states, pathToTransitionCountingJarFile)

#### FUNCTIONS - transition counts ####

countTransitionsOnPosteriorTrees <- function(countsFile, logFile, treesFile, states, pathToJarFile){
  
  # Check that transition counts file doesn't exist
  if(file.exists(countsFile) == FALSE){
    
    # Collapse states into string
    states <- paste(states, collapse=",")
    
    # Build the command to run the CountTransitions_DATE.jar file
    command <- paste0("java -jar ", pathToJarFile, " ", treesFile, " ", logFile, " ", states)
    
    # Count the number of transitions between states for the current posterior tree distribution
    system(command, wait=TRUE)
  }
}

#### FUNCTIONS - log file ####

summarisePosteriorLogTable <- function(path, logTable, code=2, arrowFactor, nBootstraps=1000,
                                       genomeSize, date){
    
  # Open a pdf file
  pdf(paste0(directory, analysis, "_ResultsSummary.pdf"))
    
  # Get the deme structure
  demeStructure <- strsplit(analysis, "_")[[1]][1]
    
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
