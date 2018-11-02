#### Read in the simulation output ####

# Set the path variable
path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/Daniel_SimulationOutput/"

# Read in a simulation output table
file <- paste0(path, "sim_05-01.csv")
simOutput <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Read in the herd and sett coordinates
something needs to be here!!!!
  
  
#### Examine the simulation output ####

# Convert the date columns to dates
simOutput$DateSampleTaken <- as.Date(simOutput$DateSampleTaken, format="%Y-%m-%d")
simOutput$InfectionDate <- as.Date(simOutput$InfectionDate, format="%Y-%m-%d")
simOutput$DetectionDate <- as.Date(simOutput$DetectionDate, format="%Y-%m-%d")
simOutput$LastSnpGeneration <- as.Date(simOutput$LastSnpGeneration, format="%Y-%m-%d")

# Briefly summarise the data
summariseSimulationOutput(simOutput)

# Plot the temporal sampling range
plotTemporalRange(simOutput)

# Check that the dates in the simulation output make sense
checkSimulationOutputDates(simOutput)

# Report how many seeds were used in simulation
countSeeds(simOutput)

# Plot the infection dynamics over time
simulationDynamics <- plotInfectionDynamics(simOutput)

# Examine spatial sampling
something needs to be here!!!!

#### Build the sequences ####

# Add a nucleotide sequence, based on the SNPs, for each individual
simOutput <- generateSequences(simOutput)

#### Build BASTA xml file ####

# Select the sequences

#############
# FUNCTIONS #
#############

plotInfectionDynamics <- function(simOutput){
  
  # Get an ordered list of all dates in simulation output
  dates <- sort(unique(c(simOutput$InfectionDate, simOutput$DateSampleTaken)))
  
  # Index the dates
  indexedDates <- indexArrayOfDates(dates)
  
  # Initialise a table to store the population dynamics
  simulationDynamics <- data.frame("Date" = dates, 
                                   "NumberBadgersInfected"=0, "NumberBadgersSampled"=0,
                                   "NumberCattleInfected"=0, "NumberCattleSampled"=0)
  
  # Examine the simulation output and count how many infected and sampled animals there are at any given date
  for(row in seq_len(nrow(simOutput))){
    
    # Progress information
    cat(paste0("\rReading row: ", row, " of ", nrow(simOutput)))
    
    # Skip the root
    if(simOutput[row, "animal_ID"] == "ROOT"){
      next
    }
    
    # Get index of infection date
    simulationDynamics <- countDaysInfectedAndSampled(simulationDynamics,
                                                      infectionDate=simOutput[row, "InfectionDate"],
                                                      samplingDate=simOutput[row, "DateSampleTaken"],
                                                      indexedDates, isCow=simOutput[row, "isCow"])
  }
  cat("\rFinished!\t\t\t\t\t\t\t\t\t\t\t\n")
  
  # Plot the dynamics
  plot(x=simulationDynamics$Date, y=simulationDynamics$NumberBadgersInfected, 
       ylim=c(0, max(simulationDynamics[, 2:5])), 
       type="l", col="red", main="Simulation dynamics", ylab="Number", xlab="Date",
       bty="n", las=1)
  points(x=simulationDynamics$Date, y=simulationDynamics$NumberBadgersSampled,
         type="l", lty=2, col="red")
  points(x=simulationDynamics$Date, y=simulationDynamics$NumberCattleInfected,
         type="l", col="blue")
  points(x=simulationDynamics$Date, y=simulationDynamics$NumberCattleSampled,
         type="l", lty=2, col="blue")
  
  # Add a legend
  legend("topleft", legend=c("N. badgers infected", "N. badgers sampled", 
                             "N. cattle infected", "N. cattle sampled"),
         lty=c(1,2,1,2), col=c("red", "red", "blue", "blue"), bty="n",
         text.col=c("red", "red", "blue", "blue"))
  
  return(simulationDynamics)
}

countDaysInfectedAndSampled <- function(simulationDynamics, infectionDate, samplingDate, indexedDates, isCow){
  
  # Get the index of the infection date
  start <- 1 # Seeds don't have infection date but are infected from day 1
  if(is.na(infectionDate) == FALSE){
    start <- indexedDates[[as.character(infectionDate)]]
  }
  
  # Get the index of the sampling date
  end <- nrow(simulationDynamics) # If no sampling date - infected until last date
  if(is.na(samplingDate) == FALSE){
    end <- indexedDates[[as.character(samplingDate)]]
  }
  
  # Increment the days when animal was infected
  if(isCow == TRUE){
    simulationDynamics[start:end, "NumberCattleInfected"] <- simulationDynamics[start:end, "NumberCattleInfected"] + 1
  }else{
    simulationDynamics[start:end, "NumberBadgersInfected"] <- simulationDynamics[start:end, "NumberBadgersInfected"] + 1
  }
  
  # Increment the days when animal was sampled
  if(end != nrow(simulationDynamics) && isCow == TRUE){
    simulationDynamics[end:nrow(simulationDynamics), "NumberCattleSampled"] <- 
      simulationDynamics[end:nrow(simulationDynamics), "NumberCattleSampled"] + 1
  }else if(end != nrow(simulationDynamics) && isCow == FALSE){
    simulationDynamics[end:nrow(simulationDynamics), "NumberBadgersSampled"] <- 
      simulationDynamics[end:nrow(simulationDynamics), "NumberBadgersSampled"] + 1
  }

  return(simulationDynamics)
}

indexArrayOfDates <- function(array){
  output <- list()
  for(i in 1:length(array)){
    output[[as.character(array[i])]] <- i
  }
  
  return(output)
}

countSeeds <- function(simOutput){
  
  nBadgerSeeds <- length(which(grepl(simOutput$animal_ID, pattern="Badger_seed") == TRUE))
  nCattleSeeds <- length(which(grepl(simOutput$animal_ID, pattern="Cow_seed") == TRUE))
  
  cat(paste("Number of badger seeds =", nBadgerSeeds, 
            "\nNumber of cattle seeds =", nCattleSeeds, "\n"))
}

checkSimulationOutputDates <- function(simOutput){
  
  # Compare sampling dates and last SNP dates - should be EQUAL
  plot(x=simOutput$DateSampleTaken, y=simOutput$LastSnpGeneration, 
       xlab="Date sample taken", ylab="Date last SNP occurred", las=1, bty="n",
       main="Comparing sampling and last SNP dates")
  points(x=simOutput$DateSampleTaken, y=simOutput$DateSampleTaken, lty=2, col="red", type="l")
  
  # Compare sampling and detection dates - should be EQUAL
  plot(x=simOutput$DateSampleTaken, y=simOutput$DetectionDate, 
       xlab="Date sample taken", ylab="Date infection detected", las=1, bty="n",
       main="Comparing detection and sampling dates")
  points(x=simOutput$DateSampleTaken, y=simOutput$DateSampleTaken, lty=2, col="red", type="l")
  
  # Compare sampling and detection dates - infection should be before sampling
  plot(x=simOutput$DateSampleTaken, y=simOutput$InfectionDate, 
       xlab="Date sample taken", ylab="Date infected", las=1, bty="n",
       main="Comparing infection and sampling dates")
  points(x=simOutput$DateSampleTaken, y=simOutput$DateSampleTaken, lty=2, col="red", type="l")
}

summariseSimulationOutput <- function(simOutput){
  
  # Subset the data
  badgers <- simOutput[simOutput$isCow == "FALSE", ]
  sampledBadgers <- badgers[is.na(badgers$DateSampleTaken) == FALSE, ]
  cattle <- simOutput[simOutput$isCow == "TRUE", ]
  sampledCattle <- cattle[is.na(cattle$DateSampleTaken) == FALSE, ]
  
  # Report number of sampled cattle and badgers
  cat(paste0("Number sampled badgers: ", nrow(sampledBadgers), " of ", nrow(badgers),
             "\nNumber sampled cattle: ", nrow(sampledCattle), " of ", nrow(cattle), "\n"))
  
  # Report the sampling date range
  cat(paste0("Sampled badgers available from: ", min(sampledBadgers$DateSampleTaken), " to ", max(sampledBadgers$DateSampleTaken),
             "\nSampled cattle available from: ", min(sampledCattle$DateSampleTaken), " to ", max(sampledCattle$DateSampleTaken), "\n"))

  # Count the number of animals in each infection status category
  statusCounts <- countAnimalStatuses(badgers, cattle, sampledBadgers, sampledCattle)
}

countAnimalStatuses <- function(badgers, cattle, sampledBadgers, sampledCattle, plot=TRUE){
  
  statusCounts <- data.frame("Exposed"=c(NA, NA, NA, NA),
                             "Infectious"=c(NA, NA, NA, NA),
                             "TestSensitive"=c(NA, NA, NA, NA),
                             row.names=c("Badgers", "Cattle", "SampledBadgers", "SampledCattle"))
  statusCounts["Badgers", "Infectious"] <- nrow(badgers)
  statusCounts["SampledBadgers", "Infectious"] <- nrow(sampledBadgers)
  cattleCounts <- table(cattle$InfectionStatus)
  statusCounts["Cattle", "Exposed"] <- cattleCounts["EXPOSED"]
  statusCounts["Cattle", "Infectious"] <- cattleCounts["INFECTIOUS"]
  statusCounts["Cattle", "TestSensitive"] <- cattleCounts["TESTSENSITIVE"]
  sampledCattleCounts <- table(sampledCattle$InfectionStatus)
  statusCounts["SampledCattle", "Exposed"] <- sampledCattleCounts["EXPOSED"]
  statusCounts["SampledCattle", "Infectious"] <- sampledCattleCounts["INFECTIOUS"]
  statusCounts["SampledCattle", "TestSensitive"] <- sampledCattleCounts["TESTSENSITIVE"]
  
  if(plot){
    # Plot the counts
    currentMar <- par("mar")
    par(mar=c(4.1, 10, 4.1, 2.1))
    barplot(c(statusCounts["Badger", "Infectious"], statusCounts["SampledBadger", "Infectious"],
              statusCounts["Cattle", "Exposed"], statusCounts["SampledCattle", "Exposed"],
              statusCounts["Cattle", "Infectious"], statusCounts["SampledCattle", "Infectious"],
              statusCounts["Cattle", "TestSensitive"], statusCounts["SampledCattle", "TestSensitive"]),
            names.arg=c("Badger - I", "Sampled Badger - I",
                        "Cattle - E", "Sampled Cattle - E",
                        "Cattle - I", "Sampled Cattle - I",
                        "Cattle - T", "Sampled Cattle - T"), 
            horiz=TRUE, las=1, xlab="Number recorded", main="Numbers of samples from infection categories",
            col=c("red", "red", "blue", "blue", "blue", "blue", "blue", "blue"))
    par(mar=currentMar)
  }
  
  return(statusCounts)
}

plotTemporalRange <- function(simOutput){
  
  # Create two histograms of the sampling dates for badgers and cattle
  histBadgers <- hist(simOutput[simOutput$isCow == "FALSE", "DateSampleTaken"], breaks=10, freq=TRUE, plot=FALSE)
  histCattle <- hist(simOutput[simOutput$isCow == "TRUE", "DateSampleTaken"], breaks=10, freq=TRUE, plot=FALSE)
  
  # Get the y axis limits
  yLim <- c(0, max(c(histBadgers$counts, histCattle$counts)))
  
  # Plot the histograms
  # Note I removed the Y axis as it doesn't plot dates properly
  plot(histBadgers, 
       col=rgb(1,0,0, 0.5), las=1, xlab="Year",
       ylim=yLim, main="Temporal range of sampling", xaxt="n")
  plot(histCattle, col=rgb(0,0,1, 0.5), add=TRUE)
  
  # Add the X axis back in
  Axis(simOutput[simOutput$isCow == "FALSE", "DateSampleTaken"], col="black", side=1)
  
  # Add a legend
  legend("topleft", legend=c("Badgers", "Cattle"), text.col=c("red", "blue"), bty="n")
}

generateSequences <- function(simOutput){
  
  # Get the maximum SNP ID - equal to the number SNPs that occured in the population
  maxSNPID <- getMaxSNPID(simOutput)
  cat(paste0(maxSNPID, " SNPs found in simulated data\n"))
  
  # Create a reference and alternate allele for each SNP
  snpInfo <- createSNPs(maxSNPID)
  
  # Create each individual's sequence
  simOutput$Sequence <- NA
  
  # Examine each individual and create its sequence
  for(row in seq_len(nrow(simOutput))){
    
    # Skip individuals with no SNPs
    if(is.na(simOutput[row, "SNPs_detected"])){
      next
    }
    
    # First assign the reference sequence to the current individual
    sequence <- snpInfo$Reference
    
    # Get the SNPs associated with the current individual
    snps <- as.numeric(strsplit(simOutput[row, "SNPs_detected"], split=";")[[1]])
    
    # Assign each of the current individual's SNPs to its sequence
    for(snp in snps){
      sequence[snp] <- snpInfo$Alternate[snp]
    }
    
    # Convert the nucleotide array to string and store it
    simOutput[row, "Sequence"] <- paste(sequence, collapse="")
  }
  
  return(simOutput)
}

getMaxSNPID <- function(simOutput){
  
  # Initialise a variable to record the maximum SNP ID - the number of SNPs in the population
  maxSNPID <- -Inf
  
  # Examine the simulation data for each animal and store its SNP IDs
  for(row in seq_len(nrow(simOutput))){
    
    # Get the maximum SNP ID from the current individual
    individualMaxSNPID <- max(as.numeric(strsplit(simOutput[row, "SNPs_detected"], split=";")[[1]]))
    
    if(is.na(individualMaxSNPID)){
      cat(paste("No SNPs found in row:", row, "with animal_ID:", simOutput[row, "animal_ID"], "\n"))
      next
    }
    
    # Check if found new maximum
    if(individualMaxSNPID > maxSNPID){
      maxSNPID <- individualMaxSNPID
    }
  }
  
  return(maxSNPID)
}

createSNPs <- function(maxSNPID){
  
  # Initialise a matrix to store the reference alleles and alternates associated with each snp
  snpInfo <- list("Reference"=c(), "Alternate"=c())
  
  # Create the reference allele for each event
  nucleotides <- c("A", "C", "G", "T")
  snpInfo[["Reference"]] <- sample(nucleotides, size=maxSNPID, replace=TRUE)
  
  # Create the alternate allele for each event
  for(i in seq_len(maxSNPID)){
    snpInfo[["Alternate"]][i] <- sample(nucleotides[nucleotides != snpInfo[["Reference"]][i]], size=1)
  }
  
  return(snpInfo)
}
