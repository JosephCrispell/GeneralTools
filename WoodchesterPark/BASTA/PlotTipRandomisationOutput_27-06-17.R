####################
# General settings #
####################

# Create a path variable
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/BASTA/TemporalSignal/"

# Note FASTA alignment size
seqLength <- 9464 # Taken from: sequences_Prox-10_plusRef_rmResequenced_SNPCov-0.1_28-10-16.fasta

# Note the burn in proportion
burnIn <- 0.1

#############################################
# Get the actual substitution rate estimate #
#############################################

# Read in the log table
file <- paste(path, "RateEstimation_HKY_StrictClock_ConstantSize_16-06-17/",
              "RateEstimation_HKY_StrictClock_ConstantSize_16-06-17.log",
              sep="")
logTable <- read.table(file, header=TRUE, stringsAsFactors=FALSE)

# Remove the burn in
logTable <- logTable[(burnIn*nrow(logTable)):nrow(logTable), ]

# Store the substitution rate estimate
actualRate <- logTable$clockRate * seqLength

################################################################
# Get the substution rate estimates from date randomisaed data #
################################################################

# Note date and number of replicates
path <- paste(path, "old/", sep="")
date <- "19-06-17"
nReplicates <- 10

# Read in the substitution rate estimates based on the randomised tip dates and store
rateEstimatesFromRandomised <- getRateEstimatesFromDateRandomisationRuns(burnIn, 
                                                                         path, date,
                                                                         nReplicates,
                                                                         seqLength,
                                                                         10001)

#################################################################
# Compare the rates estimated on the actual and randomised data #
#################################################################

# Open a PDF
file <- paste(path, "TipRandomisationResults_27-06-17.pdf", sep="")
pdf(file)

# Create a table with all the rate estimates
rateEstimates <- data.frame("Actual" = actualRate, stringsAsFactors=FALSE)
rateEstimates[, 2:(ncol(rateEstimatesFromRandomised)+1)] <- rateEstimatesFromRandomised

# Plot the rate estimates
boxplot(rateEstimates, las=1, border=c("red", rep("black", nReplicates)), 
        main="Substitution rate estimates based\n upon actual and randomised dates", pch=20, 
        outcol=setAlphaOfColours(c("red", rep("black", nReplicates)), alpha=0.5),
        names=c("Actual", 1:nReplicates),
        ylab="Substitution Rate Estimate (per genome per year)")
quantiles <- quantile(actualRate, probs=c(0.025, 0.975))
points(x=c(0.5, nReplicates+1.5), y=c(quantiles[1], quantiles[1]), type="l", lty=2, 
       col=2, lwd=2)
points(x=c(0.5, nReplicates+1.5), y=c(quantiles[2], quantiles[2]), type="l", lty=2, 
       col=2, lwd=2)

# Add points to show 95% bounds
for(i in 1:ncol(rateEstimates)){
  
  # Get the upper and lowerr bounds
  quantiles <- quantile(rateEstimates[, i], probs=c(0.025, 0.975))
  
  # Add points to figure
  points(x=c(i,i), y=quantiles, pch=16, col="blue")
}




# Quantify the difference between the rates estimated on actual and randomised data
nSamples <- 10000

differences <- sampleRateEstimateDistributionsAndCalculateDifferences(nSamples,
                                                                      rateEstimatesFromRandomised,
                                                                      actualRate)

# Plot the differences
boxplot(differences, las=1, pch=20, outcol=rgb(0,0,0, 0.5), names=1:nReplicates,
        main=paste("Difference between substitution rate estimates based on actual
and randomised dates using ", nSamples, " paired samples", sep=""))
points(x=c(0.5, nReplicates+1.5), y=c(0,0), type="l", col="red", lty=2, lwd=2)

# Add points to show 95% bounds
for(i in 1:ncol(differences)){
  
  # Get the upper and lowerr bounds
  quantiles <- quantile(differences[, i], probs=c(0.025, 0.975))
  
  # Add points to figure
  points(x=c(i,i), y=quantiles, pch=16, col="blue")
}

dev.off()

#############
# Functions #
#############

sampleRateEstimateDistributionsAndCalculateDifferences <- function(nSamples,
                                                                   rateEstimatesFromRandomised,
                                                                   actualRate){
  
  differences <- data.frame(matrix(nrow=nSamples, ncol=nReplicates), stringsAsFactors=FALSE)
  
  for(col in 1:ncol(rateEstimatesFromRandomised)){
    
    # Compare n samples from current distribution based upon random dates and that 
    # on the true dates
    samples1 <- sample(actualRate, nSamples, replace=TRUE)
    samples2 <- sample(rateEstimatesFromRandomised[, col], nSamples, replace=TRUE)
    
    # Note the difference between the samples
    differences[, col] <- samples1 - samples2
  }
  
  return(differences)
}

setAlphaOfColours <- function(colours, alpha){
  
  output <- c()
  for(i in 1:length(colours)){
    output[i] <- setAlpha(colours[i], alpha)
  }
  
  return(output)
}
setAlpha <- function(colour, alpha){
  
  rgbInfo <- col2rgb(colour)
  
  output <- rgb(rgbInfo["red", 1], rgbInfo["green", 1], rgbInfo["blue", 1], alpha=alpha*255,
                maxColorValue=255)
  
  return(output)
}

getRateEstimatesFromDateRandomisationRuns <- function(burnIn, path, date, nReplicates,
                                                      seqLength, nSamples){
 
  # Initialise a table to the substitution rate estimates from each date randomisation
  rateEstimates <- data.frame(matrix(nrow=nSamples - 
                                       round((burnIn*nSamples), digits=0),
                                     ncol=nReplicates), stringsAsFactors=FALSE)
  
  # Get the substitution rate estimates
  for(i in 1:nReplicates){
    
    # Read in the log table
    file <- paste(path, "TipRandomisation_", date, "_", i, "/",
                  "TipRandomisation_", date, "_", i, ".log",
                  sep="")
    logTable <- read.table(file, header=TRUE, stringsAsFactors=FALSE)
    
    # Remove the burn in
    logTable <- logTable[(burnIn*nrow(logTable)):nrow(logTable), ]
    
    rateEstimates[, i] <- logTable$clockRate * seqLength
  }
  
  return(rateEstimates)
}

