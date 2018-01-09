
### Remember to load the functions found at the bottom of this script before running 
### code below!

#######################
# Load BEAST log file #
#######################

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/"

# Read in the log file
logFile <- paste(path, "test.log", sep="")
logTable <- readLogFile(logFile, burnInProp=0.1)

# Note the columns to examine
colNamesToPlot <- c("posterior", "treeLikelihood1", "ucedMean")

# Open a pdf
pdfFile <- paste(path, substr(logFile, 1, nchar(string)-3), "pdf", sep="")
pdf(pdfFile)

##############################
# Calculate model AICM score #
##############################

aicmScore <- calculateAICM(logLikelihoodValues=logTable$treeLikelihood1,
                           nBootstraps=1000)

##################################
# Calculate parameter ESS values #
##################################

parameterESSValues <- calculateParameterESSValues(logTable, colNamesToPlot, plot=TRUE)

#########################
# Plot parameter traces #
#########################

plotParameterTraces(logTable, colNamesToPlot)

####################################
# Plot each parameter distribution #
####################################

# Proportion used to remove outliers if distorting plot
plotParameterDistributions(logTable, colNamesToPlot, propDistributionToPlot=0.999)

#############################
# Examine the banana effect #
#############################

# Checking to see how correlated the estimates of root height and substitution are
# Should be a random blob, if not much information available then will form a banana
# shape
# Notes:
# - Genome size used to provide substitution rate at the genome level
# - Mutation rate can be stored under: ucedMean or mutationRate
plotBananaEffect(mutationRateColumn="ucedMean", logTable, genomeSize=4349904)

# Close the pdf
dev.off()

#############
# FUNCTIONS #
#############

plotBananaEffect <- function(mutationRateColumn, logTable, genomeSize){
  
  # Get the substitution rate estimates
  rateEstimates <- logTable[, mutationRateColumn] * genomeSize

  # Define a colour palette based upon the likelihood
  rbPal <- colorRampPalette(c('red','blue'))
  colours <- rbPal(10)[as.numeric(cut(logTable$posterior,breaks = 10))]
  colours <- setAlphas(colours, 0.5)
  
  # Plot the substitution rate distribution versus the root height
  plot(x=rateEstimates, y=logTable$tree.height, pch=20, 
       xlab="Substitution Rate (per Genome per Year)", ylab="Root Height (years)",
       main="Substitution Rate versus Root Height",
       col=colours, las=1, bty="n")
  
  legend("topright", legend=c("Likelihood:", "High", "Low"),
         pch=20, col=c("white", "red", "blue"), bty="n")
}

plotParameterDistributions <- function(logTable, colNamesToPlot, propDistributionToPlot){

  # Set the margins
  par(mar=c(5.1, 4.1, 0.5, 0.5))

  # Note the upper and lower bound value
  boundValue <- 0.5 * (1 - propDistributionToPlot)
  
  # Plot a trace for each parameter
  for(col in colNamesToPlot){
    
    # Calculate the mean, median and sd
    mean <- mean(logTable[, col])
    median <- median(logTable[, col])
    sd <- sd(logTable[, col])
    
    # Remove outliers if requested
    quantiles <- quantile(logTable[, col], probs=c(boundValue, 1-boundValue))
    
    # Get the values within the upper and lower bounds
    distribution <- logTable[, col]
    distribution <- distribution[distribution >= quantiles[[1]] & 
                                 distribution <= quantiles[[2]]]

    # Plot the distribution
    hist(distribution, breaks=50, las=1, xlab=col, main="")
    
    # Add median line
    axisLimits <- par("usr") 
    points(x=c(median, median), y=c(0, axisLimits[4]), type="l", lty=2, col="red", lwd=2)
    
    # Add a note about proportion distribution plotted
    legend("topleft", 
           legend=c(paste("Prop plotted = ", propDistributionToPlot),
                    paste("mean = ", round(mean, digits=2)),
                    paste("median = ", round(median, digits=2)),
                    paste("sd = ", round(sd, digits=2))),
           bty="n", cex=0.75)
  }
  
  # Reset the margins
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

readLogFile <- function(file, burnInProp){
  
  # Read in the log table
  read.table(logFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  
  # Replace "N" in table with NAs
  # Note that BASTA in combination with rate flasg can sometimes introduce "N"s
  logTable[logTable == "N"] <- NA
  
  # Remove the burn-in
  burnIn <- round(burnInProp * nrow(logTable), digits=0)
  logTable <- logTable[burnIn:nrow(logTable), ]
  
  return(logTable)
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

calculateParameterESSValues <- function(logTable, colNamesToPlot, plot=TRUE){
  
  # Create a table to store the ESS values
  essValues <- data.frame(Parameter=colNamesToPlot, ESS=rep(NA, length(colNamesToPlot)))
  
  # Examine each parameter
  for(i in 1:length(colNamesToPlot)){
    
    # Remove NAs, if present
    values <- as.numeric(logTable[, colNamesToPlot[i]])
    values <- values[is.na(values) == FALSE]
    
    # Note and skip those where the value is always the same
    range <- range(values)
    if(range[1] == range[2]){
      essValues[i, "ESS"] <- NA
    }
    
    # Calculate the ESS value
    essValues[i, "ESS"] <- calculateEffectiveSampleSize(values)
  }
  
  # Plot ESS values if requested
  if(plot == TRUE){
    
    # Set the margins
    par(mar=c(0,6,2,0.5)) # bottom, left, top, right
    
    # Plot the ESS values
    barplot(essValues[, 2], las=1, names=colNamesToPlot,
            horiz=TRUE, xaxt="n", cex.names=0.8,
            main="Parameter Effective Sample Sizes")
    abline(v=100, lty=2, col="red")
    abline(v=1000, lty=2, col="blue")
    
    # Get the axis limits of the above plot
    # Add median line
    axisLimits <- par("usr")
    text(x=c(300, 1300), y=c(axisLimits[3], axisLimits[3]),
         labels=c("100", "1000"), cex=0.5, col=c("red", "blue"),
         pos=3)
    
    # Reset Margins
    par(mar=c(5.1, 4.1, 4.1, 2.1))
  }
  
  return(essValues)
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