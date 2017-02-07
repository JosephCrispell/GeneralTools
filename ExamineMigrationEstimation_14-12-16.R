# Read in the log file
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/BASTA/"

folderName <- "BASTA_Test_01-02-17/"

file <- paste(path, folderName,
              "BASTA_Test_01-02-17.log", sep="")

logTable <- read.table(file, header=TRUE, sep="\t")

# Remove the burn-in
burnIn <- round(0.1 * nrow(logTable), digits=0)
logTable <- logTable[burnIn:nrow(logTable), ]

# Open a PDF
file <- paste(path, folderName, "ResultsSummary_01-02-17.pdf", sep="")
pdf(file)

###################
# Migration Rates #
###################

par(mfrow=c(1,1))

# Define line weights
lineWeights_Counts <- c(
  mean(logTable[, "treePrior.Count0to1"]),    # 0 -> 1    innerBadger -> innerCow
  mean(logTable[, "treePrior.Count0to3"]),    # 0 -> 3    innerBadger -> outerBadger
  mean(logTable[, "treePrior.Count1to0"]),    # 1 -> 0    innerCow -> innerBadger
  mean(logTable[, "treePrior.Count1to2"]),    # 1 -> 2    innerCow -> outerCow
  mean(logTable[, "treePrior.Count2to1"]),    # 2 -> 1    outerCow -> innerCow
  mean(logTable[, "treePrior.Count2to3"]),    # 2 -> 3    outerCow -> outerBadger
  mean(logTable[, "treePrior.Count3to0"]),    # 3 -> 0    outerBadger -> innerBadger
  mean(logTable[, "treePrior.Count3to2"])     # 3 -> 2    outerBadger -> outerCow
)

lineWeights_Rates <- c(
  mean(logTable[, "migModel.rateMatrix_0_1"]),    # 0 -> 1    innerBadger -> innerCow
  mean(logTable[, "migModel.rateMatrix_0_3"]),    # 0 -> 3    innerBadger -> outerBadger
  mean(logTable[, "migModel.rateMatrix_1_0"]),    # 1 -> 0    innerCow -> innerBadger
  mean(logTable[, "migModel.rateMatrix_1_2"]),    # 1 -> 2    innerCow -> outerCow
  mean(logTable[, "migModel.rateMatrix_2_1"]),    # 2 -> 1    outerCow -> innerCow
  mean(logTable[, "migModel.rateMatrix_2_3"]),    # 2 -> 3    outerCow -> outerBadger
  mean(logTable[, "migModel.rateMatrix_3_0"]),    # 3 -> 0    outerBadger -> innerBadger
  mean(logTable[, "migModel.rateMatrix_3_2"])     # 3 -> 2    outerBadger -> outerCow
)

# Set arrow direction
code <- 1 # Backward in time (1 for forward)
factor <- 1000

direction <- "BACKWARDS"
if(code == 1){
  direction <- "FORWARDS"
}

title <- paste("Migration Counts ", direction, " in time", sep="")
plot(x=NULL,
     y=NULL,
     xlim=c(1, 10),
     ylim=c(1,10),
     xaxt="n", yaxt="n", xlab="", ylab="", main=title)

arrows(x0=2, y0=4, x1=2, y1=7, code=code, lwd=lineWeights_Counts[1]) # innerBadger -> innerCow
arrows(x0=3, y0=7, x1=3, y1=4, code=code, lwd=lineWeights_Counts[3]) # innerCow -> innerBadger

arrows(x0=8, y0=4, x1=8, y1=7, code=code, lwd=lineWeights_Counts[8]) # outerBadger -> outerCow
arrows(x0=9, y0=7, x1=9, y1=4, code=code, lwd=lineWeights_Counts[6]) # outerCow -> outerBadger

arrows(x0=4, y0=3, x1=7.5, y1=3, code=code, lwd=lineWeights_Counts[2]) # innerBadger -> outerBadger
arrows(x0=7.5, y0=4, x1=4, y1=4, code=code, lwd=lineWeights_Counts[7]) # outerBadger -> innerBadger

arrows(x0=4, y0=7, x1=7.5, y1=7, code=code, lwd=lineWeights_Counts[4]) # innerCow -> outerCow
arrows(x0=7.5, y0=8, x1=4, y1=8, code=code, lwd=lineWeights_Counts[5]) # outerCow -> innerCow

text(x=c(1.5, 1.5, 7.5, 7.5), y=c(3.5, 7.5, 3.5, 7.5), 
     labels=c("innerBadger", "innerCow", "outerBadger", "outerCow"), pos=4, col="red")

title <- paste("Migration Rates ", direction, " in time", sep="")
plot(x=NULL,
     y=NULL,
     xlim=c(1, 10),
     ylim=c(1,10),
     xaxt="n", yaxt="n", xlab="", ylab="", main=title)

arrows(x0=2, y0=4, x1=2, y1=7, code=code, lwd=lineWeights_Rates[1] * factor) # innerBadger -> innerCow
arrows(x0=3, y0=7, x1=3, y1=4, code=code, lwd=lineWeights_Rates[3] * factor) # innerCow -> innerBadger

arrows(x0=8, y0=4, x1=8, y1=7, code=code, lwd=lineWeights_Rates[8] * factor) # outerBadger -> outerCow
arrows(x0=9, y0=7, x1=9, y1=4, code=code, lwd=lineWeights_Rates[6] * factor) # outerCow -> outerBadger

arrows(x0=4, y0=3, x1=7.5, y1=3, code=code, lwd=lineWeights_Rates[2] * factor) # innerBadger -> outerBadger
arrows(x0=7.5, y0=4, x1=4, y1=4, code=code, lwd=lineWeights_Rates[7] * factor) # outerBadger -> innerBadger

arrows(x0=4, y0=7, x1=7.5, y1=7, code=code, lwd=lineWeights_Rates[4] * factor) # innerCow -> outerCow
arrows(x0=7.5, y0=8, x1=4, y1=8, code=code, lwd=lineWeights_Rates[5] * factor) # outerCow -> innerCow

text(x=c(1.5, 1.5, 7.5, 7.5), y=c(3.5, 7.5, 3.5, 7.5), 
     labels=c("innerBadger", "innerCow", "outerBadger", "outerCow"), pos=4, col="red")

####################
# Population sizes #
####################

par(mfrow=c(1,1))

values <- c(logTable$migModel.popSize_0,logTable$migModel.popSize_1,logTable$migModel.popSize_2,logTable$migModel.popSize_3)
xLim <- c(min(values), max(values))

breaks <- seq(from=xLim[1], to=xLim[2] + 2, by=2)

pop0 <- hist(logTable$migModel.popSize_0, plot=FALSE, breaks=breaks)
pop1 <- hist(logTable$migModel.popSize_1, plot=FALSE, breaks=breaks)
pop2 <- hist(logTable$migModel.popSize_2, plot=FALSE, breaks=breaks)
pop3 <- hist(logTable$migModel.popSize_3, plot=FALSE, breaks=breaks)

# Define the Y-axis limits
yLim <- c(0, max(pop0$counts, pop1$counts, pop2$counts, pop3$counts))

# Define the colours
colour0 <- rgb(1,0,0, 0.5) 
colour1 <- rgb(0,0,1, 0.5)
colour2 <- rgb(0,1,0, 0.5) 
colour3 <- rgb(0,0,0, 0.5)

# Plot the histograms
plot(pop0, col=colour0, las=1, main="Deme Effective Population Sizes", xlab="Effective Size",
     ylim=yLim, xlim=xLim)
plot(pop1, col=colour1, add=TRUE)
plot(pop2, col=colour2, add=TRUE)
plot(pop3, col=colour3, add=TRUE)

legend("topright", legend=c("innerBadger", "innerCow", "outerCow", "outerBadger"), bty="n",
       text.col=c(rgb(1,0,0, 1), rgb(0,0,1, 1), rgb(0,1,0, 1), rgb(0,0,0, 1)))

#################
# Mutation Rate #
#################

par(mfrow=c(1,1))
genomeSize <- 700377 + 1322380 + 1316990 + 699795 + 9464
hist(logTable$mutationRate * genomeSize, las=1, main="Mutation Rate", xlab="Rate (per Genome per Year)",
     breaks=50)

plot(logTable$mutationRate * genomeSize, type="l", xlab="", ylab="Rate (per Genome per Year)",
     main="Mutation Rate Trace", las=1)

#################
# Banana Effect #
#################

rbPal <- colorRampPalette(c('red','blue'))
levels(cut(logTable$posterior,breaks = 10))
colours <- rbPal(10)[as.numeric(cut(logTable$posterior,breaks = 10))]

plot(x=logTable$mutationRate * genomeSize, y=logTable$tree.height, pch=20, 
     xlab="Substitution Rate (per Genome per Year)", ylab="Root Height (years)",
     main="Substitution Rate versus Root Height",
     col=colours)

legend("topright", legend=c("High", "Medium", "Low"), pch=20, col=c("red", "purple", "blue"),
       bty="n")

##########################
# Effective Sample Sizes #
##########################

parameterNames <- c()
essValues <- c()
index <- 0

for(column in 2:ncol(logTable)){
  
  if(mean(logTable[, column]) == 0){
    next;
  }
  
  index <- index + 1
  parameterNames[index] <- colnames(logTable)[column]
  essValues[index] <- calculateEffectiveSampleSize(logTable[, column])
}

par(mar=c(0,11,2,0.5)) # bottom, left, top, right
barplot(essValues, las=1, names=parameterNames, horiz=TRUE, xaxt="n", cex.names=0.9,
        main="Parameter Effective Sample Sizes")
abline(v=100, lty=2, col="red")
abline(v=1000, lty=2, col="blue")
text(x=c(250, 1200), y=c(0, 0), labels=c("100", "1000"), cex=0.5, col=c("red", "blue"))

dev.off()

#############
# FUNCTIONS #
#############

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


