################################
# Read in the BASTA output log #
################################

# Read in the log file
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/BASTA/"

type <- "Equal" # Equal or Varying

date <- "12-06-17"

folderName <- paste("6Demes_", type, "PopSizes_", date, "/", sep="")

file <- paste(path, folderName,
              "6Demes_", type, "PopSizes_", date, ".log", sep="")

logTable <- read.table(file, header=TRUE, sep="\t")

# Remove the burn-in
burnIn <- round(0.1 * nrow(logTable), digits=0)
logTable <- logTable[burnIn:nrow(logTable), ]

# Open a PDF
prefix <- paste("6Demes_", type, "PopSizes_", date, sep="")
file <- paste(path, folderName, prefix, "_ResultsSummary_", date, ".pdf", sep="")
pdf(file)

###################
# Migration Rates #
###################

# BadgerInner       0     SAMPLED
# CowInner          1     SAMPLED
# CowOuterEast      2     SAMPLED
# CowOuterWest      3     SAMPLED
# BadgerOuterEast   4     UNSAMPLED
# BadgerOuterWest   5     UNSAMPLED

# Set arrow direction
direction <- "FORWARDS" 
factor <- 250
colour <- rgb(0,0,0, 1)

code <- 2
if(direction == "FORWARDS"){
  code <- 1
}

par(mfrow=c(1,1))

# Define the line style
par(lend=1) # 0 = round, 1 = butt, 2 = square

# Define line weights
lineWeights_Counts <- c(
  mean(logTable[, "treePrior.Count0to1"]), # 1  BadgerInner      -> CowInner
  mean(logTable[, "treePrior.Count0to4"]), # 2  BadgerInner      -> BadgerOuterEast
  mean(logTable[, "treePrior.Count0to5"]), # 3  BadgerInner      -> BadgerOuterWest
  
  mean(logTable[, "treePrior.Count1to0"]), # 4  CowInner         -> BadgerInner
  mean(logTable[, "treePrior.Count1to2"]), # 5  CowInner         -> CowOuterEast
  mean(logTable[, "treePrior.Count1to3"]), # 6  CowInner         -> CowOuterWest
  
  mean(logTable[, "treePrior.Count2to1"]), # 7  CowOuterEast     -> CowInner
  mean(logTable[, "treePrior.Count2to3"]), # 8  CowOuterEast     -> CowOuterWest
  mean(logTable[, "treePrior.Count2to4"]), # 9  CowOuterEast     -> BadgerOuterEast
  
  mean(logTable[, "treePrior.Count3to1"]), # 10 CowOuterWest     -> CowInner
  mean(logTable[, "treePrior.Count3to2"]), # 11 CowOuterWest     -> CowOuterEast
  mean(logTable[, "treePrior.Count3to5"]), # 12 CowOuterWest     -> BadgerOuterWest
  
  mean(logTable[, "treePrior.Count4to0"]), # 13 BadgerOuterEast  -> BadgerInner
  mean(logTable[, "treePrior.Count4to2"]), # 14 BadgerOuterEast  -> CowOuterEast
  mean(logTable[, "treePrior.Count4to5"]), # 15 BadgerOuterEast  -> BadgerOuterWest
  
  mean(logTable[, "treePrior.Count5to0"]), # 16 BadgerOuterWest  -> BadgerInner
  mean(logTable[, "treePrior.Count5to3"]), # 17 BadgerOuterWest  -> CowOuterWest
  mean(logTable[, "treePrior.Count5to4"])  # 18 BadgerOuterWest  -> BadgerOuterEast
)

lineWeights_Rates <- c(
  mean(logTable[, "migModel.rateMatrix_0_1"]), # 1  BadgerInner      -> CowInner
  mean(logTable[, "migModel.rateMatrix_0_4"]), # 2  BadgerInner      -> BadgerOuterEast
  mean(logTable[, "migModel.rateMatrix_0_5"]), # 3  BadgerInner      -> BadgerOuterWest
  
  mean(logTable[, "migModel.rateMatrix_1_0"]), # 4  CowInner         -> BadgerInner
  mean(logTable[, "migModel.rateMatrix_1_2"]), # 5  CowInner         -> CowOuterEast
  mean(logTable[, "migModel.rateMatrix_1_3"]), # 6  CowInner         -> CowOuterWest
  
  mean(logTable[, "migModel.rateMatrix_2_1"]), # 7  CowOuterEast     -> CowInner
  mean(logTable[, "migModel.rateMatrix_2_3"]), # 8  CowOuterEast     -> CowOuterWest
  mean(logTable[, "migModel.rateMatrix_2_4"]), # 9  CowOuterEast     -> BadgerOuterEast
  
  mean(logTable[, "migModel.rateMatrix_3_1"]), # 10 CowOuterWest     -> CowInner
  mean(logTable[, "migModel.rateMatrix_3_2"]), # 11 CowOuterWest     -> CowOuterEast
  mean(logTable[, "migModel.rateMatrix_3_5"]), # 12 CowOuterWest     -> BadgerOuterWest
  
  mean(logTable[, "migModel.rateMatrix_4_0"]), # 13 BadgerOuterEast  -> BadgerInner
  mean(logTable[, "migModel.rateMatrix_4_2"]), # 14 BadgerOuterEast  -> CowOuterEast
  mean(logTable[, "migModel.rateMatrix_4_5"]), # 15 BadgerOuterEast  -> BadgerOuterWest
  
  mean(logTable[, "migModel.rateMatrix_5_0"]), # 16 BadgerOuterWest  -> BadgerInner
  mean(logTable[, "migModel.rateMatrix_5_3"]), # 17 BadgerOuterWest  -> CowOuterWest
  mean(logTable[, "migModel.rateMatrix_5_4"])  # 18 BadgerOuterWest  -> BadgerOuterEast
)

## Migration Counts

title <- paste("Migration Counts ", direction, " in time", sep="")
par(mar=c(0,0,4.1,0)) # Bottom, Left, Top, Right
plot(x=NULL,
     y=NULL,
     xlim=c(0,10),
     ylim=c(0,10),
     xaxt="n", yaxt="n", xlab="", ylab="", main=title, bty="n")

text(x=c(4, 4, 0.5, 7.5, 7.5, 0.5), 
     y=c(4, 6, 8, 8, 2, 2), 
     labels=c("BadgerInner", "CowInner", "CowOuterWest", "CowOuterEast",
              "BadgerOuterEast", "BadgerOuterWest"),
     pos=4, col=c("red", "blue", "blue", "blue", "red", "red"), cex=1)

arrows(x0=4.5, y0=4.5, x1=4.5, y1=5.5, col=colour,
       code=code, lwd=lineWeights_Counts[1])    # 1  BadgerInner      -> CowInner
arrows(x0=6.5, y0=3.5, x1=8, y1=2.5, col=colour,
       code=code, lwd=lineWeights_Counts[3])    # 2  BadgerInner      -> BadgerOuterEast
arrows(x0=4, y0=3.5, x1=2.5, y1=2.5, col=colour,
       code=code, lwd=lineWeights_Counts[2])    # 3  BadgerInner      -> BadgerOuterWest

arrows(x0=5.5, y0=5.5, x1=5.5, y1=4.5, col=colour,
       code=code, lwd=lineWeights_Counts[4])    # 4  CowInner         -> BadgerInner
arrows(x0=6.5, y0=6.5, x1=8, y1=7.5, col=colour,
       code=code, lwd=lineWeights_Counts[5])    # 5  CowInner         -> CowOuterEast
arrows(x0=4, y0=6.5, x1=2.5, y1=7.5, col=colour,
       code=code, lwd=lineWeights_Counts[6])    # 6  CowInner         -> CowOuterWest

arrows(x0=7, y0=7.5, x1=5.5, y1=6.5, col=colour,
       code=code, lwd=lineWeights_Counts[5])    # 7  CowOuterEast     -> CowInner
arrows(x0=7, y0=9, x1=3.5, y1=9, col=colour,
       code=code, lwd=lineWeights_Counts[8])    # 8  CowOuterEast     -> CowOuterWest
arrows(x0=9.5, y0=7.5, x1=9.5, y1=2.5, col=colour,
       code=code, lwd=lineWeights_Counts[9])    # 9  CowOuterEast     -> BadgerOuterEast

arrows(x0=3.5, y0=7.5, x1=5, y1=6.5, col=colour,
       code=code, lwd=lineWeights_Counts[6])   # 10 CowOuterWest     -> CowInner
arrows(x0=3.5, y0=8, x1=7, y1=8, col=colour,
       code=code, lwd=lineWeights_Counts[8])   # 11 CowOuterWest     -> CowOuterEast
arrows(x0=1, y0=7.5, x1=1, y1=2.5, col=colour,
       code=code, lwd=lineWeights_Counts[9])   # 12 CowOuterWest     -> BadgerOuterWest

arrows(x0=7, y0=2.5, x1=5.5, y1=3.5, col=colour,
       code=code, lwd=lineWeights_Counts[3])   # 13 BadgerOuterEast  -> BadgerInner
arrows(x0=8.5, y0=2.5, x1=8.5, y1=7.5, col=colour,
       code=code, lwd=lineWeights_Counts[9])    # 14 BadgerOuterEast  -> CowOuterEast
arrows(x0=7, y0=1, x1=3.5, y1=1, col=colour,
       code=code, lwd=lineWeights_Counts[8])    # 15 BadgerOuterEast  -> BadgerOuterWest

arrows(x0=3.5, y0=2.5, x1=5, y1=3.5, col=colour,
       code=code, lwd=lineWeights_Counts[2])   # 16 BadgerOuterWest  -> BadgerInner
arrows(x0=2, y0=2.5, x1=2, y1=7.5, col=colour,
       code=code, lwd=lineWeights_Counts[9])  # 17 BadgerOuterWest  -> CowOuterWest
arrows(x0=3.5, y0=2, x1=7, y1=2, col=colour,
       code=code, lwd=lineWeights_Counts[8])   # 18 BadgerOuterWest  -> BadgerOuterEast

### Migration Rates
title <- paste("Migration Rates ", direction, " in time", sep="")
par(mar=c(0,0,4.1,0)) # Bottom, Left, Top, Right
plot(x=NULL,
     y=NULL,
     xlim=c(0,10),
     ylim=c(0,10),
     xaxt="n", yaxt="n", xlab="", ylab="", main=title, bty="n")

text(x=c(4, 4, 0.5, 7.5, 7.5, 0.5), 
     y=c(4, 6, 8, 8, 2, 2), 
     labels=c("BadgerInner", "CowInner", "CowOuterWest", "CowOuterEast",
              "BadgerOuterEast", "BadgerOuterWest"),
     pos=4, col=c("red", "blue", "blue", "blue", "red", "red"), cex=1)

arrows(x0=4.5, y0=4.5, x1=4.5, y1=5.5, col=colour,
       code=code, lwd=lineWeights_Rates[1] * factor)    # 1  BadgerInner      -> CowInner
arrows(x0=6.5, y0=3.5, x1=8, y1=2.5, col=colour,
       code=code, lwd=lineWeights_Rates[3] * factor)    # 2  BadgerInner      -> BadgerOuterEast
arrows(x0=4, y0=3.5, x1=2.5, y1=2.5, col=colour,
       code=code, lwd=lineWeights_Rates[2] * factor)    # 3  BadgerInner      -> BadgerOuterWest

arrows(x0=5.5, y0=5.5, x1=5.5, y1=4.5, col=colour,
       code=code, lwd=lineWeights_Rates[4] * factor)    # 4  CowInner         -> BadgerInner
arrows(x0=6.5, y0=6.5, x1=8, y1=7.5, col=colour,
       code=code, lwd=lineWeights_Rates[5] * factor)    # 5  CowInner         -> CowOuterEast
arrows(x0=4, y0=6.5, x1=2.5, y1=7.5, col=colour,
       code=code, lwd=lineWeights_Rates[6] * factor)    # 6  CowInner         -> CowOuterWest

arrows(x0=7, y0=7.5, x1=5.5, y1=6.5, col=colour,
       code=code, lwd=lineWeights_Rates[5] * factor)    # 7  CowOuterEast     -> CowInner
arrows(x0=7, y0=9, x1=3.5, y1=9, col=colour,
       code=code, lwd=lineWeights_Rates[8] * factor)    # 8  CowOuterEast     -> CowOuterWest
arrows(x0=9.5, y0=7.5, x1=9.5, y1=2.5, col=colour,
       code=code, lwd=lineWeights_Rates[9] * factor)    # 9  CowOuterEast     -> BadgerOuterEast

arrows(x0=3.5, y0=7.5, x1=5, y1=6.5, col=colour,
       code=code, lwd=lineWeights_Rates[6] * factor)   # 10 CowOuterWest     -> CowInner
arrows(x0=3.5, y0=8, x1=7, y1=8, col=colour,
       code=code, lwd=lineWeights_Rates[8] * factor)   # 11 CowOuterWest     -> CowOuterEast
arrows(x0=1, y0=7.5, x1=1, y1=2.5, col=colour,
       code=code, lwd=lineWeights_Rates[9] * factor)   # 12 CowOuterWest     -> BadgerOuterWest

arrows(x0=7, y0=2.5, x1=5.5, y1=3.5, col=colour,
       code=code, lwd=lineWeights_Rates[3] * factor)   # 13 BadgerOuterEast  -> BadgerInner
arrows(x0=8.5, y0=2.5, x1=8.5, y1=7.5, col=colour,
       code=code, lwd=lineWeights_Rates[9] * factor)    # 14 BadgerOuterEast  -> CowOuterEast
arrows(x0=7, y0=1, x1=3.5, y1=1, col=colour,
       code=code, lwd=lineWeights_Rates[8] * factor)    # 15 BadgerOuterEast  -> BadgerOuterWest

arrows(x0=3.5, y0=2.5, x1=5, y1=3.5, col=colour,
       code=code, lwd=lineWeights_Rates[2] * factor)   # 16 BadgerOuterWest  -> BadgerInner
arrows(x0=2, y0=2.5, x1=2, y1=7.5, col=colour,
       code=code, lwd=lineWeights_Rates[9] * factor)  # 17 BadgerOuterWest  -> CowOuterWest
arrows(x0=3.5, y0=2, x1=7, y1=2, col=colour,
       code=code, lwd=lineWeights_Rates[8] * factor)   # 18 BadgerOuterWest  -> BadgerOuterEast

####################
# Population sizes #
####################

par(mfrow=c(1,1))
par(mar=c(5.4, 4.1, 4.1, 2.1)) # Bottom, Left, Top, Right

values <- c(logTable$migModel.popSize_0,
            logTable$migModel.popSize_1,
            logTable$migModel.popSize_2,
            logTable$migModel.popSize_3,
            logTable$migModel.popSize_4,
            logTable$migModel.popSize_5)
xLim <- c(min(values), max(values))

breaks <- seq(from=xLim[1], to=xLim[2] + 2, by=2)

pop0 <- hist(logTable$migModel.popSize_0, plot=FALSE, breaks=breaks)
pop1 <- hist(logTable$migModel.popSize_1, plot=FALSE, breaks=breaks)
pop2 <- hist(logTable$migModel.popSize_2, plot=FALSE, breaks=breaks)
pop3 <- hist(logTable$migModel.popSize_3, plot=FALSE, breaks=breaks)
pop4 <- hist(logTable$migModel.popSize_4, plot=FALSE, breaks=breaks)
pop5 <- hist(logTable$migModel.popSize_5, plot=FALSE, breaks=breaks)

# Define the Y-axis limits
yLim <- c(0, max(pop0$counts, pop1$counts, pop2$counts, pop3$counts,
                 pop4$counts, pop5$counts))

# Define the colours
alpha <- 0.5
colour0 <- rgb(1,0,0, alpha) 
colour1 <- rgb(0,0,1, alpha)
colour2 <- rgb(0,1,0, alpha) 
colour3 <- rgb(0,0,0, alpha)
colour4 <- rgb(1,1,0, alpha)
colour5 <- rgb(0,1,1, alpha)

# Plot the histograms
plot(pop0, col=colour0, las=1, main="Deme Effective Population Sizes", xlab="Effective Size",
     ylim=yLim, xlim=xLim)
plot(pop1, col=colour1, add=TRUE)
plot(pop2, col=colour2, add=TRUE)
plot(pop3, col=colour3, add=TRUE)
plot(pop4, col=colour4, add=TRUE)
plot(pop5, col=colour5, add=TRUE)

alpha <- 1
colour0 <- rgb(1,0,0, alpha) 
colour1 <- rgb(0,0,1, alpha)
colour2 <- rgb(0,1,0, alpha) 
colour3 <- rgb(0,0,0, alpha)
colour4 <- rgb(1,1,0, alpha)
colour5 <- rgb(0,1,1, alpha)
legend("topright",  bty="n",
       legend=c("BadgerInner", "CowInner", "CowOuterWest", "CowOuterEast",
                "BadgerOuterEast", "BadgerOuterWest"),
       text.col=c(colour0, colour1, colour2, colour3,
                  colour4, colour5))

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
barplot(essValues, las=1, names=parameterNames, horiz=TRUE, xaxt="n", cex.names=0.5,
        main="Parameter Effective Sample Sizes")
abline(v=100, lty=2, col="red")
abline(v=1000, lty=2, col="blue")
text(x=c(250, 1200), y=c(0, 0), labels=c("100", "1000"), cex=0.5, col=c("red", "blue"))

########################
# Plot Root Allocation #
########################

# Reset Margins
par(mar=c(3,0,3,0))

rootStateCounts <- table(logTable$treePrior.rootColor)

barplot(rootStateCounts, yaxt="n", 
        names=c("BadgerInner", "CowInner", "CowOuterWest", "CowOuterEast",
                "BadgerOuterEast", "BadgerOuterWest"),
        main="Assignment of Demes to Root State", cex.names=0.5)


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


