# Read in the log file
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/BASTA/"

type <- "Equal"

folderName <- paste("8Demes_", type, "PopSizes_03-03-17/", sep="")

file <- paste(path, folderName,
              "8Demes_", type, "PopSizes_03-03-17.log", sep="")

logTable <- read.table(file, header=TRUE, sep="\t")

# Remove the burn-in
burnIn <- round(0.1 * nrow(logTable), digits=0)
logTable <- logTable[burnIn:nrow(logTable), ]

# Open a PDF
prefix <- paste("8Demes_", type, "PopSizes_07-02-17", sep="")
file <- paste(path, folderName, prefix, "_ResultsSummary_09-02-17.pdf", sep="")
pdf(file)

###################
# Migration Rates #
###################

# BadgerInnerEast   0     SAMPLED
# BadgerInnerWest   1     SAMPLED
# CowInnerEast      2     SAMPLED
# CowInnerWest      3     SAMPLED
# CowOuterEast      4     SAMPLED
# CowOuterWest      5     SAMPLED
# BadgerOuterEast   6     UNSAMPLED
# BadgerOuterWest   7     UNSAMPLED

par(mfrow=c(1,1))

# Define line weights
lineWeights_Counts <- c(
  mean(logTable[, "treePrior.Count0to1"]), # BadgerInnerEast  -> BadgerInnerWest
  mean(logTable[, "treePrior.Count0to2"]), # BadgerInnerEast  -> CowInnerEast
  mean(logTable[, "treePrior.Count0to6"]), # BadgerInnerEast  -> BadgerOuterEast
  mean(logTable[, "treePrior.Count1to0"]), # BadgerInnerWest  -> BadgerInnerEast
  mean(logTable[, "treePrior.Count1to3"]), # BadgerInnerWest  -> CowInnerWest
  mean(logTable[, "treePrior.Count1to7"]), # BadgerInnerWest  -> BadgerOuterWest
  mean(logTable[, "treePrior.Count2to0"]), # CowInnerEast     -> BadgerInnerEast
  mean(logTable[, "treePrior.Count2to3"]), # CowInnerEast     -> CowInnerWest
  mean(logTable[, "treePrior.Count2to4"]), # CowInnerEast     -> CowOuterEast
  mean(logTable[, "treePrior.Count3to1"]), # CowInnerWest     -> BadgerInnerWest
  mean(logTable[, "treePrior.Count3to2"]), # CowInnerWest     -> CowInnerEast
  mean(logTable[, "treePrior.Count3to5"]), # CowInnerWest     -> CowOuterWest
  mean(logTable[, "treePrior.Count4to2"]), # CowOuterEast     -> CowInnerEast
  mean(logTable[, "treePrior.Count4to5"]), # CowOuterEast     -> CowOuterWest
  mean(logTable[, "treePrior.Count4to6"]), # CowOuterEast     -> BadgerOuterEast
  mean(logTable[, "treePrior.Count5to3"]), # CowOuterWest     -> CowInnerWest
  mean(logTable[, "treePrior.Count5to4"]), # CowOuterWest     -> CowOuterEast
  mean(logTable[, "treePrior.Count5to7"]), # CowOuterWest     -> BadgerOuterWest
  mean(logTable[, "treePrior.Count6to0"]), # BadgerOuterEast  -> BadgerInnerEast
  mean(logTable[, "treePrior.Count6to4"]), # BadgerOuterEast  -> CowOuterEast
  mean(logTable[, "treePrior.Count6to7"]), # BadgerOuterEast  -> BadgerOuterWest
  mean(logTable[, "treePrior.Count7to1"]), # BadgerOuterWest  -> BadgerInnerWest
  mean(logTable[, "treePrior.Count7to5"]), # BadgerOuterWest  -> CowOuterWest
  mean(logTable[, "treePrior.Count7to6"])  # BadgerOuterWest  -> BadgerOuterEast
)

lineWeights_Rates <- c(
  mean(logTable[, "migModel.rateMatrix_0_1"]), # 1  BadgerInnerEast  -> BadgerInnerWest
  mean(logTable[, "migModel.rateMatrix_0_2"]), # 2  BadgerInnerEast  -> CowInnerEast
  mean(logTable[, "migModel.rateMatrix_0_6"]), # 3  BadgerInnerEast  -> BadgerOuterEast
  mean(logTable[, "migModel.rateMatrix_1_0"]), # 4  BadgerInnerWest  -> BadgerInnerEast
  mean(logTable[, "migModel.rateMatrix_1_3"]), # 5  BadgerInnerWest  -> CowInnerWest
  mean(logTable[, "migModel.rateMatrix_1_7"]), # 6  BadgerInnerWest  -> BadgerOuterWest
  mean(logTable[, "migModel.rateMatrix_2_0"]), # 7  CowInnerEast     -> BadgerInnerEast
  mean(logTable[, "migModel.rateMatrix_2_3"]), # 8  CowInnerEast     -> CowInnerWest
  mean(logTable[, "migModel.rateMatrix_2_4"]), # 9  CowInnerEast     -> CowOuterEast
  mean(logTable[, "migModel.rateMatrix_3_1"]), # 10 CowInnerWest     -> BadgerInnerWest
  mean(logTable[, "migModel.rateMatrix_3_2"]), # 11 CowInnerWest     -> CowInnerEast
  mean(logTable[, "migModel.rateMatrix_3_5"]), # 12 CowInnerWest     -> CowOuterWest
  mean(logTable[, "migModel.rateMatrix_4_2"]), # 13 CowOuterEast     -> CowInnerEast
  mean(logTable[, "migModel.rateMatrix_4_5"]), # 14 CowOuterEast     -> CowOuterWest
  mean(logTable[, "migModel.rateMatrix_4_6"]), # 15 CowOuterEast     -> BadgerOuterEast
  mean(logTable[, "migModel.rateMatrix_5_3"]), # 16 CowOuterWest     -> CowInnerWest
  mean(logTable[, "migModel.rateMatrix_5_4"]), # 17 CowOuterWest     -> CowOuterEast
  mean(logTable[, "migModel.rateMatrix_5_7"]), # 18 CowOuterWest     -> BadgerOuterWest
  mean(logTable[, "migModel.rateMatrix_6_0"]), # 19 BadgerOuterEast  -> BadgerInnerEast
  mean(logTable[, "migModel.rateMatrix_6_4"]), # 20 BadgerOuterEast  -> CowOuterEast
  mean(logTable[, "migModel.rateMatrix_6_7"]), # 21 BadgerOuterEast  -> BadgerOuterWest
  mean(logTable[, "migModel.rateMatrix_7_1"]), # 22 BadgerOuterWest  -> BadgerInnerWest
  mean(logTable[, "migModel.rateMatrix_7_5"]), # 23 BadgerOuterWest  -> CowOuterWest
  mean(logTable[, "migModel.rateMatrix_7_6"])  # 24 BadgerOuterWest  -> BadgerOuterEast
)

# Set arrow direction
code <- 2 # Backward in time (1 for forward)
factor <- 250
colour <- rgb(0,0,0, 1)

direction <- "BACKWARDS"
if(code == 1){
  direction <- "FORWARDS"
}

## Migration Counts

title <- paste("Migration Counts ", direction, " in time", sep="")
plot(x=NULL,
     y=NULL,
     xlim=c(1, 10),
     ylim=c(1,10),
     xaxt="n", yaxt="n", xlab="", ylab="", main=title)

arrows(x0=6, y0=9, x1=5, y1=9, col=colour,
       code=code, lwd=lineWeights_Counts[1])    # 1  BadgerInnerEast  -> BadgerInnerWest
arrows(x0=6.5, y0=8, x1=6.5, y1=3, col=colour,
       code=code, lwd=lineWeights_Counts[2])    # 2  BadgerInnerEast  -> CowInnerEast
arrows(x0=7.5, y0=8, x1=8, y1=7, col=colour,
       code=code, lwd=lineWeights_Counts[3])    # 3  BadgerInnerEast  -> BadgerOuterEast

arrows(x0=5, y0=8, x1=6, y1=8, col=colour,
       code=code, lwd=lineWeights_Counts[4])    # 4  BadgerInnerWest  -> BadgerInnerEast
arrows(x0=4.5, y0=8, x1=4.5, y1=3, col=colour,
       code=code, lwd=lineWeights_Counts[5])    # 5  BadgerInnerWest  -> CowInnerWest
arrows(x0=3, y0=8, x1=2.5, y1=7, col=colour,
       code=code, lwd=lineWeights_Counts[6])    # 6  BadgerInnerWest  -> BadgerOuterWest

arrows(x0=7, y0=3, x1=7, y1=8, col=colour,
       code=code, lwd=lineWeights_Counts[7])    # 7  CowInnerEast     -> BadgerInnerEast
arrows(x0=6, y0=3, x1=5, y1=3, col=colour,
       code=code, lwd=lineWeights_Counts[8])    # 8  CowInnerEast     -> CowInnerWest
arrows(x0=7.5, y0=3, x1=8, y1=4, col=colour,
       code=code, lwd=lineWeights_Counts[9])    # 9  CowInnerEast     -> CowOuterEast

arrows(x0=4, y0=3, x1=4, y1=8, col=colour,
       code=code, lwd=lineWeights_Counts[10])   # 10 CowInnerWest     -> BadgerInnerWest
arrows(x0=5, y0=2, x1=6, y1=2, col=colour,
       code=code, lwd=lineWeights_Counts[11])   # 11 CowInnerWest     -> CowInnerEast
arrows(x0=3, y0=3, x1=2.5, y1=4, col=colour,
       code=code, lwd=lineWeights_Counts[12])   # 12 CowInnerWest     -> CowOuterWest

arrows(x0=9, y0=4, x1=8.5, y1=3, col=colour,
       code=code, lwd=lineWeights_Counts[13])   # 13 CowOuterEast     -> CowInnerEast
arrows(x0=8, y0=5, x1=3, y1=5, col=colour,
       code=code, lwd=lineWeights_Counts[14])   # 14 CowOuterEast     -> CowOuterWest
arrows(x0=8.5, y0=5, x1=8.5, y1=6, col=colour,
       code=code, lwd=lineWeights_Counts[15])   # 15 CowOuterEast     -> BadgerOuterEast

arrows(x0=1.5, y0=4, x1=2, y1=3, col=colour,
       code=code, lwd=lineWeights_Counts[16])   # 16 CowOuterWest     -> CowInnerWest
arrows(x0=3, y0=4.5, x1=8, y1=4.5, col=colour,
       code=code, lwd=lineWeights_Counts[17])   # 17 CowOuterWest     -> CowOuterEast
arrows(x0=2, y0=5, x1=2, y1=6, col=colour,
       code=code, lwd=lineWeights_Counts[18])   # 18 CowOuterWest     -> BadgerOuterWest

arrows(x0=9, y0=7, x1=8.5, y1=8, col=colour,
       code=code, lwd=lineWeights_Counts[19])   # 19 BadgerOuterEast  -> BadgerInnerEast
arrows(x0=9.5, y0=6, x1=9.5, y1=5, col=colour,
       code=code, lwd=lineWeights_Counts[20])   # 20 BadgerOuterEast  -> CowOuterEast
arrows(x0=8, y0=6, x1=3, y1=6, col=colour,
       code=code, lwd=lineWeights_Counts[21])   # 21 BadgerOuterEast  -> BadgerOuterWest

arrows(x0=1.5, y0=7, x1=2, y1=8, col=colour,
       code=code, lwd=lineWeights_Counts[22])   # 22 BadgerOuterWest  -> BadgerInnerWest
arrows(x0=1, y0=6, x1=1, y1=5, col=colour,
       code=code, lwd=lineWeights_Counts[23])   # 23 BadgerOuterWest  -> CowOuterWest
arrows(x0=3, y0=6.5, x1=8, y1=6.5, col=colour,
       code=code, lwd=lineWeights_Counts[24])   # 24 BadgerOuterWest  -> BadgerOuterEast

text(x=c(2.5, 6, 8, 8, 6, 2.5, 0.5, 0.5), 
     y=c(8.5, 8.5, 6.5, 4.5, 2.5, 2.5, 4.5, 6.5), 
     labels=c("BadgerInnerWest", "BadgerInnerEast", "BadgerOuterEast",
              "CowOuterEast", "CowInnerEast", "CowInnerWest", "CowOuterWest",
              "BadgerOuterWest"),
     pos=4, col="red", cex=1)


### Migration Rates

title <- paste("Migration Rates ", direction, " in time", sep="")
plot(x=NULL,
     y=NULL,
     xlim=c(1, 10),
     ylim=c(1,10),
     xaxt="n", yaxt="n", xlab="", ylab="", main=title)

arrows(x0=6, y0=9, x1=5, y1=9, col=colour,
       code=code, lwd=lineWeights_Rates[1] * factor)    # 1  BadgerInnerEast  -> BadgerInnerWest
arrows(x0=6.5, y0=8, x1=6.5, y1=3, col=colour,
       code=code, lwd=lineWeights_Rates[2] * factor)    # 2  BadgerInnerEast  -> CowInnerEast
arrows(x0=7.5, y0=8, x1=8, y1=7, col=colour,
       code=code, lwd=lineWeights_Rates[3] * factor)    # 3  BadgerInnerEast  -> BadgerOuterEast

arrows(x0=5, y0=8, x1=6, y1=8, col=colour,
       code=code, lwd=lineWeights_Rates[4] * factor)    # 4  BadgerInnerWest  -> BadgerInnerEast
arrows(x0=4.5, y0=8, x1=4.5, y1=3, col=colour,
       code=code, lwd=lineWeights_Rates[5] * factor)    # 5  BadgerInnerWest  -> CowInnerWest
arrows(x0=3, y0=8, x1=2.5, y1=7, col=colour,
       code=code, lwd=lineWeights_Rates[6] * factor)    # 6  BadgerInnerWest  -> BadgerOuterWest

arrows(x0=7, y0=3, x1=7, y1=8, col=colour,
       code=code, lwd=lineWeights_Rates[7] * factor)    # 7  CowInnerEast     -> BadgerInnerEast
arrows(x0=6, y0=3, x1=5, y1=3, col=colour,
       code=code, lwd=lineWeights_Rates[8] * factor)    # 8  CowInnerEast     -> CowInnerWest
arrows(x0=7.5, y0=3, x1=8, y1=4, col=colour,
       code=code, lwd=lineWeights_Rates[9] * factor)    # 9  CowInnerEast     -> CowOuterEast

arrows(x0=4, y0=3, x1=4, y1=8, col=colour,
       code=code, lwd=lineWeights_Rates[10] * factor)   # 10 CowInnerWest     -> BadgerInnerWest
arrows(x0=5, y0=2, x1=6, y1=2, col=colour,
       code=code, lwd=lineWeights_Rates[11] * factor)   # 11 CowInnerWest     -> CowInnerEast
arrows(x0=3, y0=3, x1=2.5, y1=4, col=colour,
       code=code, lwd=lineWeights_Rates[12] * factor)   # 12 CowInnerWest     -> CowOuterWest

arrows(x0=9, y0=4, x1=8.5, y1=3, col=colour,
       code=code, lwd=lineWeights_Rates[13] * factor)   # 13 CowOuterEast     -> CowInnerEast
arrows(x0=8, y0=5, x1=3, y1=5, col=colour,
       code=code, lwd=lineWeights_Rates[14] * factor)   # 14 CowOuterEast     -> CowOuterWest
arrows(x0=8.5, y0=5, x1=8.5, y1=6, col=colour,
       code=code, lwd=lineWeights_Rates[15] * factor)   # 15 CowOuterEast     -> BadgerOuterEast

arrows(x0=1.5, y0=4, x1=2, y1=3, col=colour,
       code=code, lwd=lineWeights_Rates[16] * factor)   # 16 CowOuterWest     -> CowInnerWest
arrows(x0=3, y0=4.5, x1=8, y1=4.5, col=colour,
       code=code, lwd=lineWeights_Rates[17] * factor)   # 17 CowOuterWest     -> CowOuterEast
arrows(x0=2, y0=5, x1=2, y1=6, col=colour,
       code=code, lwd=lineWeights_Rates[18] * factor)   # 18 CowOuterWest     -> BadgerOuterWest

arrows(x0=9, y0=7, x1=8.5, y1=8, col=colour,
       code=code, lwd=lineWeights_Rates[19] * factor)   # 19 BadgerOuterEast  -> BadgerInnerEast
arrows(x0=9.5, y0=6, x1=9.5, y1=5, col=colour,
       code=code, lwd=lineWeights_Rates[20] * factor)   # 20 BadgerOuterEast  -> CowOuterEast
arrows(x0=8, y0=6, x1=3, y1=6, col=colour,
       code=code, lwd=lineWeights_Rates[21] * factor)   # 21 BadgerOuterEast  -> BadgerOuterWest

arrows(x0=1.5, y0=7, x1=2, y1=8, col=colour,
       code=code, lwd=lineWeights_Rates[22] * factor)   # 22 BadgerOuterWest  -> BadgerInnerWest
arrows(x0=1, y0=6, x1=1, y1=5, col=colour,
       code=code, lwd=lineWeights_Rates[23] * factor)   # 23 BadgerOuterWest  -> CowOuterWest
arrows(x0=3, y0=6.5, x1=8, y1=6.5, col=colour,
       code=code, lwd=lineWeights_Rates[24] * factor)   # 24 BadgerOuterWest  -> BadgerOuterEast

text(x=c(2.5, 6, 8, 8, 6, 2.5, 0.5, 0.5), 
     y=c(8.5, 8.5, 6.5, 4.5, 2.5, 2.5, 4.5, 6.5), 
     labels=c("BadgerInnerWest", "BadgerInnerEast", "BadgerOuterEast",
              "CowOuterEast", "CowInnerEast", "CowInnerWest", "CowOuterWest",
              "BadgerOuterWest"),
     pos=4, col="red", cex=1)

####################
# Population sizes #
####################

par(mfrow=c(1,1))

values <- c(logTable$migModel.popSize_0,
            logTable$migModel.popSize_1,
            logTable$migModel.popSize_2,
            logTable$migModel.popSize_3,
            logTable$migModel.popSize_4,
            logTable$migModel.popSize_5,
            logTable$migModel.popSize_6,
            logTable$migModel.popSize_7)
xLim <- c(min(values), max(values))

breaks <- seq(from=xLim[1], to=xLim[2] + 2, by=2)

pop0 <- hist(logTable$migModel.popSize_0, plot=FALSE, breaks=breaks)
pop1 <- hist(logTable$migModel.popSize_1, plot=FALSE, breaks=breaks)
pop2 <- hist(logTable$migModel.popSize_2, plot=FALSE, breaks=breaks)
pop3 <- hist(logTable$migModel.popSize_3, plot=FALSE, breaks=breaks)
pop4 <- hist(logTable$migModel.popSize_4, plot=FALSE, breaks=breaks)
pop5 <- hist(logTable$migModel.popSize_5, plot=FALSE, breaks=breaks)
pop6 <- hist(logTable$migModel.popSize_6, plot=FALSE, breaks=breaks)
pop7 <- hist(logTable$migModel.popSize_7, plot=FALSE, breaks=breaks)

# Define the Y-axis limits
yLim <- c(0, max(pop0$counts, pop1$counts, pop2$counts, pop3$counts,
                 pop4$counts, pop5$counts, pop6$counts, pop7$counts))

# Define the colours
alpha <- 0.5
colour0 <- rgb(1,0,0, alpha) 
colour1 <- rgb(0,0,1, alpha)
colour2 <- rgb(0,1,0, alpha) 
colour3 <- rgb(0,0,0, alpha)
colour4 <- rgb(1,1,0, alpha)
colour5 <- rgb(0,1,1, alpha)
colour6 <- rgb(1,0,1, alpha)
colour7 <- rgb(1,0.5,0.5, alpha)

# Plot the histograms
plot(pop0, col=colour0, las=1, main="Deme Effective Population Sizes", xlab="Effective Size",
     ylim=yLim, xlim=xLim)
plot(pop1, col=colour1, add=TRUE)
plot(pop2, col=colour2, add=TRUE)
plot(pop3, col=colour3, add=TRUE)
plot(pop4, col=colour4, add=TRUE)
plot(pop5, col=colour5, add=TRUE)
plot(pop6, col=colour6, add=TRUE)
plot(pop7, col=colour7, add=TRUE)


alpha <- 1
colour0 <- rgb(1,0,0, alpha) 
colour1 <- rgb(0,0,1, alpha)
colour2 <- rgb(0,1,0, alpha) 
colour3 <- rgb(0,0,0, alpha)
colour4 <- rgb(1,1,0, alpha)
colour5 <- rgb(0,1,1, alpha)
colour6 <- rgb(1,0,1, alpha)
colour7 <- rgb(1,0.5,0.5, alpha)
legend("topright",  bty="n",
       legend=c("BadgerInnerEast", "BadgerInnerWest", "CowInnerEast",
                "CowInnerWest", "CowOuterEast", "CowOuterWest",
                "BadgerOuterEast", "BadgerOuterWest"),
       text.col=c(colour0, colour1, colour2, colour3,
                  colour4, colour5, colour6, colour7))

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
        names=c("BadgerInnerEast", "BadgerInnerWest", "CowInnerEast",
                "CowInnerWest", "CowOuterEast", "CowOuterWest",
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


