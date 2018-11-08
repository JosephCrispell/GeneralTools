##########
# Set up #
##########

# Set the path
path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/"

#################################
# Note the Genome Size examined #
#################################

# Get the constant site counts
constantSiteCountsFile <- paste(path, "vcfFiles/", "constantSiteCounts_24-03-2018.txt",
                                sep="")
constantSiteCounts <- getConstantSiteCounts(constantSiteCountsFile)
genomeSize <- sum(constantSiteCounts)

# Get number of sites used in FASTA file
fastaFile <- paste(path, "vcfFiles/", "sequences_withoutHomoplasies_27-03-18.fasta", sep="")
nSites <- getNSitesInFASTA(fastaFile)
genomeSize <- genomeSize + nSites

#####################
# Read in the files #
#####################

# Move the path
path <- paste(path, "BASTA/", sep="")

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

# Store each of the log tables in a list
logTables <- readInBASTALogTables(date, demeStructures, popEstimationTypes,
                                  clockEstimateTypes, path, nReplicates=3)

####################
# Examine each run #
####################

# Create summary plots, note migration rates and calculate AICM scores for each model
nBootstraps <- 1000
migrationRateEstimates <- summarisePosteriorLogTables(path, logTables, code=2,
                                                      arrowFactor=20, nBootstraps,
                                                      genomeSize, date)

#################################################################################################
# Examine the transition counts recorded on the posterior tree distributions from each analysis #
#################################################################################################

# Read in the transition counts for all the trees
#   Counts created by Java: MyWork:examineBASTAPosterior:CountTransitions.java
file <- paste0(path, "TransitionCounts_", date, ".txt")
transitionCounts <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep="\t")

# Remove the burn-in period from each analysis
transitionCountsWithoutBurnIn <- removeBurnIn(transitionCounts, burnInProp=0.1)

# Note each model's AICM score
aicmScores <- getAICMScores(migrationRateEstimates)

# Sample the transition counts and associated branch length sums based on the AICM model scores
weightedSampleOfTransitionRates <- getWeightedSampleOfTransitionCountsBasedOnAICMScores(transitionCountsWithoutBurnIn,
                                                                                        aicmScores)

#################################
# Plot Model comparison results #
#################################

# Open a PDF
file <- paste(path, "SummaryFiguresOfModelEstimations_", date, ".pdf", sep="")
pdf(file, height=11, width=10)

# Define the plotting window layout
layout(matrix(c(1,1,1,2,2,2,
                1,1,1,2,2,2,
                1,1,1,2,2,2,
                4,4,3,3,4,4,
                4,4,3,3,4,4), nrow=5, ncol=6, byrow=TRUE))

# Examine the model likelihoods
plotModelAICMScores(migrationRateEstimates, nBootstraps)

######################################################################
# Produce single summary plot of rate estimation: BADGERS <-> CATTLE #
######################################################################

# Calculate the weighted mean estimated transition rates between cattle and badger demes
# weighted by the AICM model scores
weightedMeanEstimates <- 
  calculateMeanEstimatedTransitionRatesBetweenCattleAndBadgerPopulationsWeightedByAICM(
    migrationRateEstimates)


##################################################################
# Plot a summary of the transition counts on the posterior trees #
##################################################################
plotSummaryOfTransitionRatesBasedOnPosteriorTrees(weightedSampleOfTransitionRates)

# Reset the plotting window dimensions
par(mfrow=c(1,1))

# Plot the rates
# plotTransitionRatesBetweenBadgerAndCow(badgerToCow=weightedMeanEstimates[1],
#                                        cowToBadger=weightedMeanEstimates[2],
#                                        code=2, arrowFactor=20)


# Close the pdf
dev.off()

#############
# FUNCTIONS #
#############

plotSummaryOfTransitionRatesBasedOnPosteriorTrees <- function(weightedSampleOfTransitionRates){
  
  # Set the margin sizes
  currentMar <- par("mar")
  par(mar=c(10, 5.5, 4, 0.1)) # bottom, left, top, right
  
  # Note the columns of interest
  columns <- c("Count_BB", "Count_BC", "Count_CB", "Count_CC")
  
  # Calculate the range of the transition rates
  yLim <- range(weightedSampleOfTransitionRates[, columns])

  # Create an empty plot
  plot(x=NULL, y=NULL, bty="n", 
       xlim=c(1, 4), ylim=yLim, las=2, xaxt="n", ylab="Number of transitions", xlab="",
       main="Counting transitions on\n phylogenies", cex.main=2, cex.lab=2, cex.axis=1.25)
  axis(side=1, at=1:4, labels=c("Badger-to-badger", "Badger-to-cow", "Cow-to-badger", "Cow-to-cow"), las=2, cex.axis=1.25)
  
  # Add the bars
  for(i in 1:length(columns)){

    # Calculate the quantiles for the current rate
    quantiles <- quantile(weightedSampleOfTransitionRates[, columns[i]], probs=c(0.025, 0.975))
    
    # Plot a line from the lower to the upper thresholds
    lines(x=c(i, i), y=c(quantiles[1], quantiles[2]))
    points(x=i, y=median(weightedSampleOfTransitionRates[, columns[i]]), pch=19, col="black")
  }
  
  # Add a plot label
  mtext("C", side=3, line=1, at=-0.35, cex=2.5)
  
  # Reset the plotting margins
  par(mar=currentMar)
}

getWeightedSampleOfTransitionCountsBasedOnAICMScores <- function(transitionCountsWithoutBurnIn, aicmScores){
  
  # Calculating the weighted mean rates of transitions between badger and cattle demes
  # across BASTA models
  # Weighting is by the AICM score
  #
  # 1. Each model: Get the posterior distributions for each rate
  # 2. Each model: Calculate sum of posteriors between cattle and badger demes
  #                   sum of cattle -> badger posteriors
  #                   sum of badger -> cattle posteriors
  # 3. Each Model: Calculate AICM
  # 4. Convert model AICM scores to weights:
  #       weight = exp((AIC_min - AIC_i)/2)
  # 5. Normalise weights to sum to 1: 
  #       normalisedWeights = weight / sum(weights)
  # 6. Calculate weighted average rates:
  #      - Sample the posterior sums for cattle->badgers relative to AICM weight
  #      - Sample the posterior sums for badgers->cattle relative to AICM weight
  #       c->b = median(sampleCB, na.rm=TRUE)
  #       b->c = median(sampleBC, na.rm=TRUE)
  #      - NOTE: NAs introduced when flag = 0, these are removed after relative sampling
  #
  # Above method for weighting AICM suggested by Paul Johnson
  # I think this is equivalent to Ensemble Bayesian Model Averaging
  # AICM weights calculating described in: Wagenmakers and Farrell - 2004 - AIC model selection using Akaike weights
  
  # Get the analyses names
  analyses <- names(aicmScores)
  names <- c()
  
  # Get the AICM values
  scores <- c()
  shortenedNames <- c()
  for(i in 1:length(analyses)){
    scores[i] <- aicmScores[[analyses[i]]][1]
    
    parts <- strsplit(analyses[i], split="_")[[1]]
    shortenedNames[i] <- paste(parts[1], parts[2], parts[3], sep="_")
  }
  
  # Convert the AICM scores into weights
  # Exponential penalises those that are bad
  # Flips the score
  modelAICMWeights <- exp((min(scores) - scores)/2)
  
  # Normalise the AICM weights
  normalisedModelAICMWeights <- modelAICMWeights / sum(modelAICMWeights)                           
  
  # Initialise a table to store the weighted samples of the transition counts and associated branch length sums
  weightedSampleOfTransitionRates <- data.frame("Count_BB"=NA, "Count_BC"=NA, "Count_CB"=NA, "Count_CC"=NA,
                                                "SumBranchLengths_BB"=NA, "SumBranchLengths_BC"=NA,
                                                "SumBranchLengths_CB"=NA, "SumBranchLengths_CC"=NA)
  
  # Examine each analysis
  for(i in seq_along(analyses)){
    
    # Get the current analysis name
    name <- shortenedNames[i]
    
    # Get the subset of the transition counts table associated with the current analysis (all replicates)
    analysisTransitionInfo <- transitionCountsWithoutBurnIn[transitionCountsWithoutBurnIn$Analysis == name, ]
    
    # Randomly select row indices based upon the AICM weight associated with the current analysis
    sampledIndices <- sample(1:nrow(analysisTransitionInfo),
                            size=round(normalisedModelAICMWeights[i] * nrow(analysisTransitionInfo), digits=0))
    
    # Store the sample
    weightedSampleOfTransitionRates <- rbind(weightedSampleOfTransitionRates,
                                             analysisTransitionInfo[sampledIndices, 6:13])
  }
  
  # Remove first empty row
  weightedSampleOfTransitionRates <- weightedSampleOfTransitionRates[-1, ]
  
  return(weightedSampleOfTransitionRates)
}

getAICMScores <- function(migrationRateEstimates){
  
  # Initialise a list to store the analysis AICM scores
  scores <- list()
  
  # Examine each analysis
  for(key in names(migrationRateEstimates)){
    
    scores[[key]] <- migrationRateEstimates[[key]]$AICM
  }
  
  return(scores)
}

removeBurnIn <- function(transitionCounts, burnInProp=0.1, nReplicates=3){
  
  # Note the names of the analyses considered
  analyses <- unique(transitionCounts$Analysis)
  
  # Examine each analysis
  for(analysis in analyses){
    
    # Examine each replicate
    for(rep in 1:nReplicates){
      
      # Get the indices assocated with the current analysis
      rows <- which(transitionCounts$Analysis == analysis & transitionCounts$Replicate == rep)
      
      # Note the first burn-in indices
      burnIn <- round(burnInProp * length(rows), digits=0)
      indicesToRemove <- rows[1:burnIn]
      
      # Remove the burn-in rows from the transition counts
      transitionCounts <- transitionCounts[-indicesToRemove, ]
    }
  }
  
  return(transitionCounts)
}

getNSitesInFASTA <- function(fastaFile){
  
  # Open a connection to a file to read (open="r")
  connection <- file(fastaFile, open="r")
  
  # Get first line of file
  firstLine <- readLines(connection, n=1)
  
  # Close file connection
  close(connection)
  
  # Get the number of sites used in the FASTA file from first line
  nSites <- as.numeric(strsplit(firstLine, " ")[[1]][2])
  
  return(nSites)
}

getConstantSiteCounts <- function(constantSiteCountsFile){
 
  # Open the file and store the file lines
  connection <- file(constantSiteCountsFile, "r")
  lines <- readLines(connection)
  close(connection)
  
  # Calculate the constant site counts
  counts <- c(0,0,0,0)
  for(line in lines){
    
    # Skip lines without counts
    if(grepl(line, pattern="Counts|A") == TRUE || line == ""){
      next
    }
    
    # Split the line into its parts
    parts <- strsplit(line, "\t")[[1]]
    
    # Add the counts to the tally
    counts <- counts + as.numeric(parts)
  }
  
  return(counts)
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

plotTransitionRatesBetweenBadgerAndCow <- function(badgerToCow, cowToBadger, code, 
                                                   arrowFactor){
  
  # Create empty plot
  createEmptyPlot()

  # Get deme names and assign colours
  demeNames <- c("badgers", "cattle")
  demeColours <- c("red", "blue")
  
  # Normalise the values between 0 and arrowFactor
  values <- c(badgerToCow, cowToBadger)
  values <- (values / max(values)) * arrowFactor

  # Add labels
  x <- c(0.1, 0.9)
  y <- c(0.5, 0.5)
  text(x=x, y=y, 
       labels=demeNames,
       col=demeColours, cex=2)
  
  # badger -> cow
  arrows(x0=x[1]+0.15, x1=x[2]-0.15, y0=y[1]-0.15, y1=y[2]-0.15,
         code=code, lwd=values[1])
  
  # cow -> badger
  arrows(x0=x[2]-0.15, x1=x[1]+0.15, y0=y[2]+0.15, y1=y[1]+0.15,
         code=code, lwd=values[2])
  
  legend("topleft", legend=c("Badgers -> Cattle",
                       paste("-------------- =", round(badgerToCow/cowToBadger, digits=2)),
                             "Cattle -> Badgers"),
         bty="n", cex=1)
  
  # Reset the margins
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

calculateMeanEstimatedTransitionRatesBetweenCattleAndBadgerPopulationsWeightedByAICM <-
  function(migrationRateEstimates){
    
    # Calculating the weighted mean rates of transitions between badger and cattle demes
    # across BASTA models
    # Weighting is by the AICM score
    #
    # 1. Each model: Get the posterior distributions for each rate
    # 2. Each model: Calculate sum of posteriors between cattle and badger demes
    #                   sum of cattle -> badger posteriors
    #                   sum of badger -> cattle posteriors
    # 3. Each Model: Calculate AICM
    # 4. Convert model AICM scores to weights:
    #       weight = exp((AIC_min - AIC_i)/2)
    # 5. Normalise weights to sum to 1: 
    #       normalisedWeights = weight / sum(weights)
    # 6. Calculate weighted average rates:
    #      - Sample the posterior sums for cattle->badgers relative to AICM weight
    #      - Sample the posterior sums for badgers->cattle relative to AICM weight
    #       c->b = median(sampleCB, na.rm=TRUE)
    #       b->c = median(sampleBC, na.rm=TRUE)
    #      - NOTE: NAs introduced when flag = 0, these are removed after relative sampling
    #
    # Above method for weighting AICM suggested by Paul Johnson
    # I think this is equivalent to Ensemble Bayesian Model Averaging
    # AICM weights calculating described in: Wagenmakers and Farrell - 2004 - AIC model selection using Akaike weights
    
    # Get the analyses names
    analyses <- names(migrationRateEstimates)
    names <- c()
    
    # Get the AICM values
    aicmScores <- c()
    shortenedNames <- c()
    for(i in 1:length(analyses)){
      aicmScores[i] <- migrationRateEstimates[[analyses[i]]][["AICM"]][1]
      
      parts <- strsplit(analyses[i], split="_")[[1]]
      shortenedNames[i] <- paste(parts[1], parts[2], sep="_")
    }
    
    # Convert the AICM scores into weights
    # Exponential penalises those that are bad
    # Flips the score
    modelAICMWeights <- exp((min(aicmScores) - aicmScores)/2)
    
    # Normalise the AICM weights
    normalisedModelAICMWeights <- modelAICMWeights / sum(modelAICMWeights)                           
    
    # Initialise two arrays to store samples of the posterior sums from each model - sample size relative to AICM weight
    badgerToCowRatePosteriorSumSamples<- c()
    cowToBadgerRatePosteriorSumSamples <- c()
    
    # Initialise a vectors to store the median and range of the total rates between cattle and badgers
    modelBadgerToCowRateMedians <- c()
    modelBadgerToCowRateLower <- c()
    modelBadgerToCowRateUpper <- c()
    modelCowToBadgerRateMedians <- c()
    modelCowToBadgerRateLower <- c()
    modelCowToBadgerRateUpper <- c()
    
    # Examine each of the different model structures
    for(i in 1:length(analyses)){
      
      # Initialise arrays to store the rates from badgers to cattle and vice versa
      sumRatesBadgerToCow <- rep(0, length(migrationRateEstimates[[analyses[i]]][[1]]))
      sumRatesCowToBadger <- rep(0, length(migrationRateEstimates[[analyses[i]]][[1]]))
      
      # Get the demeStructure
      parts <- strsplit(analyses[i], "_")[[1]]
      demeStructure <- parts[1]
      names[i] <- paste(parts[c(1,2)], collapse="_")
      
      # Examine each rate
      for(key in names(migrationRateEstimates[[analyses[i]]])){
        
        # Ignore AICM and zero rates (rates never estimate e.g. outer-cattle -> inner-badger)
        if(key == "AICM"){
          next
        }
        
        # Split the key into its deme numbers
        demeNumbers <- as.numeric(strsplit(key, split="_")[[1]])
 
        # Check if current rate is between badger and cattle populations
        if(grepl(getDemeNamesForDemeStructure(demeStructure, demeNumbers[1]),
                 pattern="badger") == TRUE &&
           grepl(getDemeNamesForDemeStructure(demeStructure, demeNumbers[2]),
                 pattern="cow") == TRUE){
          values <- migrationRateEstimates[[analyses[i]]][[key]]
          values[is.na(values)] <- 0
          sumRatesBadgerToCow <- sumRatesBadgerToCow + values
          
        }else if(grepl(getDemeNamesForDemeStructure(demeStructure, demeNumbers[1]),
                       pattern="cow") == TRUE &&
                 grepl(getDemeNamesForDemeStructure(demeStructure, demeNumbers[2]),
                       pattern="badger") == TRUE){
          
          values <- migrationRateEstimates[[analyses[i]]][[key]]
          values[is.na(values)] <- 0
          sumRatesCowToBadger <- sumRatesCowToBadger + values
        }
      }
      
      # Store summary statistics of these posterior sums
      modelBadgerToCowRateMedians[i] <- median(sumRatesBadgerToCow, na.rm=TRUE)
      quantiles <- quantile(sumRatesBadgerToCow, probs=c(0.025, 0.975), na.rm=TRUE)
      modelBadgerToCowRateLower[i] <- quantiles[[1]]
      modelBadgerToCowRateUpper[i] <- quantiles[[2]]
      
      modelCowToBadgerRateMedians[i] <- median(sumRatesCowToBadger, na.rm=TRUE)
      quantiles <- quantile(sumRatesCowToBadger, probs=c(0.025, 0.975), na.rm=TRUE)
      modelCowToBadgerRateLower[i] <- quantiles[[1]]
      modelCowToBadgerRateUpper[i] <- quantiles[[2]]
      
      # Sample from these posterior sums in proportion to the current analyses AIC weight
      badgerToCowRatePosteriorSumSamples <- 
        c(badgerToCowRatePosteriorSumSamples,
          sample(sumRatesBadgerToCow, size=round(normalisedModelAICMWeights[i] * length(sumRatesBadgerToCow), digits=0)))
      cowToBadgerRatePosteriorSumSamples <- 
        c(cowToBadgerRatePosteriorSumSamples,
          sample(sumRatesCowToBadger, size=round(normalisedModelAICMWeights[i] * length(sumRatesCowToBadger), digits=0)))
    }
    
    ## Plot the interspecies transmission rates for each model
    currentMar <- par("mar")
    par(mar=c(17, 5.5, 5, 0.1)) # bottom, left, top, right
    plot(x=NULL, y=NULL, xlim=c(1, length(analyses) + 1), 
         ylim=c(0, max(c(modelBadgerToCowRateUpper, modelCowToBadgerRateUpper), na.rm=TRUE)),
         ylab="", main="Estimated inter-species transition rates",
         las=1, xaxt="n", xlab="", bty="n", cex.lab=2, cex.main=2, cex.axis=1.25)
    mtext(side=2, text="Per lineage transition rate per year", line=3.5, cex=1.25)
    
    for(i in 1:length(analyses)){

      # Badgers to cattle
      points(x=c(i-0.1, i-0.1), y=c(modelBadgerToCowRateLower[i], modelBadgerToCowRateUpper[i]), type="l", col=rgb(1,0,0, 1))
      points(x=i-0.1, y=modelBadgerToCowRateMedians[i], pch=19, col=rgb(1,0,0, 0.75))
      
      # Cattle to Badger
      points(x=c(i+0.1, i+0.1), y=c(modelCowToBadgerRateLower[i], modelCowToBadgerRateUpper[i]), type="l", col=rgb(0,0,1, 1))
      points(x=i+0.1, y=modelCowToBadgerRateMedians[i], pch=19, col=rgb(0,0,1, 0.75))
    }
    
    # Model average - badgers to cattle
    quantiles <- quantile(badgerToCowRatePosteriorSumSamples, probs=c(0.025, 0.975), na.rm=TRUE)
    median <- median(badgerToCowRatePosteriorSumSamples, na.rm=TRUE)
    points(x=c(length(analyses) + 0.9, length(analyses) + 0.9), y=c(quantiles[[1]], quantiles[[2]]), type="l", col="red")
    points(x=length(analyses) + 0.9, y=median, pch=19, col=rgb(1,0,0, 0.75))
    
    # Model average - cattle to badgers
    quantiles <- quantile(cowToBadgerRatePosteriorSumSamples, probs=c(0.025, 0.975), na.rm=TRUE)
    median <- median(cowToBadgerRatePosteriorSumSamples, na.rm=TRUE)
    points(x=c(length(analyses) + 1.1, length(analyses) + 1.1), y=c(quantiles[[1]], quantiles[[2]]), type="l", col="blue")
    points(x=length(analyses) + 1.1, y=median, pch=19, col=rgb(0,0,1, 0.75))
    
    # Axis
    axis(side=1, at=1:(length(analyses)+1), labels=c(names, ""), las=2, cex.axis=1.25)
    axis(side=1, at=length(analyses)+1, labels="Weighted model average", las=2, cex.axis=1.5, col.axis="black", tick=FALSE)
    legend("top", legend=c("Badgers-to-Cattle", "Cattle-to-Badgers"), text.col=c("red", "blue"), bty="n")
    
    # Add a plot label
    mtext("B", side=3, line=1, at=-2, cex=2.5)
    
    # ## Plot badger -> cattle versus cattle -> badger
    # par(mar=c(5.1, 4.1, 0.5, 2.1))
    # max <- max(c(modelBadgerToCowRateMedians, modelCowToBadgerRateMedians), na.rm=TRUE)
    # plot(x=modelBadgerToCowRateMedians, y=modelCowToBadgerRateMedians, las=1,
    #      ylim=c(0, max), xlim=c(0, max), pch=19, bty="n",
    #      col=rgb(1,0,0, 0.75), cex=normalisedModelAICMWeights * 10,
    #      xlab="Badgers -> Cattle lineage transition rate per year",
    #      ylab="Cattle -> Badgers lineage transition rate per year")
    # points(x=c(0, max), y=c(0, max), type="l", lty=2, col="blue")
    # overlayText(x=modelBadgerToCowRateMedians, y=modelCowToBadgerRateMedians,
    #             labels=names, cex=0.5)
    # legend("topleft", legend=c("High", "", "", "Low"), pch=19,
    #        title="AICM Weight", bty="n",
    #        pt.cex=c(2.5, 1.5, 0.5, 0.1))
    # 
    # # Calculate the weighted means for the rate sums between badgers and cattle
    weightedMedianBadgerToCow <- median(badgerToCowRatePosteriorSumSamples, na.rm=TRUE)
    weightedMedianCowToBadger <- median(cowToBadgerRatePosteriorSumSamples, na.rm=TRUE)
    # points(x=weightedMedianBadgerToCow, y=weightedMedianCowToBadger,
    #        pch=15, col=rgb(0,0,0, 0.5), cex=2)
    # text(x=weightedMedianBadgerToCow, y=weightedMedianCowToBadger,
    #      labels="Weighted Median", cex=0.75, pos=2, offset=0.6, col=rgb(0,0,0, 0.75))
    
    par(mar=currentMar)

    return(c(weightedMedianBadgerToCow, weightedMedianCowToBadger))
  }

plotModelAICMScores <- function(migrationRateEstimates, nBootstraps){
  
  # Get the deme structure names
  analyses <- names(migrationRateEstimates)
  
  # Initialise an array to store the shortened analyses names
  names <- c()
  
  # Get the AICM values
  aicmScores <- c()
  bootstrapMins <- c()
  bootstrapMaxs <- c()
  for(i in 1:length(analyses)){
    
    values <- migrationRateEstimates[[analyses[i]]][["AICM"]]
    
    aicmScores[i] <- values[1]
    bootstrapMins[i] <- values[2]
    bootstrapMaxs[i] <- values[3]
    
    parts <- strsplit(analyses[i], split="_")[[1]]
    names[i] <- paste(parts[1], parts[2], sep="_")
  }
  
  # Calculate the range of the AICM scores - y axis limits
  range <- range(aicmScores, bootstrapMins, bootstrapMaxs)
  yLim <- c(range[1] - (0.1 * (range[2] - range[1])),
            range[2] + (0.1 * (range[2] - range[1])))
  
  # Set the margins
  par(mar=c(17, 8, 5, 0)) # bottom, left, top, right
  
  # Set bar colours to highlight min
  colours <- rep("black", length(aicmScores))
  min <- min(aicmScores)
  minIndex <- which(aicmScores == min)
  colours[minIndex] <- "red"
  
  # Create an empty plot
  plot(x=NULL, y=NULL, bty="n", 
       xlim=c(1, length(analyses)), ylim=yLim,
       xaxt="n", yaxt="n", ylab="", xlab="", main="Model AICM score", cex.main=2)
  axis(side=2, at=seq(yLim[1], yLim[2], by=(yLim[2] - yLim[1])/5), las=2, cex.axis=1.25)
  mtext(side=2, text="AICM score", line=6, cex=1.5)
  axis(side=1, at=1:length(analyses),
       labels=names, las=2, cex.axis=1.25)
  
  # Add the bars
  for(i in 1:length(analyses)){
    # polygon(x=c(i-0.4, i-0.4, i+0.4, i+0.4),
    #         y=c(0, aicmScores[i], aicmScores[i], 0),
    #         col=colours[i])
    lines(x=c(i, i),
          y=c(bootstrapMins[i], bootstrapMaxs[i]), col=colours[i])
    points(x=i, y=aicmScores[i], pch=19, col=colours[i])
  }
  
  # Add some horizontal lines
  points(x=c(1, length(analyses)),
         y=c(min, min), lty=2, col="red", type="l")
  points(x=c(1, length(analyses)),
         y=c(bootstrapMins[minIndex], bootstrapMins[minIndex]),
         lty=2, col=rgb(0,0,0, 0.5), type="l")
  points(x=c(1, length(analyses)),
         y=c(bootstrapMaxs[minIndex], bootstrapMaxs[minIndex]),
         lty=2, col=rgb(0,0,0, 0.5), type="l")

  # Add a plot label
  mtext("A", side=3, line=1, at=-2, cex=2.5)
  
  # Reset the margins
  par(mar=c(5.1, 4.1, 4.1, 2.1))
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
  }
  cat("\rFinished examining log tables...\t\t\t\t\t\tn")
  
  # Reset the margins
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  
  return(migrationRateEstimates)
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

getValues <- function(list){
  values <- c()
  keys <- names(list)
  for(i in 1:length(keys)){
    values[i] <- list[[keys[i]]]
  }
  
  return(values)
}

setRangeOfValues <- function(values, minToSet, maxToSet){
  output <- c()
  range <- range(values)
  for(i in 1:length(values)){
    output[i] <- ((values[i] - range[1]) / (range[2] - range[1])) * (maxToSet - minToSet)
    output[i] <- output[i] + minToSet
  }
  
  return(output)
}

timesValuesInListByValue <- function(list, value){
  
  output <- list()
  for(key in names(list)){
    output[[key]] <- list[[key]] * value
  }
  
  return(output)
}

divideValuesInListByMax <- function(arrowRates){
    max <- max(getValues(arrowRates), na.rm=TRUE)
    output <- timesValuesInListByValue(arrowRates, 1/max)
    
  return(output)
}

setRangeOfValuesInList <- function(arrowRates, minToSet, maxToSet){
  
  range <- range(getValues(arrowRates))
  for(key in names(arrowRates)){
    arrowRates[[key]] <- ((arrowRates[[key]] - range[1]) / (range[2] - range[1]) * maxToSet) + minToSet
  }
  
  return(arrowRates)
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

plotMigrationRates <- function(logTable, demeStructure, code, arrowFactor){
  
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
  
  return(output)
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
  if(arrowWeights[["4_5"]] != 0){
    arrows(x0=x[5], x1=x[6], y0=y[5]-0.05, y1=y[6]+0.05,
           code=code, lwd=arrowWeights[["5_6"]])
  }
  # cow-outer-south -> cow-outer-north
  if(arrowWeights[["5_4"]] != 0){
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

createEmptyPlot <- function(){
  par(mar=c(0,0,0,0))
  plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", bty="n",
       ylab="", xlab="")
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

checkDemeSamplingForDemeStructure <- function(demeStructure, number=NULL){
  
  # Species-InOrOut-Location
  # badger  inner   east
  # cow     outer   west
  #                 north
  #                 south
  demeNames <- list(
    # badger cow
    "2Deme"=c(TRUE, TRUE),
    
    # badger-inner cow-inner outer
    "3Deme-outerIsBoth"=c(TRUE, TRUE, TRUE),
    
    # badger-inner cow badger-outer
    "3Deme-outerIsBadger"=c(TRUE, TRUE, FALSE),
    
    # badger cow-inner cow-outer
    "3Deme-outerIsCattle"=c(TRUE, TRUE, TRUE),
    
    # badger-inner cow-inner cow-outer badger-outer
    "4Deme"=c(TRUE, TRUE, TRUE, FALSE),
    
    # badger-inner cow-inner cow-outer-east cow-outer-west badger-outer-east
    # badger-outer-west"),
    "6Deme-EastWest"=c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
    
    # badger-inner cow-inner cow-outer-north cow-outer-south badger-outer-north
    # badger-outer-south"),
    "6Deme-NorthSouth"=c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
    
    # badger-inner-east badger-inner-west cow-inner-east cow-inner-west cow-outer-east
    # cow-outer-west badger-outer-east badger-outer-west"),
    "8Deme-EastWest"=c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
    
    # badger-inner-north badger-inner-south cow-inner-north cow-inner-south 
    # cow-outer-north cow-outer-south badger-outer-north badger-outer-south")
    "8Deme-NorthSouth"=c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
  )
  
  output = demeNames[[demeStructure]]
  if(is.null(number) == FALSE){
    output = demeNames[[demeStructure]][number + 1]
  }
  
  return(output)
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

readInBASTALogTables <- function(date, demeStructures, popEstimationTypes,
                                 clockEstimateTypes, path, burnInProp=0.1,
                                 nReplicates=NULL, ignoreIfFlagged=FALSE){
  
  # Store each of the log tables in a list
  logTables <- list()
  
  for(demeStructure in names(demeStructures)){

    for(popEstimationType in popEstimationTypes){
      
      for(clockEstimationType in clockEstimateTypes){
        
        # Build run defining prefix
        prefix <- paste(demeStructure, "_", popEstimationType, "_",
                        clockEstimationType, "_", date, sep="")
        
        # Initilialise a list to store the log tables for the current run's replicates
        logTablesForReplicates <- list()
        
        # Check if repliocates available
        if(is.null(nReplicates) == FALSE){
          
          # Retrieve the data from each replicate
          for(replicate in 1:nReplicates){
            
            # Print progress information
            cat(paste("\rReading: ", prefix, ".log\tReplicate: ", replicate,
                      "\t\t\t\t\t", sep=""))
            
            # Create file name
            file <- paste(path, "Replicate", replicate, "_", date, "/",
                          prefix, "/", prefix, ".log", sep="")
            
            # Read in the file as table
            logTable <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
            
            # Replace "N" in table with NAs
            logTable[logTable == "N"] <- NA
            
            # Remove the burn-in
            burnIn <- round(burnInProp * nrow(logTable), digits=0)
            logTable <- logTable[burnIn:nrow(logTable), ]
            
            # Calculate the forward rates
            logTable <- calculateForwardMigrationRates(logTable)
            
            # Store the logTable
            logTablesForReplicates[[as.character(replicate)]] <- logTable
          }
          
          # Store the tables as a single log table
          logTables[[prefix]] <- combineLogTables(logTablesForReplicates)
          
        }else{
          
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
          logTable <- calculateForwardMigrationRates(logTable)

          # Store the tables as a single log table
          logTables[[prefix]] <- logTable
        }
      }
    }
  }
  cat("\rFinished reading in log tables...\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n")
  
  return(logTables)
}

combineLogTables <- function(logTablesForReplicates){
  
  output <- logTablesForReplicates[[1]]
  if(length(logTablesForReplicates) > 1){
    for(i in 2:length(logTablesForReplicates)){
      output <- rbind(output, logTablesForReplicates[[i]])
    }
  }
  
  return(output)
}

# Text overlay functions
overlayText <- function(x, y, labels, cex){
  
  ###########################################
  # Produce a list of alternative locations #
  ###########################################
  
  # Get the axis limits
  axisLimits <- par("usr")
  
  # Get the text label heights and lengths
  textHeights <- strheight(labels) * cex
  textWidths <- strwidth(labels) * cex
  
  # Define the spacer for each axis
  spacerX <- 0.01 * (axisLimits[2] - axisLimits[1])
  spacerY <- 0.01 * (axisLimits[4] - axisLimits[3])
  
  # Generate the set of points based upon the spacer
  altX <- c()
  altY <- c()
  for(i in seq(axisLimits[1], axisLimits[2], spacerX)){
    for(j in seq(axisLimits[3], axisLimits[4], spacerY)){
      
      altX[length(altX) + 1] <- i
      altY[length(altY) + 1] <- j
    }
  }
  #points(altX, altY, col=rgb(0,0,0, 0.5), pch=20, xpd=TRUE)
  
  # Remove points that are too close to actual values
  remove <- c()
  for(i in 1:length(altX)){
    
    for(j in 1:length(x)){
      
      if(abs(altX[i] - x[j]) < textWidths[j] &&
         abs(altY[i] - y[j]) < textHeights[j]){
        remove[length(remove) + 1] <- i
        break
      }
    }
  }
  #points(altX[remove], altY[remove], col=rgb(1,1,1), pch=20, xpd=TRUE)
  if(length(remove) > 0){
    altX <- altX[-remove]
    altY <- altY[-remove]
  }
  
  ##############################################################
  # Add labels to plot assigning new locations where necessary #
  ##############################################################
  
  # Plot the point label
  for(i in 1:length(x)){
    
    # Is the current point too close to others?
    if(tooClose(x, y, i, textHeights[i], textWidths[i]) == TRUE && length(altX) != 0){
      
      # Get a new location
      newLocationIndex <- chooseNewLocation(x[i], y[i], altX, altY)
      
      # Add label
      text(x=altX[newLocationIndex], y=altY[newLocationIndex],
           labels=labels[i], xpd=TRUE, cex=cex)
      
      # Add line back to previous location
      points(x=c(altX[newLocationIndex], x[i]),
             y=c(altY[newLocationIndex], y[i]),
             type="l", col=rgb(0,0,0, 0.5))
      
      # Remove new location and any locations too close to it
      output <- removeLocationAndThoseCloseToItFromAlternatives(
        altX, altY, newLocationIndex, textHeights[i], textWidths[i])
      altX <- output[["X"]]
      altY <- output[["Y"]]
      
    }else{
      text(x=x[i], y=y[i],
           labels=labels[i], xpd=TRUE, cex=cex)
    }
  }
}

removeLocationAndThoseCloseToItFromAlternatives <- function(altX, altY, index, textHeight, textWidth){
  remove <- c(index)
  for(i in 1:length(altX)){
    
    if(i == index){
      next
    }
    
    if(abs(altX[index] - altX[i]) < textWidth &&
       abs(altY[index] - altY[i]) < textHeight){
      remove[length(remove) + 1] <- i
    }
  }
  
  altX <- altX[-remove]
  altY <- altY[-remove]
  
  return(list("X" = altX, "Y" = altY))
}

chooseNewLocation <- function(x, y, altXs, altYs){
  
  # Calculate the distance from point to all alternatives
  distances <- c()
  for(i in 1:length(altXs)){
    distances[i] <- euclideanDistance(x, y, altXs[i], altYs[i])
  }
  
  return(which.min(distances))
}

tooClose <- function(x, y, index, textHeight, textWidth){
  
  result <- FALSE
  for(i in 1:length(x)){
    
    if(i == index){
      next
    }else if(abs(x[index] - x[i]) < textWidth &&
             abs(y[index] - y[i]) < textHeight){
      result <- TRUE
      break
    }
  }
  
  return(result) 
}

euclideanDistance <- function(x1, y1, x2, y2){
  return(sqrt(sum((x1 - x2)^2 + (y1 - y2)^2)))
}