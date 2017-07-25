##################
# Load Libraries #
##################
library(ape) # Phylogenetic tree
library(geiger) # Phylogenetic tree tools

# Set directory path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/"

#######################################
# Open the tip to root distances file #
#######################################

# Read in the file
file <- paste(path, "allVCFs-IncludingPoor/vcfFiles/",
              "DistancesToRoot_mlTree_28-10-16.txt", sep="")
tipToRootDistances <- read.table(file, header=TRUE, stringsAsFactors=FALSE)

################################################
# Get sampling dates from sampling information #
################################################

# Read in the badger isolate metadata
fileName <- paste(path, "IsolateData/", "BadgerInfo_08-04-15_LatLongs_XY.csv", sep="")
badgerInfo <- read.table(fileName, header=TRUE, stringsAsFactors=FALSE, sep=",")

# Read in the cattle isolate metadata
file <- paste(path, "IsolateData/", 
              "CattleIsolateInfo_LatLongs_plusID_outbreakSize_Coverage_AddedTB1453-TB1456.csv", sep="")
cattleInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Note the isolate sampling dates
tipToRootDistances <- addSamplingDates(badgerInfo, cattleInfo, tipToRootDistances)

##################################################
# Note which isolates are used in BASTA analyses #
##################################################

# Read in the newick tree
file <- paste(path, "allVCFs-IncludingPoor/vcfFiles/",
              "mlTree_Prox-10_plusRef_rmResequenced_SNPCov-0.1_28-10-16.tree", sep="")
tree <- read.tree(file=file)

# Get a list of the isolates in the clade
node <- 289
cladeTips <- tips(tree, node=node)

###########################
# Root-to-tip Regressions #
###########################

# Open a PDF
# file <- paste(path, "allVCFs-IncludingPoor/vcfFiles/",
#               "DistancesToRoot_mlTree_28-10-16.pdf", sep="")
# pdf(file, width=8)

##### All isolates

# Fit linear model
fitLinearRegression(tipToRootDistances, 
                    "Tip to Root Distances vs Sampling Dates - ALL")

##### BASTA isolates only

# Select the isolates used for the BASTA analysis
bastaTipToRootDistances <- tipToRootDistances[
  tipToRootDistances$IsolateID %in% cladeTips, ]

# Fit linear model
fitLinearRegression(bastaTipToRootDistances, 
                    "Tip to Root Distances vs Sampling Dates - BASTA")

# Remove two outliers
subsetTipToRootDistances <- bastaTipToRootDistances[
  bastaTipToRootDistances$PatristicDistanceToRoot > 0.0105, ]

fitLinearRegression(subsetTipToRootDistances, 
                    "Tip to Root Distances vs Sampling Dates - Removed Outliers")


# dev.off()

cor.test(x=tipToRootDistances$PatristicDistanceToRoot, 
         y=as.numeric(tipToRootDistances$Date))

cor.test(x=bastaTipToRootDistances$PatristicDistanceToRoot, 
         y=as.numeric(bastaTipToRootDistances$Date))

cor.test(x=subsetTipToRootDistances$PatristicDistanceToRoot, 
         y=as.numeric(subsetTipToRootDistances$Date))


#############
# FUNCTIONS #
#############

compareRSquaredToNull <- function(tipToRootDistances, rSquared, n){
  
  # Generate null distribution of R squared values by permuting sampling dates
  nullDistribution <- c()
  for(i in 1:n){
    
    # Permute the dates
    randomDates <- sample(tipToRootDistances$Date)
    
    # Fit linear model
    linearModel <- lm(tipToRootDistances$PatristicDistanceToRoot~randomDates)
    summary <- summary(linearModel)
    
    # Store the R squared value
    nullDistribution[i] <- summary$adj.r.squared
  }
  
  # Plot the null distribution
  hist(nullDistribution, las=1, xlab=expression("Adjusted" ~ R^{2} ~ "Value"),
       main=expression("Null distribution of" ~ R^{2} ~ "values"))
  abline(v=rSquared, lty=2, lwd=2, col="red")
  
}

fitLinearRegression <- function(tipToRootDistances, title){
  
  linearModel <- lm(PatristicDistanceToRoot~Date, data=tipToRootDistances)
  summary <- summary(linearModel)
  
  # Get the confidence intervals
  newData <- data.frame(Date=seq(from=min(tipToRootDistances$Date),
                                 to=max(tipToRootDistances$Date), by=1))
  confidenceIntervals <- as.data.frame(predict(linearModel, interval="confidence",
                                               newdata=newData))
  
  # Plot the model output
  plot(tipToRootDistances$Date, tipToRootDistances$PatristicDistanceToRoot,
       pch=20, las=1, ylab="Root-to-tip Distance", xlab="Sampling Date",
       main=title, mgp=c(3,0.5,0))
  abline(linearModel, col="red")
  lines(x=newData$Date, y=confidenceIntervals$upr, col=rgb(1,0,0, 0.2))
  lines(x=newData$Date, y=confidenceIntervals$lwr, col=rgb(1,0,0, 0.2))
  
  legend("bottomleft", 
         legend=c(paste("R^2 = ", round(summary$adj.r.squared, 2)),
                  paste("p-value = ", round(anova(linearModel)$Pr[[1]], 2))),
         bty='n')
  
  # Get the R Squared value
  rSquared <- summary$adj.r.squared
  
  # Create a null distribution
  compareRSquaredToNull(tipToRootDistances, rSquared, 1000)
}

addSamplingDates <- function(badgerInfo, cattleInfo, tipToRootDistances){
  
  # Initialise list to store each isolates sampling date
  samplingDates <- list()
  
  # Note the format of the dates
  dateFormat <- "%d/%m/%Y"
  
  # Read the badger info
  for(row in 1:nrow(badgerInfo)){
    samplingDates[[badgerInfo[row, "WB_id"]]] <- 
      as.Date(badgerInfo[row, "date"], dateFormat)
  }
  
  # Read the cattle info
  for(row in 1:nrow(cattleInfo)){
    samplingDates[[cattleInfo[row, "StrainId"]]] <- 
      as.Date(cattleInfo[row, "DateCultured"], dateFormat)
  }
  
  # Remove Reference if present
  tipToRootDistances <- tipToRootDistances[grepl(x=tipToRootDistances$IsolateID, 
                                                 pattern="Ref-1997") == FALSE, ]
  
  # Add the dates to the tip to root distances table
  tipToRootDistances$Date <- as.Date(rep(NA, nrow(tipToRootDistances)))
  for(row in 1:nrow(tipToRootDistances)){
    tipToRootDistances[row, "Date"] <- samplingDates[[tipToRootDistances[row, "IsolateID"]]]
  }
  
  return(tipToRootDistances)
}

