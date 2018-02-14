############
# Packages #
############

library(ape)
library(geiger) # For the tips function
library(adephylo)

########
# Path #
########

path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_13-07-17/"

################
# Read in tree #
################

file <- paste(path, "vcfFiles/",
              "mlTree_BASTAClade_DatedTips_TempEstRooted_12-02-18.tree", sep="")
tree <- readTree(file)

###################################
# Calculate root-to-tip distances #
###################################

tipToRootDistances <- distRoot(tree, tree$tip.label, method="patristic")

############################
# Get sampling information #
############################

# Read in the badger isolate metadata
fileName <- paste(path, "IsolateData/", "BadgerInfo_08-04-15_LatLongs_XY_Centroids.csv", sep="")
badgerInfo <- read.table(fileName, header=TRUE, stringsAsFactors=FALSE, sep=",")

# Read in the cattle isolate metadata
file <- paste(path, "IsolateData/", 
              "CattleIsolateInfo_LatLongs_plusID_outbreakSize_Coverage_AddedStrainIDs.csv",
              sep="")
cattleInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

#########################
# Fit linear regression #
#########################

# Open a pdf
file <- paste(path, "vcfFiles/",
              "RootedBastaCladeTree_tipToRootDistances_12-02-18.pdf", sep="")
pdf(file)


# Combine the sampling dates with tip-to-root distances
tipInfo <- combineSamplingInformationAndTipToRootDistances(tipToRootDistances,
                                                           badgerInfo,
                                                           cattleInfo)

# Fit linear regression
plot(as.Date(tipInfo$SamplingDate), tipInfo$DistanceToRoot,
     las=1, ylab="Tip-to-root Distance", xlab="Sampling Date",
     main="", mgp=c(3,0.5,0), col=rgb(0,0,0, 0.5), cex=2,
     pch=ifelse(grepl(x=tipInfo$IsolateID, pattern="WB") == TRUE,
                19, 17))

linearModel <- lm(tipInfo$DistanceToRoot ~ as.Date(tipInfo$SamplingDate))
summary <- summary(linearModel)
abline(linearModel, col="red", lwd=2)

legend("bottomright", 
       legend=c(paste("R^2 = ", round(summary$adj.r.squared, 2)),
                paste("p-value = ", round(anova(linearModel)$Pr[[1]], 2))),
       bty='n')

dev.off()

#############
# FUNCTIONS #
#############

combineSamplingInformationAndTipToRootDistances <- function(tipToRootDistances,
                                                            badgerInfo, cattleInfo){
  
  # Note the format of the dates
  dateFormat <- "%d/%m/%Y"
  
  # Get the tip labels
  tipLabels <- names(tipToRootDistances)
  
  # Create a dataframe to store the tip information
  tipInfo <- data.frame(matrix(nrow=length(tipLabels), ncol=3))
  colnames(tipInfo) <- c("IsolateID", "DistanceToRoot", "SamplingDate")
  
  # Fill the data frame
  for(i in 1:length(tipLabels)){
    
    # Get the current label
    isolateID <- getIsolateID(tipLabels[i])
    
    # Store it and the distance to root
    tipInfo[i, "IsolateID"] <- isolateID
    tipInfo[i, "DistanceToRoot"] <- tipToRootDistances[[tipLabels[i]]]
    
    # Check if badger or cow
    if(grepl(x=isolateID, pattern="WB") == TRUE){
      
      # Find the isolate's row
      row <- which(badgerInfo$WB_id == isolateID)
      
      # Store the tip to root distance
      tipInfo[i, "SamplingDate"] <- as.character(
        as.Date(badgerInfo[row, "date"], dateFormat))
      
    }else if(grepl(x=isolateID, pattern="TB") == TRUE){
      
      # Find the isolate's row
      row <- which(cattleInfo$StrainId == isolateID)
      
      # Store the tip to root distance
      tipInfo[i, "SamplingDate"] <- as.character(
        as.Date(cattleInfo[row, "DateCultured"], dateFormat))
    }
  }
  
  return(tipInfo)
}

getIsolateID <- function(tipName){
  return(strsplit(tipName, split="_")[[1]][1])
}

readTree <- function(treeFile){
  
  # Get the first line of the file
  con <- file(treeFile,"r")
  firstLine <- readLines(con,n=1)
  close(con)
  
  # Check if nexus of newick tree
  if(grepl(x=firstLine, pattern="#NEXUS") == TRUE){
    tree <- read.nexus(file=file)
  }else{
    tree <- read.read(file=file)
  }
  
  return(tree)
}