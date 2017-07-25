# Load the ape package
library(ape)
library(geiger) # For the tips function

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/"

################################
# Read in sampling information #
################################

# Read in the badger isolate metadata
fileName <- paste(path, "IsolateData/", "BadgerInfo_08-04-15_LatLongs_XY.csv", sep="")
badgerInfo <- read.table(fileName, header=TRUE, stringsAsFactors=FALSE, sep=",")

# Read in the cattle isolate metadata
file <- paste(path, "IsolateData/", 
              "CattleIsolateInfo_LatLongs_plusID_outbreakSize_Coverage_AddedTB1453-TB1456.csv", sep="")
cattleInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

#######################################
# Read in the Maximum Likelihood tree #
#######################################

# Read in the newick tree
file <- paste(path, "allVCFs-IncludingPoor/vcfFiles/",
              "mlTree_Prox-10_plusRef_rmResequenced_SNPCov-0.1_28-10-16.tree", sep="")
tree <- read.tree(file=file)

# Get a list of the isolates at the tips
isolates <- tree$tip.label

############################
# Add sampling date to tip #
############################

# Get sampling times
tree$tip.label <- addSamplingTimes(isolates, badgerInfo, cattleInfo)

##################################
# Print out tree with new labels #
##################################

# Create a file
file <- paste(path, "allVCFs-IncludingPoor/vcfFiles/",
              "mlTree_DatedTips_28-10-16.tree", sep="")

write.tree(tree, file = file, append = FALSE,
           digits = 20, tree.names = FALSE)

######################
# Select BASTA clade #
######################

# Get basta clade
node <- 289
bastaClade <- extract.clade(tree, node)

# Print out tree
file <- paste(path, "allVCFs-IncludingPoor/vcfFiles/",
              "mlTree_BASTAclade_DatedTips_28-10-16.tree", sep="")

write.tree(bastaClade, file = file, append = FALSE,
           digits = 20, tree.names = FALSE)

#############
# FUNCTIONS #
#############

addSamplingTimes <- function(tipLabels, badgerInfo, cattleInfo){
  
  # Initialise an array to store the new tip labels
  newTipLabels <- c()
  
  # Note the format of the dates
  dateFormat <- "%d/%m/%Y"
  
  # Examine each tip
  for(i in 1:length(tipLabels)){
    
    # Check if reference
    if(tipLabels[i] == "Ref-1997"){
      newTipLabels[i] <- paste(tipLabels[i], "_1997-01-01", sep="")
      
    # Check if badger or cow isolate
    }else if(grepl(x=tipLabels[i], pattern="WB") == TRUE){
      
      row <- which(badgerInfo$WB_id == tipLabels[i])
      
      newTipLabels[i] <- paste(tipLabels[i], "_", 
                               as.character(as.Date(badgerInfo[row, "date"], dateFormat)),
                                            sep="")
      
    }else{
      
      row <- which(cattleInfo$StrainId == tipLabels[i])
      
      newTipLabels[i] <- paste(tipLabels[i], "_", 
                               as.character(as.Date(cattleInfo[row, "DateCultured"],
                                               dateFormat)),
                                            sep="")
    }
  }
  
  return(newTipLabels)
}

