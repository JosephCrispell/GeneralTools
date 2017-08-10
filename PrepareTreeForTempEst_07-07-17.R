# Load the ape package
library(ape)
library(geiger) # For the tips function

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_13-07-17/"

################################
# Read in sampling information #
################################

# Read in the badger isolate metadata
fileName <- paste(path, "IsolateData/", "BadgerInfo_08-04-15_LatLongs_XY_Centroids.csv", sep="")
badgerInfo <- read.table(fileName, header=TRUE, stringsAsFactors=FALSE, sep=",")

# Read in the cattle isolate metadata
file <- paste(path, "IsolateData/", 
              "CattleIsolateInfo_LatLongs_plusID_outbreakSize_Coverage_AddedTB1453-TB1456-TB1785.csv", sep="")
cattleInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

#######################################
# Read in the Maximum Likelihood tree #
#######################################

# Read in the newick tree
file <- paste(path, "vcfFiles/", "mlTree_01-08-2017.tree", sep="")
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
file <- paste(path, "vcfFiles/",
              "mlTree_DatedTips_01-08-17.tree", sep="")

write.tree(tree, file = file, append = FALSE,
           digits = 20, tree.names = FALSE)

######################
# Select BASTA clade #
######################

# Set node number
nodeDefiningBastaClade <- 209

pdf(paste(path, "vcfFiles/mlTree_BastaClade_01-08-17.pdf", sep=""))

branchColours <- defineBranchColoursOfClade(tree, nodeDefiningBastaClade,
                                            "black", "lightgrey")

plot.phylo(tree, "fan", edge.color=branchColours, edge.width=3,
           show.tip.label=FALSE)
#nodelabels()

dev.off()

# Get basta clade
bastaClade <- extract.clade(tree, nodeDefiningBastaClade)

# Print out tree
file <- paste(path, "vcfFiles/", "mlTree_BASTAClade_DatedTips_01-08-17.tree", sep="")

write.tree(bastaClade, file = file, append = FALSE,
           digits = 20, tree.names = FALSE)

#############
# FUNCTIONS #
#############

defineBranchColoursOfClade <- function(tree, nodeDefiningClade,
                                       colour, defaultColour){
  branchColours <- rep(defaultColour, dim(tree$edge)[1])
  clade <- tips(tree, node=nodeDefiningClade)
  branchesInClades <- which.edge(tree, clade)
  branchColours[branchesInClades] <- colour
  
  return(branchColours)
}

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

