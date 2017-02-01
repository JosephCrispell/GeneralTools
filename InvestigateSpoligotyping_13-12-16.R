# Load the ape package
library(ape)
library(geiger) # For the tips function

####################
# Read in the Tree #
####################

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/"

# Read in the newick tree
file <- paste(path, "InvestigatingCattleMislabelling/",
              "mlTree_Prox-10_plusRef_rmResequenced_SNPCov-0.1_28-10-16.tree", sep="")
tree <- read.tree(file=file)
tips <- tree$tip.label

# Convert Branch lengths to SNPs
fastaLength <- 9464
tree$edge.length <- tree$edge.length * fastaLength

####################################################
# Get the Spoligotype Information for each Isolate #
####################################################

# Cattle Isolates
file <- paste(path, "IsolateData/",
              "CattleIsolateInfo_LatLongs_plusID_outbreakSize_Coverage_AddedTB1453-TB1456.csv", sep="")
cattleInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Badger Isolates
file <- paste(path, "IsolateData/",
              "BadgerInfo_08-04-15_LatLongs_XY.csv", sep="")
badgerInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Initialise an array to store the tip spoligotypes
tipSpoligotypes <- c()

for(tipIndex in 1:length(tree$tip.label)){
  
  spoligotype <- "NA"
  
  # Cattle
  if(grepl(pattern="TB", x=tree$tip.label[tipIndex]) == TRUE){
    
    # Find index in table
    strainIndex <- which(cattleInfo$StrainId == tree$tip.label[tipIndex])
    spoligotype <- cattleInfo[strainIndex, "Spoligotype"]
    
    # Badgers
  }else if(grepl(pattern="WB", x=tree$tip.label[tipIndex]) == TRUE){
    
    # Find index in table
    strainIndex <- which(badgerInfo$WB_id == tree$tip.label[tipIndex])
    spoligotype <- badgerInfo[strainIndex, "AHVLASpoligo"]
  }

  # Replace blanks, NA, or RETEST with NA
  if(is.na(spoligotype) == TRUE || spoligotype == "" || spoligotype == "RETEST" || spoligotype == "NA"){
    spoligotype <- "NA"
  }
  
  tipSpoligotypes[tipIndex] <- spoligotype
  tree$tip.label[tipIndex] <- paste(tree$tip.label[tipIndex], "_", spoligotype, sep="")
}

##################################
# Get the Isolate Coverage Table #
##################################

file <- paste(path, "allVCFs-IncludingPoor/vcfFiles/",
              "isolateGenomeCoverageSummary_28-10-16.txt", sep="")
coverageInfo <- read.table(file, header=TRUE, stringsAsFactors=FALSE)

coverageInfo$IsolateID <- getIsolateIDFromFileNames(coverageInfo$IsolateID)

#####################################
# Plot the spoligotypes on the tree #
#####################################

# Define tip colours by spoligotype
uniqueSpoligotypes <- unique(tipSpoligotypes)
spoligotypeColours <- list("9"="green",  "15"="cyan", "10"="gold", "17"="blue", "NA"="grey")
tipColours <- defineTipColourBasedOnSpoligotype(tipSpoligotypes=tipSpoligotypes, colours=spoligotypeColours)

# Define tip shape by species
tipShapes <- defineTipShapesForSpecies(tipLabels=tree$tip.label, cow=17, badger=16)

# Define tip size by quality
isolateQuality <- getIsolateQuality(coverageInfo)

tipSizes <- defineTipSizeBySequencingQuality(tips, isolateQuality)

# Open an output file
file <- paste(path, "InvestigatingCattleMislabelling/",
              "mlTree_Spoligotypes_26-01-16.pdf", sep="")
pdf(file)

# Plot the phylogenetic tree
plot.phylo(tree, type="fan",
           edge.color="black", edge.width=3,
           show.node.label=TRUE,
           show.tip.label=TRUE, tip.color=tipColours, cex=0.4, label.offset=5)

# Add node labels
nodelabels(node=1:length(tree$tip.label), 
           cex=tipSizes, pch=tipShapes, col=tipColours)

# Add Legends
addLegendForQuality("bottomright", 0.8)
legend("bottomleft", legend=c("COW", "BADGER"),
       pch=c(17, 16), cex=0.65, bty='n')
legend("topleft", legend=uniqueSpoligotypes, text.col=c("green","cyan","gold","blue","grey"), bty="n", cex=0.8)

points(x=c(80, 76, 96, 105, -114, 75), y=c(105, 106, 87, 82, -62, -100), pch=16,
       cex=2.5, col=rgb(0,0,0, 0.5))



dev.off()

#############
# FUNCTIONS #
#############

addLegendForQuality <- function(position, cex){
  
  sizes <- seq(0.1, 1, 0.1)
  
  legend(position, legend=sizes, col="black", pch=16, bty='n',
         pt.cex=sizes, cex=cex)
}

defineTipSizeBySequencingQuality <- function(tipLabels, isolateQuality){
  
  tipQuality <- c()
  for(i in 1:length(tipLabels)){
    if(tipLabels[i] != "Ref-1997"){
      
      tipQuality[i] <- isolateQuality[[tipLabels[i]]]
    }else{
      tipQuality[i] <- 1
    }
  }
  
  return(tipQuality)
}

getIsolateQuality <- function(table){
  isolateQuality <- list()
  for(i in 1:nrow(table)){
    isolateQuality[[table[i, "IsolateID"]]] <- table[i, "PercentageCoverage"]
  }
  
  return(isolateQuality)
}

getIsolateIDFromFileNames <- function(fileNames){
  isolates <- c()
  for(i in 1:length(fileNames)){
    isolates[i] <- strsplit(fileNames[i], split="_")[[1]][1]
  }
  
  return(isolates)
}

defineTipShapesForSpecies <- function(tipLabels, cow, badger){
  
  shapes <- c()
  for(i in 1:length(tipLabels)){
    
    if(grepl(pattern="TB", x=tipLabels[i]) == TRUE){
      shapes[i] <- cow
    }else if(grepl(pattern="WB", x=tipLabels[i]) == TRUE){
      shapes[i] <- badger
    }else{
      shapes[i] <- cow
    }
  }
  
  return(shapes)
}

defineTipColourBasedOnSpoligotype <- function(tipSpoligotypes, colours){
  tipColours <- c()
  for(i in 1:length(tipSpoligotypes)){
    tipColours[i] <- colours[[tipSpoligotypes[i]]]
  }
  return(tipColours)
}